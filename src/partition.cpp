#include "partition.h"

#include <bits/types/time_t.h>

// #include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "app/reachabilityconfig.h"
#include "extern/KaHIP/lib/algorithms/strongly_connected_components.h"
#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/definitions.h"
#include "extern/KaHIP/lib/io/graph_io.h"
#include "extern/KaHIP/lib/tools/chronotimer.h"
#include "extern/KaHIP/lib/tools/random_functions.h"
#include "grail/Grail.h"
#include "grail/GraphUtil.h"
#include "grail/exception_list_incremental_plus.h"
#include "graph.h"
#include "lib/algorithms/bfs.h"
#include "lib/reachability/matrix/full_reach.h"
#include "parse_reachability_parameters.h"
#include "query_generator.h"
#include "reachability/oreach.h"
#include "utils.h"
#define InnerBFS
// #define PATHTREE
#define BORDERBFS
// #define NOEQCLASS // 不使用等价类
// #define MTHREAD
#define INNODE

struct Element {
    int value;
    size_t vecIndex;   // 当前元素来自哪个vector
    size_t elemIndex;  // 当前元素在其vector中的位置

    bool operator>(const Element &other) const { return value > other.value; }
};

std::vector<int> mergeAndUniqueUsingMinHeap(const std::vector<std::vector<int>> &vectors) {
    // 使用lambda函数来定义比较逻辑，构建一个最小堆
    auto comp = [](const Element &a, const Element &b) { return a.value > b.value; };
    std::priority_queue<Element, std::vector<Element>, decltype(comp)> minHeap(comp);

    // 初始化堆，每个vector的第一个元素加入堆
    for (size_t i = 0; i < vectors.size(); ++i) {
        if (!vectors[i].empty()) {
            minHeap.push({vectors[i][0], i, 0});
        }
    }

    std::vector<int> result;
    int last_added = INT_MIN;  // 用于去重，初始化为一个不可能的值

    while (!minHeap.empty()) {
        Element current = minHeap.top();
        minHeap.pop();

        // 只有当当前值与最后加入的值不同，才添加到结果中
        if (result.empty() || current.value != last_added) {
            result.push_back(current.value);
            last_added = current.value;
        }

        // 如果当前vector还有元素，把下一个元素加入堆
        if (current.elemIndex + 1 < vectors[current.vecIndex].size()) {
            minHeap.push({vectors[current.vecIndex][current.elemIndex + 1], current.vecIndex, current.elemIndex + 1});
        }
    }

    return result;
}

void Partitioner::PrintPartitionRatio(ptree::Graph &graph) const {
    if (this->blocks.empty()) {
        WARN("empty blocks");
        return;
    }
    double ce = 0;
    int partitionMaxIdx = -1;
    auto size = this->blocks.size();
    vector<vector<int>> matrix(size, vector<int>(size, 0));
    vector<int> blk_size(size, 0);
    for (auto node : graph.vertices()) {
        blk_size[node.partition]++;
        partitionMaxIdx = max(partitionMaxIdx, node.partition);
        for (auto edge : graph.out_edges(node.id)) {
            if (graph[edge].partition != node.partition) {
                matrix[node.partition][graph[edge].partition] = 1;
                ce += 1;
            }
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (matrix[i][j] == 1 and matrix[j][i] == 1) {
                WARN("{} and {} has bidirectional edge", i, j);
            }
        }
    }

    LOG("edge cut ratio = {}", ce / graph.num_edges());
    for (auto &blk : this->blocks) {
        LOG("block {} size = {}", blk->getBlockID(), blk->getSize());
    }
}

void Partitioner::runPartition(ptree::Graph &graph, const string &algo) {
    CCR::Timer timer;
    timer.ticker();
    LOG("begin Run Partition");
    auto gsize = graph.num_vertices();
    LOG("graph vertex size = {}, edge size = {}, block size = {}", gsize, graph.num_edges(), this->minBlkSize);
    this->blocks = vector<Block *>(this->minBlkSize);
    int maxVertexSize = gsize / this->minBlkSize + 1;
    for (int i = 0; i < minBlkSize; i++) {
        if (algo == "PathTree") {
            this->blocks[i] = new Block();
        } else if (algo == "OReach") {
            this->blocks[i] = new OReachBlock();
        } else if (algo == "Grail") {
            this->blocks[i] = new GrailBlock();
        } else {
            ERROR("未定义 {}", algo);
            exit(0);
        }
        this->blocks[i]->local2global.reserve(maxVertexSize + 1);
        // this->blocks[i].global2local = vector<int>(graph.num_edges() + 1, -1);
    }

    auto cmp = [](ptree::Vertex &first, ptree::Vertex &second) {
        if (first.start == second.start) {
            return first.end > second.end;
        }
        return first.start > second.start;
    };
    auto checkCondition = [&maxVertexSize, this](Block *blk) {
        if (blk->getSize() >= maxVertexSize) {
            DEBUG("partition {} 超过阈值 {}", blk->getSize(), maxVertexSize);
            return true;
        }

        return false;
    };

    priority_queue<ptree::Vertex, std::vector<ptree::Vertex>, decltype(cmp)> que(cmp);
    for (auto &node : graph.vertices()) {
        if (graph.in_degree(node.id) == 0) {
            graph[node.id].visited = true;
            que.push(node);
        }
    }

    int partitionIdx = 0;
    while (!que.empty()) {
        blocks[partitionIdx]->setBlockID(partitionIdx);
        ptree::Vertex v = que.top();
        que.pop();

        graph[v.id].partition = partitionIdx;
        blocks[partitionIdx]->addNode(graph[v.id]);
        for (auto t : graph.out_edges(v.id)) {
            blocks[partitionIdx]->increaseCrossOutedge();
            if (not graph[t].visited) {
                que.push(graph[t]);
                graph[t].visited = true;
            }
        }
        for (auto t : graph.in_edges(v.id)) {
            if (graph[t].partition != partitionIdx) {
               
            } else {
                blocks[partitionIdx]->decreaseCrossOutedge();
            }
        }
        auto checkRes = checkCondition(blocks[partitionIdx]);
        if (checkRes) {
            partitionIdx++;
        }
    }

    double cutedge = 0;
    // 处理边界节点
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (auto t : graph.out_edges(i)) {
            int idx = graph[i].partition;
            if (idx == -1 or graph[t].partition == -1) {
                ERROR("自环");
                exit(1);
            }
            if (graph[i].partition == graph[t].partition) {
                blocks[idx]->addEdge(i, t);
            } else {
                cutedge++;
                // 添加边界节点的连接关系
                borderGraph.addVertex(graph[i]);
                int id1 = borderGraph.getVertexId(graph[i]);
                borderGraph.addVertex(graph[t]);
                int id2 = borderGraph.getVertexId(graph[t]);
                // 检查对同一个Vertex是否一致
                if (this->borderGlobal2Local.find(i) != this->borderGlobal2Local.end()) {
                    if (borderGlobal2Local[i] != id1) {
                        ERROR(
                            "边界节点存在两个对"
                            "应值,global = {}, "
                            "origin local "
                            "= {}, now = {}",
                            i, borderGlobal2Local[i], id1);
                        exit(-1);
                    }
                } else {
                    this->borderGlobal2Local[i] = id1;
                }

                if (this->borderGlobal2Local.find(t) != this->borderGlobal2Local.end()) {
                    if (borderGlobal2Local[t] != id2) {
                        ERROR(
                            "边界节点存在两个对"
                            "应值,global = {}, "
                            "origin local "
                            "= {}, now = {}",
                            t, borderGlobal2Local[t], id1);
                        exit(-1);
                    }
                } else {
                    this->borderGlobal2Local[t] = id2;
                }

                // id1 --> id2
                borderGraph.addEdge(id1, id2);
                blocks[idx]->setBorderType(i, 1);
                blocks[graph[t].partition]->setBorderType(t, 2);
            }
        }
    }

    double sum = 0;
 
    timer.ticker();
    LOG("partition end, time = {}, cut size = {}, cut ratio = {}", timer.get_last_consuming(), cutedge,
        cutedge / graph.num_edges());
}

void Partitioner::runLocalReachability() {
#ifdef MTHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < blocks.size(); i++) {
        blocks[i]->runReachability();
    }
}

void Partitioner::runBorderReachability() {
    LOG("[BEFORE] borderGraph vertex size = {}, edge size = {}", this->borderGraph.num_vertices(),
        this->borderGraph.num_edges());
    // 预先分配好节点
    for (int i = 0; i < blocks.size(); i++) {
        blocks[i]->allocateID(this->borderGraph);
    }

    // 1. 构建块内部的联通关系
#ifdef MTHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < blocks.size(); i++) {
        blocks[i]->buildBorderLink(this->borderGraph, this->borderGlobal2Local);
    }
    LOG("[AFTER] borderGraph vertex size = {}, edge size = {}", borderGraph.num_vertices(), borderGraph.num_edges());
    this->index_size += borderGraph.num_edges();
    this->index_size += borderGraph.num_vertices();
#ifdef PathTreeBorder
    // 2. 用paththee计算边界点的可达性
    CCR::Timer timer;
    auto &g = this->borderGraph;
    int gsize = g.num_vertices();
    vector<int> reverse_topo_sort;
    ptree::GraphUtil::topological_sort(g, reverse_topo_sort);
    cout << "#DAG vertex size:" << g.num_vertices() << "\t#DAG edges size:" << g.num_edges() << endl;

    pt = make_unique<ptree::PathTree>(g, reverse_topo_sort);
    timer.ticker();
    ifstream cfile;
    pt->createLabels(4, cfile, false);
    timer.ticker();
    auto labeling_time = timer.get_last_consuming();
    LOG("Border #construction time:{} (ms)", labeling_time);
#endif
}

void Block::allocateID(ptree::Graph &borderGraph) {
#ifdef INNODE
    int osize = borderClassIn.size();
#else
    int osize = borderClassOut.size();
#endif
    if (osize == 0) {
        return;
    }
    int borderGraphSize = borderGraph.num_vertices();
    int lastOneID = borderGraphSize + osize - 1;
    borderGraph.addVertex(lastOneID);
    for (int i = 0; i < osize; i++) {
        this->sout.emplace(make_pair(i, i + borderGraphSize));
    }
}

void Block::buildBorderLink(ptree::Graph &borderGraph, std::unordered_map<int, int> &borderGlobal2Local) {
    std::unordered_set<int> s;
    int size = localGraph.num_vertices();
#ifdef INNODE
    for (int i = 0; i < size; i++) {
        // 连接该节点与相关的out节点
        if (localGraph[i].border_out) {
#ifdef NOEQCLASS
            int gid = this->local2global[i];
            for (auto &node : this->borderOut) {
                if (this->pt->reach(i, node)) {
                    borderGraph.addEdge(borderGlobal2Local[gid], borderGlobal2Local[this->local2global[node]]);
                }
            }
#else

            auto oidx = localGraph[i].borderInIdx;
            if (oidx < 0) {
                continue;
            }
            if (s.count(oidx) == 0) {
                // 加入集合、新建虚拟节点的连接关系
                s.insert(oidx);
                auto idx = this->sout[oidx];

                // 新增一个虚拟节点
                // auto idx = borderGraph.num_vertices();
                s.insert(oidx);
                // borderGraph.addVertex(idx);
                // cout << idx << "  " << borderGraph.num_vertices() << "\n";
                int gid = this->local2global[i];
                // 和虚拟节点相连
                borderGraph.addEdge(idx, borderGlobal2Local[gid]);
                // 虚拟节点再和out边界点相连
                auto out_nodes = this->borderClassIn[oidx];
                for (auto dst_bid : out_nodes) {
                    borderGraph.addEdge(dst_bid, idx);
                }
            } else {
                // auto idx = s[oidx];
                auto idx = this->sout[oidx];
                int gid = this->local2global[i];
                // 和虚拟节点相连
                borderGraph.addEdge(idx, borderGlobal2Local[gid]);
            }
#endif
        }
    }
#else

    for (int i = 0; i < size; i++) {
        // 连接该节点与相关的out节点
        if (localGraph[i].border_in) {
#ifdef NOEQCLASS
            int gid = this->local2global[i];
            for (auto &node : this->borderOut) {
                if (this->pt->reach(i, node)) {
                    borderGraph.addEdge(borderGlobal2Local[gid], borderGlobal2Local[this->local2global[node]]);
                }
            }
#else

            auto oidx = localGraph[i].borderOutIdx;
            if (oidx < 0) {
                continue;
            }
            if (s.count(oidx) == 0) {
                // 加入集合、新建虚拟节点的连接关系
                s.insert(oidx);
                auto idx = this->sout[oidx];

                // 新增一个虚拟节点
                // auto idx = borderGraph.num_vertices();
                s.insert(oidx);
                // borderGraph.addVertex(idx);
                // cout << idx << "  " << borderGraph.num_vertices() << "\n";
                int gid = this->local2global[i];
                // 和虚拟节点相连
                borderGraph.addEdge(borderGlobal2Local[gid], idx);
                // 虚拟节点再和out边界点相连
                auto out_nodes = this->borderClassOut[oidx];
                for (auto dst_bid : out_nodes) {
                    borderGraph.addEdge(idx, dst_bid);
                }
            } else {
                // auto idx = s[oidx];
                auto idx = this->sout[oidx];
                int gid = this->local2global[i];
                // 和虚拟节点相连
                borderGraph.addEdge(borderGlobal2Local[gid], idx);
            }
            // auto out_nodes = this->borderClassOut[oidx];
            // int gid = this->local2global[i];
            // for (auto dst_bid : out_nodes) {
            //     borderGraph.addEdge(borderGlobal2Local[gid], dst_bid);
            // }
#endif
        }
    }
#endif
}

void Block::runReachability() {
    CCR::Timer timer;
    auto &g = localGraph;
    int gsize = localGraph.num_vertices();
    vector<int> reverse_topo_sort;
    ptree::GraphUtil::topological_sort(g, reverse_topo_sort);
    cout << "#DAG vertex size:" << g.num_vertices() << "\t#DAG edges size:" << g.num_edges() << endl;

    pt = make_unique<ptree::PathTree>(g, reverse_topo_sort);
    timer.ticker();
    ifstream cfile;
    pt->createLabels(4, cfile, false);
    timer.ticker();
    auto labeling_time = timer.get_last_consuming();
    LOG("#construction time:{} (ms)", labeling_time);
}

void GrailBlock::runReachability() {
    bool SKIPSCC = true;
    bool BIDIRECTIONAL = false;
    int LABELINGTYPE = 0;
    bool UseExceptions = true;
    bool UsePositiveCut = false;
    bool POOL = false;
    int POOLSIZE = 100;
    int DIM = 2;
    int query_num = 100000;
    char *filename = NULL;
    char *testfilename = NULL;
    bool debug = false;
    bool LEVEL_FILTER = false;
    bool LEVEL_FILTER_I = false;

    float labeling_time, query_time, query_timepart, exceptionlist_time;
    int alg_type = 1;

    auto &g = localGraph;
    cout << "#vertex size:" << g.num_vertices() << "\t#edges size:" << g.num_edges() << endl;
    int s, t;
    int left = 0;
    int gsize = g.num_vertices();

    bool r;
    struct timeval after_time, before_time, after_timepart, before_timepart;

    int *sccmap;
    if (!SKIPSCC) {
        sccmap = new int[gsize];  // store pair of orignal vertex and corresponding
                                  // vertex in merged graph
        vector<int> reverse_topo_sort;

        // merge strongly connected component
        cout << "merging strongly connected component..." << endl;
        gettimeofday(&before_time, NULL);
        GraphUtil::mergeSCC(g, sccmap, reverse_topo_sort);
        gettimeofday(&after_time, NULL);
        query_time = (after_time.tv_sec - before_time.tv_sec) * 1000.0 +
                     (after_time.tv_usec - before_time.tv_usec) * 1.0 / 1000.0;
        cout << "merging time:" << query_time << " (ms)" << endl;
        cout << "#DAG vertex size:" << g.num_vertices() << "\t#DAG edges size:" << g.num_edges() << endl;
    }

    grail::GraphUtil::topo_leveler(g);
    gettimeofday(&before_time, NULL);

    int dimension;
    grailinst = new grail::Grail(g, DIM, LABELINGTYPE, POOL, POOLSIZE);
    // Grail grail(g,DIM,LABELINGTYPE,POOL,POOLSIZE);

    grailinst->set_level_filter(LEVEL_FILTER);
    gettimeofday(&after_time, NULL);

    labeling_time =
        (after_time.tv_sec - before_time.tv_sec) * 1000.0 + (after_time.tv_usec - before_time.tv_usec) * 1.0 / 1000.0;
    cout << "#construction time:" << labeling_time << " (ms)" << endl;

    if (UseExceptions) {
        gettimeofday(&before_time, NULL);

        el = new grail::ExceptionListIncrementalPlus(g, DIM, 0);
        gettimeofday(&after_time, NULL);
        exceptionlist_time = (after_time.tv_sec - before_time.tv_sec) * 1000.0 +
                             (after_time.tv_usec - before_time.tv_usec) * 1.0 / 1000.0;
        cout << "exceptionlist time:" << exceptionlist_time << " (ms)" << endl;
    }
}

bool GrailBlock::query(int src, int dst) {
    auto s = global2local[src];
    auto d = global2local[dst];
    return this->grailinst->reach(s, d, this->el);
}

void Partitioner::runQuery(ptree::Graph &graph, std::vector<CCR::queryInfo> &vec) {
    auto Merge = [](int a, int b) {
        uint64_t result = ((uint64_t)a << 32) | b;
        return result;
    };
    auto Split = [](uint64_t value, uint32_t &int1, uint32_t &int2) {
        int1 = value >> 32;         // 右移32位得到高32位int
        int2 = value & 0xFFFFFFFF;  // 与操作得到低32位int
    };
    std::unordered_set<int64_t> s;

    for (auto const &info : vec) {
        // 1. get partition
        auto src = info.u;
        auto dst = info.w;
        auto mr = Merge(src, dst);
        if (s.count(mr) == 0) {
            s.insert(mr);
        } else {
            continue;
        }
        if (graph[src].partition == graph[dst].partition) {
            auto res = this->blocks[graph[src].partition]->query(src, dst);
            // DEBUG("{} --> {}, 可达 = {}", src, dst, res);
            continue;
        }
        // 1. 对于src找到所有可达出边
        // 2. 对于dst找到所有可达入边
        // 3. 两两直接查
        std::unordered_set<int> srcCandi;
        std::unordered_set<int> dstCandi;

        // 1.
        auto &tmp = blocks[graph[src].partition]->borderOut;
        for (auto t : tmp) {
            if (blocks[graph[src].partition]->query(src, t)) {
                srcCandi.insert(t);
            }
        }
        // 2.
        auto &tmp2 = blocks[graph[dst].partition]->borderIn;
        for (auto t : tmp2) {
            if (blocks[graph[dst].partition]->query(t, dst)) {
                dstCandi.insert(t);
            }
        }

        // 3.
        bool res = false;
        for (auto s : srcCandi) {
            for (auto d : dstCandi) {
                res = pt->reach(this->borderGlobal2Local[s], this->borderGlobal2Local[d]);
                if (res) {
                    break;
                }
            }
            if (res) {
                break;
            }
        }
        DEBUG("{} --> {}, 可达 = {}", src, dst, res);
    }
}

bool file_exists(const std::string &filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

void queryTransform(BiGraph &bg, ptree::Graph &g, int u, int w, time_t start, time_t end,
                    pair<vector<int>, vector<int>> &querypairs) {
    vector<int> v1s;
    vector<int> v2s;

    for (int i = 0; i < bg.timeSection_u[u].size(); i++)
        if (bg.timeSection_u[u][i].first <= start && bg.timeSection_u[u][i].second > start)
            v1s.push_back(bg.adj_matrix_u[u][i]);

    if (v1s.empty()) {
        time_t max_end = 0;
        int v1 = -1;
        for (int i = 0; i < bg.timeSection_u[u].size(); i++)
            if (bg.timeSection_u[u][i].first > start) {
                max_end = bg.timeSection_u[u][i].second;
                v1 = i;
                break;
            }
        if (v1 >= 0) {
            for (int i = v1; i < bg.timeSection_u[u].size(); i++)
                if (bg.timeSection_u[u][i].first < max_end) v1s.push_back(bg.adj_matrix_u[u][i]);
        }
    }

    // cout<<bg.timeSection_u[u].size()<<"\n";
    for (int i = 0; i < bg.timeSection_u[u].size(); i++)
        if (bg.timeSection_u[u][i].first < end && bg.timeSection_u[u][i].second >= end) {
            v2s.push_back(bg.adj_matrix_u[u][i]);
        }

    if (v2s.empty()) {
        time_t min_start = numeric_limits<time_t>::max();
        int v2 = -1;
        for (int i = bg.timeSection_u[u].size() - 1; i >= 0; i--)
            if (bg.timeSection_u[u][i].second < end) {
                min_start = bg.timeSection_u[u][i].first;
                v2 = i;
                break;
            }
        if (v2 >= 0) {
            for (int i = v2; i >= 0; i--)
                if (bg.timeSection_u[u][i].second > min_start) {
                    v2s.push_back(bg.adj_matrix_u[u][i]);
                }
        }
    }

    for (auto v1 : v1s) {
        int v1_idx = -1;
        time_t mm = std::numeric_limits<time_t>::max();

        for (int i = 0; i < bg.timeTable_u[v1].size(); i++)
            if (bg.timeTable_u[v1][i].second > start)
                if (bg.timeTable_u[v1][i].second < mm) {
                    mm = bg.timeTable_u[v1][i].second;
                    v1_idx = i;
                }
        if (v1_idx >= 0) {
            int src = g.cnum(make_tuple(v1, bg.timeTable_u[v1][v1_idx].first, bg.timeTable_u[v1][v1_idx].second));
            querypairs.first.push_back(src);
        }
    }

    for (auto v2 : v2s) {
        int v2_idx = -1;
        time_t mm = 0;
        for (int i = 0; i < bg.timeTable_u[v2].size(); i++)
            if (bg.timeTable_u[v2][i].first < end)
                if (bg.timeTable_u[v2][i].first > mm) {
                    mm = bg.timeTable_u[v2][i].first;
                    v2_idx = i;
                }
        if (v2_idx >= 0) {
            int dst = g.cnum(make_tuple(v2, bg.timeTable_u[v2][v2_idx].first, bg.timeTable_u[v2][v2_idx].second));
            querypairs.second.push_back(v2_idx);
        }
    }
}

void Partitioner::runQueryWithBfs(BiGraph &bg, ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo,
                                  std::vector<int> &queryRes) {
    // Bi BFS
    LOG("query start");
    std::unordered_set<int64_t> s;

    this->que1 = vector<int>(this->borderGraph.num_vertices());
    this->que2 = vector<int>(this->borderGraph.num_vertices());
    this->visited_inner = vector<int>(borderGraph.num_vertices());
    this->visited_outter = vector<int>(borderGraph.num_vertices());
    int size = 0;

    for (auto &t : this->blocks) {
        size = max(size, t->getSize());
    }
    this->inner_que1 = vector<int>(size);
    this->inner_que2 = vector<int>(size);
    int idx = 0;
    for (auto const &info : queryInfo) {
        idx++;
        auto u = info.u;
        auto w = info.w;
        auto start = info.start;
        auto end = info.end;
        // 开始时间晚于结束时间，查询结果为假
        if (start > end) {
            queryRes[idx - 1] = false;
            continue;
        }
        pair<vector<int>, vector<int>> querypairs;
        queryTransform(bg, graph, u, w, start, end, querypairs);
        // 如果不存在v1或v2，则查询结果为false
        auto &srcs = querypairs.first;
        auto &dsts = querypairs.second;

        // LOG("srcs size = {}, dsts = {}",srcs.size(),dsts.size());
        if (querypairs.first.empty() || querypairs.second.empty()) {
            queryRes[idx - 1] = false;
            continue;
        }
        for (auto &src : srcs) {
            if (queryRes[idx - 1]) {
                break;
            }
            for (auto &dst : dsts) {
                if (graph[src].partition == graph[dst].partition) {
                    if (this->blocks[graph[src].partition]->query(src, dst)) {
                        queryRes[idx - 1] = true;
                        break;
                    }
                }
            }
        }

        if (queryRes[idx - 1]) {
            continue;
        }

        int q1lo = 0, q1hi = 0, q2lo = 0, q2hi = 0;
#ifdef InnerBFS
        bool res = getInBorder(q2hi, dsts, graph, idx);
        res |= getOutBorder(q1hi, srcs, graph, idx);
        if (res) {
            queryRes[idx - 1] = res;
            continue;
        }
        // LOG("");

#else
        auto &tmp = blocks[graph[src].partition].borderOut;

        for (auto t : tmp) {
            if (blocks[graph[src].partition].query(src, t)) {
                // srcCandi.insert(this->borderGlobal2Local[t]);
                this->que1[q1hi++] = this->borderGlobal2Local[t];
            }
        }
        // 2.
        auto &tmp2 = blocks[graph[dst].partition].borderIn;
        for (auto t : tmp2) {
            if (blocks[graph[dst].partition].query(t, dst)) {
                // dstCandi.insert(this->borderGlobal2Local[t]);
                this->que2[q2hi++] = this->borderGlobal2Local[t];
            }
        }
#endif
        int sum = 0;
        queryRes[idx - 1] = runBiBFS(q1hi, q2hi, sum, idx);
    }
}

void Partitioner::runQueryWithOReach(ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo) {
    auto Merge = [](int a, int b) {
        uint64_t result = ((uint64_t)a << 32) | b;
        return result;
    };
    auto Split = [](uint64_t value, uint32_t &int1, uint32_t &int2) {
        int1 = value >> 32;         // 右移32位得到高32位int
        int2 = value & 0xFFFFFFFF;  // 与操作得到低32位int
    };
    bool runSReach{true};
    ReachabilityConfig config;
    graph_access input_graph;
    std::string graph_filename{""};
    ChronoTimer t;
    t.restart();

    // if (graph_io::readDiGraphFromPT(input_graph, graph)) {
    //     return;
    // }

    if (graph_io::readDiGraphFromPT(input_graph, this->borderGraph)) {
        return;
    }
    std::cout << "input IO took " << t.elapsed<std::milli>() << " ms\n"
              << "n(input): " << input_graph.number_of_nodes() << "\n"
              << "m(input): " << input_graph.number_of_edges() << std::endl;

    random_functions::setSeed(config.seed);
    if (config.writeCondensationGraph) {
        graph_access condensation;
        strongly_connected_components scc;
        std::vector<int> component_assignment(input_graph.number_of_nodes());
        scc.build_scc_graph(input_graph, component_assignment, condensation);
        graph_io::writeGraph(condensation, config.outputCondensation, false);
        std::cout << "done writing condensation to " << config.outputCondensation << "\n";
        graph_filename = config.outputCondensation;
    }
    double time_scc_init = 0;

    oreach r(input_graph, config);

    t.restart();
    r.initialize();
    time_scc_init = t.elapsed<std::milli>();
    std::cout << "init took " << time_scc_init << " ms\n";
    FullReach fullreach(input_graph);
    std::vector<std::vector<Query>> all_queries;
    std::vector<std::string> query_types;
    size_t N = queryInfo.size();

    // QueryGenerator generator(config, input_graph, 100000, r);
    // generator.setWriteOutput(false);
    // generator.generate_random();
    // queries = generator.getQueries();
    std::unordered_set<int64_t> s;
    std::vector<std::pair<int, int>> nodes;

    for (size_t i = 0; i < queryInfo.size(); i++) {
        auto src = queryInfo[i].u;
        auto dst = queryInfo[i].w;
        // auto mr = Merge(src, dst);
        // if (s.count(mr) == 0) {
        //     s.insert(mr);
        // } else {
        //     continue;
        // }

        if (graph[src].partition == graph[dst].partition) {
            auto res = this->blocks[graph[src].partition]->query(src, dst);
            DEBUG("{} --> {}, 可达 = {}", src, dst, res);
            continue;
        }
        nodes.emplace_back(std::make_pair(src, dst));
        std::unordered_set<int> srcCandi;
        std::unordered_set<int> dstCandi;
        CCR::Timer t1;
        t1.start();
        auto &tmp = blocks[graph[src].partition]->borderOut;
        for (auto t : tmp) {
            if (blocks[graph[src].partition]->query(src, t)) {
                srcCandi.insert(t);
            }
        }
        // 2.
        auto &tmp2 = blocks[graph[dst].partition]->borderIn;
        for (auto t : tmp2) {
            if (blocks[graph[dst].partition]->query(t, dst)) {
                dstCandi.insert(t);
            }
        }

        std::vector<Query> queries;
        queries.reserve(srcCandi.size() * dstCandi.size());
        if (srcCandi.size() == 0 || dstCandi.size() == 0) {
            continue;
        }
        for (auto s : srcCandi) {
            for (auto d : dstCandi) {
                // res = pt->reach(this->borderGlobal2Local[s],
                //                 this->borderGlobal2Local[d]);
                Query q(this->borderGlobal2Local[s], this->borderGlobal2Local[d]);
                queries.push_back(q);
            }
        }

        all_queries.push_back(queries);
        t1.ticker();
        LOG("src query = {}, dst query = {}, total query = {}, generate time = "
            "{}",
            srcCandi.size(), dstCandi.size(), queries.size(), t1.get_last_consuming());
    }

    query_types.push_back("random");
    std::vector<FallbackMode> fallbacks;

    fallbacks.push_back(FB_NONE);

    for (std::size_t q = 0; q < all_queries.size(); q++) {
        double time_scc_query_pbibfs = 0;
        double time_scc_query_pdfs = 0;
        double time_cur_query = 0;
        double time_fullreach_query = 0;
        auto num_queries = all_queries[q].size();

        // std::cout << "\n======= " << query_types[q] << " ========\n";
        int successful_queries = 0;
        bool res = false;
        if (runSReach) {
            // double elapsed = 0;
            for (auto fb : fallbacks) {
                r.setFallbackMode(fb);
                std::string algorithm;
                std::vector<Query> queries(all_queries[q]);
                switch (fb) {
                    case FB_PIBFS:
                        algorithm = "oreach-pbibfs";
                        std::cout << "\n\nRunning "
                                     "O'Reach with "
                                     "fallback "
                                     "pBiBFS ..."
                                  << std::endl;
                        t.restart();
                        // r.query<FB_PIBFS>(all_queries[q]);
                        r.query<FB_PIBFS>(queries);
                        time_scc_query_pbibfs = t.elapsed<std::micro>();
                        time_cur_query = time_scc_query_pbibfs;
                        break;
                    case FB_PDFS:
                        algorithm = "oreach-pdfs";
                        std::cout << "\n\nRunning "
                                     "O'Reach with "
                                     "fallback pDFS "
                                     "..."
                                  << std::endl;
                        t.restart();
                        // r.query<FB_PDFS>(all_queries[q]);
                        r.query<FB_PDFS>(queries);
                        time_scc_query_pdfs = t.elapsed<std::micro>();
                        time_cur_query = time_scc_query_pdfs;
                        break;
                    case FB_NONE:
                    default:
                        algorithm = "none";
                        t.restart();
                        res = r.query<FB_NONE>(queries);
                        time_cur_query = t.elapsed<std::micro>();

                        break;
                }

              
                r.reset();
            }
        }
        // TODO: 去掉重复query
        DEBUG("{} --> {}, 可达 = {}", nodes[q].first, nodes[q].second, res);
    }
}

bool Partitioner::getInBorder(int &q2hi, vector<int> dsts, ptree::Graph &graph, int idx) {
#ifdef NOEQCLASS
    for (auto &dst : dsts) {
        auto src_partition = graph[dst].partition;
        auto &tmp2 = blocks[src_partition]->borderIn;
        for (auto t : tmp2) {
            if (blocks[graph[dst].partition]->query(t, dst)) {
                // dstCandi.insert(this->borderGlobal2Local[t]);
                this->que2[q2hi++] = this->borderGlobal2Local[t];
            }
        }
    }
    return false;
#else
    for (auto &dst : dsts) {
        auto src_partition = graph[dst].partition;
        bool res = this->blocks[src_partition]->getInBorder2(dst, this->que2, this->visited_outter, idx, q2hi);
        if (res) {
            return true;
        }
    }
    return false;
#endif
}

bool Partitioner::getOutBorder(int &q1hi, vector<int> srcs, ptree::Graph &graph, int idx) {
#ifdef NOEQCLASS
    for (auto &src : srcs) {
        auto src_partition = graph[src].partition;
        auto &tmp2 = blocks[src_partition]->borderOut;
        for (auto t : tmp2) {
            if (blocks[graph[src].partition]->query(src, t)) {
                // dstCandi.insert(this->borderGlobal2Local[t]);
                this->que1[q1hi++] = this->borderGlobal2Local[t];
            }
        }
    }
    return false;
#else
    for (auto &src : srcs) {
        auto src_partition = graph[src].partition;
        bool res = this->blocks[src_partition]->getOutBorder2(src, this->que1, this->visited_outter, idx, q1hi);
        if (res) {
            return true;
        }
    }
    return false;
#endif
}

void Partitioner::getInBorder(int &q2hi, int dst, ptree::Graph &graph, int idx) {
    auto src_partition = graph[dst].partition;
    q2hi = this->blocks[src_partition]->getInBorder(dst, this->que2);
}

void Partitioner::getOutBorder(int &q1hi, int src, ptree::Graph &graph, int idx) {
    auto src_partition = graph[src].partition;
    q1hi = this->blocks[src_partition]->getOutBorder(src, this->que1);
}

bool Partitioner::runBiBFS(int q1hi, int q2hi, int &sum, int idx) {
    auto &g = this->borderGraph;
    int q1lo = 0;
    int q2lo = 0;
    sum = 0;
    while (q1lo < q1hi and q2lo < q2hi) {
        sum++;
        if (q1hi - q1lo < q2hi - q2lo) {
            auto origin = q1hi;
            for (int i = q1lo; i < origin; i++) {
                auto cur = que1[q1lo++];
                if (visited_outter[cur] == -idx) {
                    return true;
                }
                visited_outter[cur] = idx;
                for (auto edge : g.out_edges(cur)) {
                    if (visited_outter[edge] == -idx) {
                        return true;
                    } else if (visited_outter[edge] != idx) {
                        visited_outter[edge] = idx;
                        this->que1[q1hi++] = edge;
                    }
                }
            }
        } else {
            auto origin = q2hi;
            for (int i = q2lo; i < origin; i++) {
                auto cur = que2[q2lo++];
                if (visited_outter[cur] == idx) {
                    return true;
                }
                visited_outter[cur] = -idx;
                for (auto edge : g.in_edges(cur)) {
                    if (visited_outter[edge] == idx) {
                        return true;
                    } else if (visited_outter[edge] != -idx) {
                        visited_outter[edge] = -idx;
                        this->que2[q2hi++] = edge;
                    }
                }
            }
        }
    }
    return false;
}

void Partitioner::ComputeBorder() {
#ifdef NOEQCLASS
    return;
#endif

#ifdef MTHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < blocks.size(); i++) {
#ifdef PATHTREE
        blocks[i].fetchClassByPT(this->borderGlobal2Local);
#elif defined BORDERBFS
        blocks[i]->fetchClassByBorderBFS(this->borderGlobal2Local);
#else
        blocks[i].fetchClass(this->borderGlobal2Local);
#endif
    }
}

void Block::fetchClass(std::unordered_map<int, int> &borderGlobal2Local) {
    std::unordered_map<size_t, int> table;
    int size = localGraph.num_vertices();
    LOG("total node = {}", size);
    for (int i = 0; i < size; i++) {
        std::deque<int> deq;
        std::set<int> visited;
        deq.push_back(i);
        std::unordered_set<int> s;
        while (!deq.empty()) {
            auto node = deq.front();
            deq.pop_front();
            if (localGraph[node].border_out) {
                s.insert(borderGlobal2Local[local2global[node]]);
            }
            for (auto edge : localGraph.out_edges(node)) {
                if (visited.count(edge) == 0) {
                    deq.push_back(edge);
                    visited.insert(edge);
                }
            }
        }
        if (s.empty()) {
            continue;
        }
        std::vector<int> candidate(s.begin(), s.end());
        std::sort(candidate.begin(), candidate.end());
        auto hash_res = hash_vector(candidate);
        if (table.find(hash_res) == table.end()) {
            int idx = this->borderClassOut.size();
            this->borderClassOut.emplace_back(candidate);
            table[hash_res] = idx;
        }
        localGraph[i].borderOutIdx = table[hash_res];

        if (i % 10000 == 0) {
            LOG("out block id = {}, 遍历了{}个节点, 新图节点 = {}", this->blockID, i, table.size());
        }
    }
    LOG("after compact, out node = {}, table size = {}", this->borderClassOut.size(), table.size());

    table = std::unordered_map<size_t, int>();
    for (int i = 0; i < size; i++) {
        std::deque<int> deq;
        std::set<int> visited;
        deq.push_back(i);
        std::unordered_set<int> s;
        while (!deq.empty()) {
            auto node = deq.front();
            deq.pop_front();
            if (localGraph[node].border_in) {
                s.insert(borderGlobal2Local[local2global[node]]);
            }
            for (auto edge : localGraph.in_edges(node)) {
                if (visited.count(edge) == 0) {
                    deq.push_back(edge);
                    visited.insert(edge);
                }
            }
        }
        if (s.empty()) {
            continue;
        }
        std::vector<int> candidate(s.begin(), s.end());
        std::sort(candidate.begin(), candidate.end());
        auto hash_res = hash_vector(candidate);
        if (table.find(hash_res) == table.end()) {
            int idx = this->borderClassIn.size();
            this->borderClassIn.emplace_back(candidate);
            table[hash_res] = idx;
        }
        localGraph[i].borderInIdx = table[hash_res];
        // table[ss.str()]++;
        if (i % 10000 == 0) {
            LOG("in block id = {}, 遍历了{}个节点, 新图节点 = {}", this->blockID, i, table.size());
        }
    }
    LOG("after compact, in node = {}", this->borderClassIn.size());
}

void Block::fetchClassByPT(std::unordered_map<int, int> &borderGlobal2Local) {
    std::unordered_map<size_t, int> table;
    int size = localGraph.num_vertices();
    LOG("total node = {}", size);
    for (int i = 0; i < size; i++) {
        std::unordered_set<int> s;
        for (auto border_node : this->borderOut) {
            if (pt->reach(i, border_node)) {
                s.insert(borderGlobal2Local[local2global[border_node]]);
            }
        }
        if (s.empty()) {
            continue;
        }
        std::vector<int> candidate(s.begin(), s.end());
        std::sort(candidate.begin(), candidate.end());
        auto hash_res = hash_vector(candidate);
        if (table.find(hash_res) == table.end()) {
            int idx = this->borderClassOut.size();
            this->borderClassOut.emplace_back(candidate);
            table[hash_res] = idx;
        }
        localGraph[i].borderOutIdx = table[hash_res];
        if (i % 10000 == 0) {
            LOG("out block id = {}, 遍历了{}个节点, 新图节点 = {}", this->blockID, i, table.size());
        }
    }
    LOG("after compact, out node = {}, table size = {}", this->borderClassOut.size(), table.size());

    table = std::unordered_map<size_t, int>();
    for (int i = 0; i < size; i++) {
        std::unordered_set<int> s;
        for (auto border_node : this->borderIn) {
            if (pt->reach(border_node, i)) {
                s.insert(borderGlobal2Local[local2global[border_node]]);
            }
        }
        if (s.empty()) {
            continue;
        }
        std::vector<int> candidate(s.begin(), s.end());
        std::sort(candidate.begin(), candidate.end());
        auto hash_res = hash_vector(candidate);
        if (table.find(hash_res) == table.end()) {
            int idx = this->borderClassIn.size();
            this->borderClassIn.emplace_back(candidate);
            table[hash_res] = idx;
        }
        localGraph[i].borderInIdx = table[hash_res];
        if (i % 10000 == 0) {
            LOG("in block id = {}, 遍历了{}个节点, 新图节点 = {}", this->blockID, i, table.size());
        }
    }
    LOG("after compact, in node = {}, table size = {}", this->borderClassIn.size(), table.size());
}

void Block::fetchClassByBorderBFS(std::unordered_map<int, int> &borderGlobal2Local) {
    // 1 反向遍历graph,记录每个节点的边界出度
    std::unordered_map<size_t, int> table;
    int size = localGraph.num_vertices();
    LOG("fetchClassByBorderBFS, total node = {}", size);
    vector<int> border_degree(size, 0);
    vector<int> que(size, 0);
    vector<int> visited(size, 0);
    vector<vector<int>> buffer(size);

    int lo = 0;
    int hi = 0;
    for (auto border_node : this->borderOut) {
        que[hi++] = border_node;
        visited[border_node] = 1;
    }
    // 1.2 收集one hop有哪些节点
    // 1.3 遍历整图收集边界degree
    while (lo < hi) {
        auto node = que[lo++];
        for (auto inEdge : localGraph.in_edges(node)) {
            border_degree[inEdge]++;
            if (visited[inEdge] != 1) {
                visited[inEdge] = 1;
                que[hi++] = inEdge;
            }
        }
    }
    if (not this->borderOut.empty()) {
        for (int i = 0; i < size; i++) {
            if (border_degree[i] > 0) {
                buffer[i].reserve(border_degree[i]);
            }
        }
    }

    // 3. 反向bfs
    // 3.1 先把source找好
    lo = 0;
    hi = 0;
    for (auto border_node : this->borderOut) {
        if (border_degree[border_node] == 0) {
            que[hi++] = border_node;
        }
    }
    while (lo < hi) {
        auto node = que[lo++];
        // 收集自己出边传来的class信息
        std::vector<std::vector<int>> allCandidate;
        allCandidate.reserve(1 + border_degree[node]);
        if (localGraph[node].border_out) {
            allCandidate.push_back({borderGlobal2Local[local2global[node]]});
        }
        // 归并来自子节点的

        for (auto &kid : buffer[node]) {
            allCandidate.push_back(this->borderClassOut[kid]);
        }
        // 归并且去重
        auto res = mergeAndUniqueUsingMinHeap(allCandidate);
        auto hash_res = hash_vector(res);
        if (table.find(hash_res) == table.end()) {
            this->index_size += res.size();
            this->index_size += res.size();
            int idx = this->borderClassOut.size();
            this->borderClassOut.emplace_back(res);
            table[hash_res] = idx;
        }
        localGraph[node].borderOutIdx = table[hash_res];

        // 向前传递信息
        for (auto inEdge : localGraph.in_edges(node)) {
            border_degree[inEdge]--;
            buffer[inEdge].push_back(table[hash_res]);

            if (border_degree[inEdge] == 0) {
                que[hi++] = inEdge;
            }
        }
    }
    this->eqClassVisitedO = vector<int>(table.size());
    LOG("after compact, out node = {}, table size = {}", this->borderClassOut.size(), table.size());

    // 从入边界点进行处理
    table = std::unordered_map<size_t, int>();
    lo = 0;
    hi = 0;
    for (auto border_node : this->borderIn) {
        que[hi++] = border_node;
        visited[border_node] = -1;
    }
    while (lo < hi) {
        auto node = que[lo++];
        for (auto inEdge : localGraph.out_edges(node)) {
            border_degree[inEdge]++;
            if (visited[inEdge] != -1) {
                visited[inEdge] = -1;
                que[hi++] = inEdge;
            }
        }
    }
    if (not this->borderIn.empty()) {
        for (int i = 0; i < size; i++) {
            buffer[i].clear();
            if (border_degree[i] > 0) {
                buffer[i].reserve(border_degree[i]);
            }
        }
    }

    lo = 0;
    hi = 0;
    for (auto border_node : this->borderIn) {
        if (border_degree[border_node] == 0) {
            que[hi++] = border_node;
        }
    }

    while (lo < hi) {
        auto node = que[lo++];
        std::vector<std::vector<int>> allCandidate;
        allCandidate.reserve(1 + border_degree[node]);
        if (localGraph[node].border_in) {
            allCandidate.push_back({borderGlobal2Local[local2global[node]]});
        }
        // 归并来自子节点的

        for (auto &kid : buffer[node]) {
            allCandidate.push_back(this->borderClassIn[kid]);
        }
        // 归并且去重
        auto res = mergeAndUniqueUsingMinHeap(allCandidate);
        auto hash_res = hash_vector(res);
        // LOG("res size = {}",res.size());
        if (table.find(hash_res) == table.end()) {
            this->index_size += res.size();
            this->index_size += res.size();
            int idx = this->borderClassIn.size();
            this->borderClassIn.emplace_back(res);
            table[hash_res] = idx;
        }
        localGraph[node].borderInIdx = table[hash_res];

        // 向前传递信息
        for (auto outEdge : localGraph.out_edges(node)) {
            border_degree[outEdge]--;
            buffer[outEdge].push_back(table[hash_res]);

            if (border_degree[outEdge] == 0) {
                que[hi++] = outEdge;
            }
        }
    }
    this->eqClassVisitedI = vector<int>(table.size());
    LOG("after compact, in node = {}, table size = {}", this->borderClassIn.size(), table.size());
}

void OReachBlock::runReachability() {
    LOG("use OReach");
    // 内部构建
    ReachabilityConfig config;
    graph_access input_graph;
    std::string graph_filename{""};
    ChronoTimer t;
    t.restart();
    if (graph_io::readDiGraphFromPT(input_graph, this->localGraph)) {
        return;
    }
    double time_scc_init = 0;
    this->r = new oreach(input_graph, config);
    t.restart();
    r->initialize();
    time_scc_init = t.elapsed<std::milli>();
    std::cout << "init took " << time_scc_init << " ms\n";
}

bool OReachBlock::query(int src, int dst) {
    auto s = global2local[src];
    auto d = global2local[dst];
    Query query(s, d);
    r->singleQuery<FB_NONE>(query);
    query.setBinaryFromResult();
    return query.reachable;
}
