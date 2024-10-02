#include "utils.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <utility>

#include "CCRTree.h"

namespace CCR {
vector<string> split(string line, const char *symbols, int lenth) {
    vector<string> res = vector<string>();
    char *m;
    char *slice = NULL;
    char *buffer_ = (char *)malloc(sizeof(char) * lenth);
    strcpy(buffer_, line.c_str());
    m = strtok_r(buffer_, symbols, &slice);
    while (m != nullptr) {
        res.push_back(m);
        m = strtok_r(nullptr, symbols, &slice);
    }
    return res;
}

Graph transformation(BiGraph &bg, int delta) {
    Graph g;
    vector<vector<pair<int, tuple<int, time_t, time_t>>>> flow(bg.adj_matrix_u.size(),
                                                               vector<pair<int, tuple<int, time_t, time_t>>>());
    for (int v = 0; v < bg.adj_matrix_l.size(); v++) {
        vector<Edge> edges;
        set<time_t> times;

        for (int u = 0; u < bg.adj_matrix_l[v].size(); u++) {
            edges.emplace_back(
                Edge(bg.adj_matrix_l[v][u], bg.timeSection_l[v][u].first, bg.timeSection_l[v][u].second));
            times.insert(bg.timeSection_l[v][u].first);
            times.insert(bg.timeSection_l[v][u].second);
        }

        sort(edges.begin(), edges.end(), [](Edge &a, Edge &b) { return a.start < b.start; });

        vector<time_t> Tlist;             // 增序排列的时间节点
        unordered_map<time_t, int> Tmap;  // 时间节点映射到时间节点序号的哈希表

        for (auto t : times) {
            Tlist.push_back(t);
            Tmap[t] = Tlist.size() - 1;
        }

        vector<pair<int, int>> partition;  // 连接区间
        for (int i = 0; i < edges.size(); i++) {
            auto edge = edges[i];

            pair<int, int> tmp;
            tmp.first = Tmap[edge.start];
            tmp.second = Tmap[edge.end];
            // 当前区间不与连接区间重合或者当前还没划分区间
            if (partition.size() == 0 || partition.back().second < tmp.first) {
                // 不需要考虑最后一条边：如果最后一条边不与之前的partition相交，那么最后一条边就不会包含可达信息
                if (i < edges.size() - 1) {
                    auto edge_ = edges[i + 1];

                    // 这里只验证后一条边是否与当前边连通，如果不连通，则当前边不包含可达信息
                    if (edge.end > edge_.start) {
                        partition.emplace_back(tmp);
                        flow[edge.node].emplace_back(partition.size(),
                                                     make_tuple(v, Tlist[tmp.first], Tlist[tmp.second]));
                        flow[edge.node].emplace_back(partition.size(),
                                                     make_tuple(v, Tlist[tmp.second - 1], Tlist[tmp.second]));
                    }
                }
            } else if (partition.back().second > tmp.first) {
                flow[edge.node].emplace_back(partition.size(), make_tuple(v, Tlist[tmp.first], Tlist[tmp.second]));
                flow[edge.node].emplace_back(partition.size(), make_tuple(v, Tlist[tmp.second - 1], Tlist[tmp.second]));
                partition.back().second = max(partition.back().second, tmp.second);
            }
        }

        for (auto part : partition) {
            g.addVertex(v, Tlist[part.first], Tlist[part.first + 1]);
            bg.timeTable_u[v].push_back(make_pair(Tlist[part.first], Tlist[part.first + 1]));
            for (int k = part.first + 1; k < part.second; k++) {
                g.addVertex(v, Tlist[k], Tlist[k + 1]);
                bg.timeTable_u[v].push_back(make_pair(Tlist[k], Tlist[k + 1]));
                g.addEdge(make_tuple(v, Tlist[k - 1], Tlist[k]), make_tuple(v, Tlist[k], Tlist[k + 1]));
            }
        }
    }

    // 降低查询转换过程的平均时间复杂度
    for (int i = 0; i < bg.timeTable_u.size(); i++)
        sort(bg.timeTable_u[i].begin(), bg.timeTable_u[i].end(),
             [](pair<time_t, time_t> &a, pair<time_t, time_t> &b) { return a.first < b.first; });

    // 构建不同下部点之间的连接边
    for (int u = 0; u < bg.adj_matrix_u.size(); u++) {
        sort(flow[u].begin(), flow[u].end(),
             [](pair<int, tuple<int, time_t, time_t>> &a, pair<int, tuple<int, time_t, time_t>> &b) {
                 return get<1>(a.second) < get<1>(b.second);
             });
        if (flow[u].size() < 3) continue;
        flow[u].erase(flow[u].begin());
        for (int i = 0; i < flow[u].size() / 2; i++) {
            // 单条边时间跨度约束
            if (
                // flow[u][2 * i + 1].second.second - flow[u][2 * i].second.second <= delta &&
                !(get<0>(flow[u][2 * i + 1].second) == get<0>(flow[u][2 * i].second) &&
                  flow[u][2 * i + 1].first == flow[u][2 * i].first))
                g.addEdge(flow[u][2 * i].second, flow[u][2 * i + 1].second);
        }
    }

    return g;
}

Graph transformation(BiGraph &bg) {
    Graph g;
    // flow[i][j] = <partition.id, <v,startTime>>
    vector<vector<pair<int, tuple<int, time_t, time_t>>>> flow(bg.adj_matrix_u.size(),
                                                               vector<pair<int, tuple<int, time_t, time_t>>>());
    for (int v = 0; v < bg.adj_matrix_l.size(); v++) {
        vector<Edge> edges;
        set<time_t> times;

        for (int u = 0; u < bg.adj_matrix_l[v].size(); u++) {
            edges.emplace_back(
                Edge(bg.adj_matrix_l[v][u], bg.timeSection_l[v][u].first, bg.timeSection_l[v][u].second));
            times.insert(bg.timeSection_l[v][u].first);
            times.insert(bg.timeSection_l[v][u].second);
        }

        sort(edges.begin(), edges.end(), [](Edge &a, Edge &b) { return a.start < b.start; });

        vector<time_t> Tlist;             // 增序排列的时间节点
        unordered_map<time_t, int> Tmap;  // 时间节点映射到时间节点序号的哈希表

        for (auto t : times) {
            Tlist.push_back(t);
            Tmap[t] = Tlist.size() - 1;

        }

        vector<pair<int, int>> partition;  // 连接区间
        for (int i = 0; i < edges.size(); i++) {
            auto edge = edges[i];

            pair<int, int> tmp;
            tmp.first = Tmap[edge.start];
            tmp.second = Tmap[edge.end];
            if (partition.size() == 0 || partition.back().second < tmp.first) {
                if (i < edges.size() - 1) {
                    auto edge_ = edges[i + 1];
                    if (edge.end > edge_.start) {
                        partition.emplace_back(tmp);
                        flow[edge.node].emplace_back(partition.size(),
                                                     make_tuple(v, Tlist[tmp.first], Tlist[tmp.second]));
                        flow[edge.node].emplace_back(partition.size(),
                                                     make_tuple(v, Tlist[tmp.second - 1], Tlist[tmp.second]));
                    }
                }
            } else if (partition.back().second > tmp.first) {
                flow[edge.node].emplace_back(partition.size(), make_tuple(v, Tlist[tmp.first], Tlist[tmp.second]));
                flow[edge.node].emplace_back(partition.size(), make_tuple(v, Tlist[tmp.second - 1], Tlist[tmp.second]));
                partition.back().second = max(partition.back().second, tmp.second);
            }
        }

        for (auto part : partition) {
            vector<int> path;
            g.addVertex(v, Tlist[part.first], Tlist[part.first + 1]);
            bg.timeTable_u[v].push_back(make_pair(Tlist[part.first], Tlist[part.first + 1]));
            for (int k = part.first + 1; k < part.second; k++) {
                g.addVertex(v, Tlist[k], Tlist[k + 1]);
                bg.timeTable_u[v].push_back(make_pair(Tlist[k], Tlist[k + 1]));
                g.addEdge(make_tuple(v, Tlist[k - 1], Tlist[k]), make_tuple(v, Tlist[k], Tlist[k + 1]));
            }
        }
    }

    // 降低查询转换过程的平均时间复杂度
    for (int i = 0; i < bg.timeTable_u.size(); i++)
        sort(bg.timeTable_u[i].begin(), bg.timeTable_u[i].end(),
             [](pair<time_t, time_t> &a, pair<time_t, time_t> &b) { return a.first < b.first; });

    // 构建不同下部点之间的连接边
    for (int u = 0; u < bg.adj_matrix_u.size(); u++) {
        sort(flow[u].begin(), flow[u].end(),
             [](pair<int, tuple<int, time_t, time_t>> &a, pair<int, tuple<int, time_t, time_t>> &b) {
                 return get<1>(a.second) < get<1>(b.second);
             });
        if (flow[u].size() < 3) continue;
        flow[u].erase(flow[u].begin());
        for (int i = 0; i < flow[u].size() / 2; i++) {
            // 单条边时间跨度约束
            if (get<1>(flow[u][2 * i + 1].second) - get<1>(flow[u][2 * i].second) <= 4 * 24 * 3600 &&
                flow[u][2 * i + 1].first != flow[u][2 * i].first)
                g.addEdge(flow[u][2 * i].second, flow[u][2 * i + 1].second);
        }
    }

    return g;
}

using Quaduple = tuple<int, time_t, time_t, int>;
struct QuadupleHasher {
    std::size_t operator()(const Quaduple &q) const {
        std::size_t h1 = (std::get<0>(q));
        std::size_t h2 = (std::get<1>(q));
        std::size_t h3 = (std::get<2>(q));
        std::size_t h4 = (std::get<3>(q));
        return ((h1 << 32) | h2) + ((h3 << 32) | h4);
    }
};


struct QuadupleEqual {
    bool operator()(const Quaduple &q1, const Quaduple &q2) const {
        return std::get<0>(q1) == std::get<0>(q2) && std::get<1>(q1) == std::get<1>(q2) &&
               std::get<2>(q1) == std::get<2>(q2) && std::get<3>(q1) == std::get<3>(q2);
    }
};

Graph newTransformation(BiGraph &bg, int delta) {
    Graph g;

    unordered_map<Quaduple, pair<time_t, time_t>, QuadupleHasher, QuadupleEqual> flow;
    for (int v = 0; v < bg.adj_matrix_l.size(); v++) {
        set<time_t> times;
        unordered_map<int, set<int>> events;
        for (int u = 0; u < bg.adj_matrix_l[v].size(); u++) {
            times.insert(bg.timeSection_l[v][u].first);
            times.insert(bg.timeSection_l[v][u].second);
            int k = u + 1;
            events[bg.timeSection_l[v][u].first].insert(k);
            events[bg.timeSection_l[v][u].second].insert(-k);
        }

        vector<time_t> Tlist;
        unordered_map<time_t, int> Tmap;

        for (auto t : times) {
            Tlist.push_back(t);
            Tmap[t] = Tlist.size() - 1;
        }

        for (int i = 0; i < bg.adj_matrix_l[v].size(); i++) {
            int u = bg.adj_matrix_l[v][i];
            time_t start = bg.timeSection_l[v][i].first;
            time_t end = bg.timeSection_l[v][i].second;
            time_t start_ = Tlist[Tmap[start] + 1];
            time_t end_ = Tlist[Tmap[end] - 1];
            flow[make_tuple(u, start, end, v)] = make_pair(start_, end_);
        }

        set<int> current_edges;
        for (int i = 0; i < Tlist.size(); i++) {
            auto set_size = current_edges.size();

            time_t t = Tlist[i];
            for (auto e : events[t]) {
                if (e >= 0) {
                    current_edges.insert(e);
                } else {
                    current_edges.erase(-e);
                }
            }

            if (current_edges.empty())
                continue;
            else {
                bg.timeTable_u[v].push_back(make_pair(Tlist[i], Tlist[i + 1]));
                if (i == 0) {
                    g.addVertex(v, Tlist[i], Tlist[i + 1]);
                } else {
                    g.addVertex(v, Tlist[i], Tlist[i + 1]);
                    if (set_size > 0) {
                        g.addEdge(make_tuple(v, Tlist[i - 1], Tlist[i]), make_tuple(v, Tlist[i], Tlist[i + 1]));
                    }
                }
            }
        }
    }

    for (int u = 0; u < bg.adj_matrix_u.size(); u++) {
        map<int, set<int>> events;
        for (int v = 0; v < bg.adj_matrix_u[u].size(); v++) {
            int k = v + 1;
            events[bg.timeSection_u[u][v].first].insert(k);
            events[bg.timeSection_u[u][v].second].insert(-k);
        }

        vector<int> ends;
        vector<int> starts;
        ends.reserve(50);
        starts.reserve(50);

        int FIND_END = 1;
        int FIND_START = 0;
        time_t max_end = 0;
        time_t min_start = INT64_MAX;

        auto iter = events.begin();
        while (iter != events.end()) {
            auto kv = *iter;
            time_t t = kv.first;
            set<int> &es = kv.second;
            if (FIND_END) {
                for (auto e : es) {
                    if (e < 0) {
                        int k = (-e) - 1;
                        ends.push_back(k);
                    }
                }
                if (!ends.empty()) {
                    max_end = max(max_end, t);
                }
                for (auto e : es) {
                    if (e > 0 && !ends.empty()) {
                        int k = e - 1;
                        starts.push_back(k);
                    }
                }
                if (!starts.empty()) {
                    FIND_START = 1;
                    FIND_END = 0;
                    min_start = min(min_start, t);
                }
                iter++;

            } else if (FIND_START) {
                for (auto e : es) {
                    if (e < 0) {
                        FIND_END = 1;
                        FIND_START = 0;
                        break;
                    }
                }

                // FIND_END=1时加边，否则继续找开始点
                if (FIND_END) {
                    if (starts.size() == 1 && ends.size() == 1) {
                        int v1 = bg.adj_matrix_u[u][starts[0]];
                        int v2 = bg.adj_matrix_u[u][ends[0]];
                        time_t t1 = bg.timeSection_u[u][starts[0]].first;
                        time_t t2 = bg.timeSection_u[u][starts[0]].second;

                        time_t t3 = bg.timeSection_u[u][ends[0]].first;
                        time_t t4 = bg.timeSection_u[u][ends[0]].second;

                        time_t t5 = flow[make_tuple(u, t1, t2, v1)].first;
                        time_t t6 = flow[make_tuple(u, t3, t4, v2)].second;

                        g.addEdge(make_tuple(v2, t6, t4), make_tuple(v1, t1, t5));

                    } else if (starts.size() == 1 && ends.size() > 1) {
                        int v1 = bg.adj_matrix_u[u][starts[0]];
                        time_t t1 = bg.timeSection_u[u][starts[0]].first;
                        time_t t2 = bg.timeSection_u[u][starts[0]].second;
                        time_t t3 = flow[make_tuple(u, t1, t2, v1)].first;

                        g.addVertex(v1, min_start, min_start);

                        g.addEdge(make_tuple(v1, min_start, min_start), make_tuple(v1, t1, t3));

                        for (auto e : ends) {
                            int v2 = bg.adj_matrix_u[u][e];
                            time_t t4 = bg.timeSection_u[u][e].first;
                            time_t t5 = bg.timeSection_u[u][e].second;
                            time_t t6 = flow[make_tuple(u, t4, t5, v2)].second;
                            g.addEdge(make_tuple(v2, t6, t5), make_tuple(v1, min_start, min_start));
                        }

                    } else if (starts.size() > 1 && ends.size() == 1) {
                        int v1 = bg.adj_matrix_u[u][ends[0]];
                        time_t t1 = bg.timeSection_u[u][ends[0]].first;
                        time_t t2 = bg.timeSection_u[u][ends[0]].second;
                        time_t t3 = flow[make_tuple(u, t1, t2, v1)].second;
                        int v_min = INT32_MAX;
                        for (auto s : starts) v_min = min(v_min, bg.adj_matrix_u[u][s]);
                        g.addVertex(v_min, min_start, min_start);

                        g.addEdge(make_tuple(v1, t3, t2), make_tuple(v_min, max_end, max_end));

                        for (auto s : starts) {
                            int v2 = bg.adj_matrix_u[u][s];
                            time_t t4 = bg.timeSection_u[u][s].first;
                            time_t t5 = bg.timeSection_u[u][s].second;
                            time_t t6 = flow[make_tuple(u, t4, t5, v2)].first;

                            g.addEdge(make_tuple(v_min, max_end, max_end), make_tuple(v2, t4, t6));
                        }
                    } else if (starts.size() > 1 && ends.size() > 1) {
                        int v_min = INT32_MAX;
                        for (auto s : starts) v_min = min(v_min, bg.adj_matrix_u[u][s]);

                        g.addVertex(v_min, min_start, min_start);
                        for (auto s : starts) {
                            int v1 = bg.adj_matrix_u[u][s];
                            time_t t1 = bg.timeSection_u[u][s].first;
                            time_t t2 = bg.timeSection_u[u][s].second;
                            time_t t3 = flow[make_tuple(u, t1, t2, v1)].first;

                            g.addEdge(make_tuple(v_min, min_start, min_start), make_tuple(v1, t1, t3));
                        }

                        for (auto e : ends) {
                            int v2 = bg.adj_matrix_u[u][e];
                            time_t t4 = bg.timeSection_u[u][e].first;
                            time_t t5 = bg.timeSection_u[u][e].second;
                            time_t t6 = flow[make_tuple(u, t4, t5, v2)].second;

                            g.addEdge(make_tuple(v2, t6, t5), make_tuple(v_min, min_start, min_start));
                        }
                    }
                    max_end = 0;
                    min_start = INT64_MAX;
                    starts.clear();
                    ends.clear();
                    starts.reserve(50);
                    ends.reserve(50);
                } else {
                    for (auto e : es) {
                        if (e > 0) {
                            int k = e - 1;
                            starts.push_back(k);
                        }
                    }
                    min_start = min(min_start, t);
                    iter++;
                }
            }
        }
    }

    for (int i = 0; i < bg.timeTable_u.size(); i++)
        sort(bg.timeTable_u[i].begin(), bg.timeTable_u[i].end(),
             [](pair<time_t, time_t> &a, pair<time_t, time_t> &b) { return a.first < b.first; });

    return g;
}

vector<queryInfo> readQuery(string filepath) {
    vector<queryInfo> querys;

    fstream fin(filepath, ios::in);
    int u, w;
    time_t start, end;
    while (fin >> u >> w >> start >> end) {
        queryInfo query;
        query.u = u;
        query.w = w;
        query.start = start;
        query.end = end;
        querys.push_back(query);
    }

    fin.close();
    return querys;
}

vector<queryInfo> readQueryForTest(string filepath) {
    vector<queryInfo> querys;

    fstream fin(filepath, ios::in);
    int u, w;
    time_t start, end;
    while (fin >> u >> w) {
        queryInfo query;
        query.u = u;
        query.w = w;
        querys.push_back(query);
    }

    fin.close();
    return querys;
}

/*
 * 通过将有向图投影为无向图寻找独立子图
 */
vector<subGraph> allSubGraph(Graph &g) {
    vector<vector<int>> connection(g.num_vertices() - 1, vector<int>());
    vector<unordered_set<int>> subDagNodes;
    // 避开最后添加的虚拟点
    for (int u = 0; u < g.num_vertices() - 1; u++)
        for (int v : g.out_edges(u)) {
            if (v == g.num_vertices() - 1) continue;
            connection[u].push_back(v);
            connection[v].push_back(u);
        }

    vector<int> nodeToSub(g.num_vertices() - 1, -1);
    unordered_set<int> cliqueID;
    for (int i = 0; i < g.num_vertices() - 1; i++) cliqueID.insert(i);
    while (cliqueID.size() > 0) {
        subDagNodes.push_back(unordered_set<int>());
        int num = subDagNodes.size() - 1;
        queue<int> q;
        q.push(*cliqueID.begin());
        cliqueID.erase(q.front());
        while (q.size() > 0) {
            int u = q.front();
            q.pop();
            subDagNodes.back().insert(u);
            nodeToSub[u] = num;
            for (int v : connection[u])
                if (cliqueID.count(v)) {
                    q.push(v);
                    cliqueID.erase(v);
                }
        }
    }

    vector<subGraph> subGraphs(subDagNodes.size(), subGraph());
    for (int i = 0; i < subDagNodes.size(); i++) {
        subGraphs[i].set_gid(i);
        for (int u : subDagNodes[i]) {
            subGraphs[i].addVertex(g[u]);
            for (int v : g.out_edges(u)) {
                subGraphs[i].addVertex(g[v]);
                subGraphs[i].addEdge(make_tuple(g[u].node, g[u].start, g[u].end),
                                     make_tuple(g[v].node, g[v].start, g[v].end));
            }
        }
    }

    return subGraphs;
}

void printGraph(Graph &graph, string prefix, int id) {
    fstream ftmp("result/" + prefix + "_" + to_string(id) + ".dot", ios::out);
    stringstream stmp;
    stmp << "digraph graphe {\ngraph[rankdir=\"LR\"];\nnodesep = 3;\nnode "
            "[fontname=\"Arial\", shape = \"record\",style=filled,size = "
            "5];\nedge [color=black];"
         << endl;

    for (int i = 0; i < graph.num_vertices(); i++) {
        stmp << "{rank=\"same\";";

        // for (int j = 0; j < depthMap[i].size(); j++) {
        stmp << graph[i].id << "[label = \"id=" << graph[i].id << "|start=" << graph[i].start
             << "|end=" << graph[i].end;
        stmp << "\"]};" << endl;
        // }
        // for (int node : depthMap[i])
        //   for (int u : graph.in_edges(node))
        //     stmp << u << "->" << node << ";" << endl;
        for (auto edge : graph.out_edges(i)) {
            stmp << i << "->" << edge << ";" << endl;
        }
    }

    stmp << "}" << endl;
    ftmp << stmp.str();
    ftmp.close();
}

void printGraph(Tree &tree, string prefix, int depth, int id) {
    vector<int> nodeDepth(tree.num_vertices(), -1);
    vector<vector<int>> depthMap(depth, vector<int>());
    int max_dep = 0;

    // for (int i = 0; i < tree.num_vertices(); i++)
    // {
    // 	int d = dfs(tree, nodeDepth, i);
    // 	max_dep = std::max(max_dep, d);
    // 	if (d < depth)
    // 		depthMap[d].push_back(i);
    // }
    // max_dep = std::min(max_dep, depth);

    fstream ftmp("result/" + prefix + "_" + to_string(id) + ".dot", ios::out);
    stringstream stmp;
    stmp << "digraph graphe {\ngraph[rankdir=\"LR\"];\nnodesep = 3;\nnode "
            "[fontname=\"Arial\", shape = \"record\",style=filled,size = "
            "5];\nedge [color=black];"
         << endl;

    for (int i = 0; i < max_dep; i++) {
        stmp << "{rank=\"same\";";
        for (int j = 0; j < depthMap[i].size() - 1; j++) stmp << to_string(depthMap[i][j]) << ",";
        stmp << depthMap[i].back() << "}";
        for (int j = 0; j < depthMap[i].size(); j++) {
            stmp << depthMap[i][j] << "[label = \"id=" << depthMap[i][j] << "|depth=" << i
                 << "|out=" << tree.out_degree(depthMap[i][j]);
            if (tree.isStart_[depthMap[i][j]] == 1)
                stmp << ";color=red";
            else if (tree.isStart_[depthMap[i][j]] == 2)
                stmp << ";color=green";
            else if (tree.isStart_[depthMap[i][j]] == 3)
                stmp << ";color=blue";
            stmp << "\"];" << endl;
        }
        for (int node : depthMap[i])
            for (int u : tree.in_edges(node)) stmp << u << "->" << node << ";" << endl;
    }

    stmp << "}" << endl;
    ftmp << stmp.str();
    ftmp.close();
}
void printDAG_metis(Graph &graph, string outputPath) {
    fstream fout(outputPath.c_str(), ios::out);
    stringstream ss;
    ss << graph.num_vertices() << " " << graph.num_edges() << endl;
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (int j : graph.out_edges(i)) {
            ss << (j + 1) << " ";
        }
        for (int j : graph.in_edges(i)) {
            ss << (j + 1) << " ";
        }
        ss << endl;
    }
    fout << ss.str();
    fout.close();
}

void printDAG_1(Graph &graph, string outputPath) {
    fstream fout(outputPath.c_str(), ios::out);
    stringstream ss;
    ss << graph.num_vertices() << " " << graph.num_edges() << endl;
    // std::cout << "edge = " << graph.num_edges() << "\n";
    int size = 0;
    for (int i = 0; i < graph.num_vertices(); i++) {
        for (int j : graph.out_edges(i)) {
            ss << (j + 1) << " ";
            size++;
        }
        ss << endl;
    }
    // std::cout << "edge size = " << size << "\n";
    fout << ss.str();
    fout.close();
}

void printDAG_2(Graph &graph, string outputPath) {
    fstream fout(outputPath.c_str(), ios::out);
    stringstream ss;
    ss << "graph_for_greach" << endl;
    ss << graph.num_vertices() << endl;
    for (int i = 0; i < graph.num_vertices(); i++) {
        ss << i << ": ";
        for (int j : graph.out_edges(i)) ss << j << " ";
        ss << "#" << endl;
    }

    fout << ss.str();
    fout.close();
}


std::string extractFileName(const std::string &filePath) {
    // 使用 basename 函数提取文件名
    char *baseName = basename(const_cast<char *>(filePath.c_str()));
    return std::string(baseName);
}
}  // namespace CCR
