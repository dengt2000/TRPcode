#ifndef PATITION_H
#define PATITION_H

#include <boost/functional/hash.hpp>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "grail/Grail.h"
#include "graph.h"
#include "lib/reachability/oreach.h"
#include "pathtree/PathTree.h"
#include "utils.h"

template <typename T>
std::size_t hash_vector(const std::vector<T> &vec) {
    std::size_t seed = 0;
    for (const auto &i : vec) {
        boost::hash_combine(seed, i);
    }
    return seed;
}

class Block {
   private:
   public:
    // key = 等价类编号, value = 虚拟节点的id
    // 第i个等价类 对应的虚拟节点id
    std::unordered_map<int, int> sout;
    void allocateID(ptree::Graph &borderGraph);
    // 标识等价类是否被遍历过
    std::vector<int> eqClassVisitedO;
    std::vector<int> eqClassVisitedI;
    void setEqVisitedO(int gid, int idx) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderOutIdx;
        this->eqClassVisitedO[eqId] = idx;
    }
    void setEqVisitedI(int gid, int idx) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderInIdx;
        this->eqClassVisitedI[eqId] = idx;
    }

    int getEqVisitedO(int gid) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderOutIdx;
        return this->eqClassVisitedO[eqId];
    }

    int getEqVisitedI(int gid) {
        auto lid = global2local[gid];
        auto eqId = localGraph[lid].borderInIdx;
        return this->eqClassVisitedI[eqId];
    }

    uint64_t index_size{0};
    std::unordered_map<int, int> global2local;
    // std::vector<int> global2local;
    std::vector<int> local2global;
    ptree::Graph localGraph;
    int crossOutEdge{0};
    int crossInEdge{0};
    int blockID{-1};
    // PathTree pt;
    unique_ptr<ptree::PathTree> pt;

    // 存的是边界节点的id, 不是gid也不是local id
    std::vector<std::vector<int>> borderClassOut;
    std::vector<std::vector<int>> borderClassIn;
    void buildBorderLink(ptree::Graph &borderGraph, std::unordered_map<int, int> &borderGlobal2Local);
    void fetchClass(std::unordered_map<int, int> &borderGlobal2Local);
    void fetchClassByPT(std::unordered_map<int, int> &borderGlobal2Local);
    void fetchClassByBorderBFS(std::unordered_map<int, int> &borderGlobal2Local);
    uint64_t getIndexSize() { return pt->getIndexSize(); }
    std::unordered_set<int> borderOut;
    std::unordered_set<int> borderIn;
    int getInBorder(int gid, std::vector<int> &que) {
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderInIdx;
        if (idx < 0) {
            return 0;
        }
        std::copy(borderClassIn[idx].begin(), borderClassIn[idx].end(), que.begin());
        return borderClassIn[idx].size();
    }

    int getOutBorder(int gid, std::vector<int> &que) {
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderOutIdx;
        if (idx < 0) {
            return 0;
        }
        std::copy(borderClassOut[idx].begin(), borderClassOut[idx].end(), que.begin());

        return borderClassOut[idx].size();
    }

    bool getInBorder2(int gid, std::vector<int> &que, std::vector<int> &visited, int IDX, int &q2hi) {
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderInIdx;
        if (idx < 0) {
            return 0;
        }
        if (eqClassVisitedI[idx] == -IDX) {
            return false;
        }
        eqClassVisitedI[idx] = -IDX;
        for (auto &t : borderClassIn[idx]) {
            // 节点已经在队列中了
            if (visited[t] == -IDX) {
                continue;
            }
            if (visited[t] == IDX) {
                return true;
            }
            visited[t] = -IDX;
            que[q2hi++] = t;
        }
        return false;
    }

    bool getOutBorder2(int gid, std::vector<int> &que, std::vector<int> &visited, int IDX, int &q1hi) {
        // 先看等价类有没有被访问过
        // 再看每个点
        auto lid = global2local[gid];
        auto idx = localGraph[lid].borderOutIdx;
        if (idx < 0) {
            return 0;
        }
        if (eqClassVisitedO[idx] == IDX) {
            return false;
        }
        eqClassVisitedO[idx] = IDX;
        for (auto &t : borderClassOut[idx]) {
            // 节点已经在队列中了
            if (visited[t] == IDX) {
                continue;
            }
            if (visited[t] == -IDX) {
                return true;
            }
            visited[t] = IDX;
            que[q1hi++] = t;
        }
        return false;
    }

    void dumpGraph() {
        for (int i = 0; i < global2local.size(); i++) {
            cout << i << " : " << global2local[i] << "\n";
        }
        std::cout << "edge = \n";
        localGraph.dumpGraph();
    }
    void setBlockID(int id) { blockID = id; }
    int getBlockID() const { return blockID; }
    int getSize() const { return localGraph.num_vertices(); }
    int getCrossOutEdge() const { return crossOutEdge; }
    void setBorderType(int gid, int type) {
        auto lid = this->global2local[gid];
        localGraph.setVertexBorderFlag(lid, type);
        if (type == 1) {
            this->borderOut.insert(lid);
        } else if (type == 2) {
            this->borderIn.insert(lid);
        }
    }
    void addNode(ptree::Vertex v) {
        localGraph.addVertex(v);
        int len = localGraph.num_vertices();
        // global2local.insert(make_pair(v.id, len - 1));
        global2local[v.id] = len - 1;
        local2global[len - 1] = v.id;
    }

    void addEdge(int gfirst, int gsecond) {
        if (global2local[gfirst] < 0 or global2local[gsecond] < 0) {
            ERROR("add edge error, not in the same blk");
            exit(-1);
        }
        int localA = global2local[gfirst];
        int localB = global2local[gsecond];
        localGraph.addEdge(localA, localB);
    }

    void increaseCrossInedge() { crossInEdge++; }
    void increaseCrossOutedge() { crossOutEdge++; }
    void decreaseCrossInedge() { crossInEdge--; }
    void decreaseCrossOutedge() { crossOutEdge--; }
    virtual bool query(int src, int dst) {
        auto s = global2local[src];
        auto d = global2local[dst];
        return pt->reach(s, d);
    }
    virtual void runReachability();
};

class PTBlock : public Block {
    // 父类默认使用pathtree，子类不做任何修改
};

class OReachBlock : public Block {
   public:
    oreach *r{nullptr};
    void runReachability() override;
    bool query(int src, int dst) override;
};

class GrailBlock : public Block {
   public:
    grail::Grail *grailinst{nullptr};
    ExceptionList *el{nullptr};
    void runReachability() override;
    bool query(int src, int dst) override;
};

class Partitioner {
   public:
    Partitioner() = default;

    Partitioner(int maxs, int mins) : maxBlkSize(maxs), minBlkSize(mins) {}
    Partitioner(int numBlk) : minBlkSize(numBlk) {}

    ~Partitioner() = default;

    void runPartition(ptree::Graph &graph, const string &algo);
    uint64_t getIndexSize() {
        uint64_t res = this->index_size;
        for (auto t : blocks) {
            res += t->getIndexSize();
        }
        return res;
    }
    // tools
    void PrintPartitionRatio(ptree::Graph &graph) const;
    void ComputeBorder();

    void runLocalReachability();
    void runBorderReachability();
    void runQuery(ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo);
    void runQueryWithBfs(BiGraph &bg, ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo,
                         std::vector<int> &queryRes);
    void runQueryWithOReach(ptree::Graph &graph, std::vector<CCR::queryInfo> &queryInfo);
    std::vector<Block *> blocks;
    std::unordered_map<int, int> borderGlobal2Local;
    // std::vector<int> borderGlobal2Local;
    ptree::Graph borderGraph;
    uint64_t index_size{0};

   private:
    vector<int> que1;
    vector<int> que2;
    vector<int> inner_que1;
    vector<int> inner_que2;
    vector<int> visited_inner;
    vector<int> visited_outter;
    int maxBlkSize{0};
    int minBlkSize{0};
    unique_ptr<ptree::PathTree> pt;
    void getInBorder(int &q2hi, int dst, ptree::Graph &graph, int idx);
    void getOutBorder(int &q1hi, int src, ptree::Graph &graph, int idx);
    bool getInBorder(int &q2hi, vector<int> dsts, ptree::Graph &graph, int idx);
    bool getOutBorder(int &q1hi, vector<int> srcs, ptree::Graph &graph, int idx);
    bool runBiBFS(int q1hi, int q2hi, int &sum, int idx);
};

// for oreach

#endif