#ifndef _GRAPH_H
#define _GRAPH_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <list>
#include <map>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#include "bigraph.h"
#define GRAIL

namespace ptree {
using namespace __gnu_cxx;
using namespace std;
struct IdAndSection {
    int id_;
    std::pair<uint32_t, uint32_t> section_;
};

struct Vertex {
    bool fat{};     // fat node
    int path_id{};  // path id
    int dfs_order{};
    int pre_order{};
    int post_order{};
    int first_visit{};  // for test
    bool visited{false};
    int node{};
    time_t start{};
    time_t end{};
    int pnum{};
    int id{};
    int topo_id{};         // topological order
    bool in_tree = false;  // 标记是否出现在树中

    // for matrix
    int partition{-1};
    bool border_in{false};
    bool border_out{false};
    int borderOutIdx{-1};
    int borderInIdx{-1};

#ifdef GRAIL
    // for grail
    int top_level{-1};     // topological level
    int min_parent_level;  // level of the highest parent in top_order
    int min_int;
    long volume;
    double adj_vol;
    double tcs;
    int mingap;
    vector<int> *pre;
    vector<int> *post;
    vector<int> *middle;
#endif
};

// only for grail
struct VertexCompare {
    bool operator()(const Vertex p1, const Vertex p2) const { return p1.id < p2.id; }
};

typedef vector<int> EdgeList;       // edge list represented by vertex id list
typedef vector<Vertex> VertexList;  // vertices list (store real vertex property) indexing by id

struct In_OutList {
    EdgeList inList;
    EdgeList outList;
};

class Hasher {
   public:
    size_t operator()(const tuple<int, time_t, time_t> &p) const {
        // size_t h1 = (get<0>(p));
        size_t h2 = (get<1>(p));
        size_t h3 = (get<2>(p));
        return (h2<<32)|h3;
    }
};

class Equaler {
   public:
    bool operator()(const tuple<int, time_t, time_t> &a, const tuple<int, time_t, time_t> &b) const {
        // 仅比较开始时间和关联的地点，因为同一个地点的划分出的clique的时间区间相互不覆盖
        return (get<0>(a) == get<0>(b)) && (get<1>(a) == get<1>(b)) && (get<2>(a) == get<2>(b));
    }
};

typedef vector<In_OutList> GRA;  // index graph
typedef unordered_map<tuple<int, time_t, time_t>, int, Hasher, Equaler> Cfunc;

class Graph {
   private:
    GRA graph;
    VertexList vl;
    Cfunc f;
    int vsize;

   public:
    Graph();

    Graph(int);

    Graph(const Graph &);

    Graph(string path) {}

    Graph(GRA &, VertexList &);

    Graph(unordered_map<int, vector<int>> &inlist, unordered_map<int, vector<int>> &outlist);

    ~Graph();

    void addVRoot();

    void addVertex(Vertex);

    int getVertexId(Vertex);

    void addVertex(int);

    void addVertex(int node, time_t start, time_t end);

    void addVertex(int node, time_t start, time_t end, int num);

    void setVertexBorderFlag(int node, int type);

    void addEdge(int, int);

    void addEdge(tuple<int, time_t, time_t>, tuple<int, time_t, time_t>);

    vector<int> topological_sort();

    int num_vertices() const;

    int num_edges();

    VertexList &vertices();

    EdgeList &out_edges(int);

    EdgeList &in_edges(int);

    int out_degree(int);

    int in_degree(int);

    int cnum(tuple<int, time_t, time_t>);

    vector<int> getRoots();

    bool hasEdge(int, int);

    Graph &operator=(const Graph &);

    Graph &operator=(Graph &&);

    Vertex &operator[](const int);

    void extract(unordered_map<int, vector<int>> &inlist, unordered_map<int, vector<int>> &outlist);

    void writeGraph(ostream &);

    // for gripp
    vector<int> *GetChild(int node);

    bool CanReach(unsigned nodeA, unsigned nodeB);

    bool search(unsigned nodeA, unsigned nodeB, std::vector<uint8_t> &vis);

    // Graph(unordered_map<int, vector<int>> &inlist, unordered_map<int,
    // vector<int>> &outlist); void extract(unordered_map<int, vector<int>>
    // &inlist, unordered_map<int, vector<int>> &outlist); void
    // printMap(unordered_map<int, vector<int>> &inlist, unordered_map<int,
    // vector<int>> &outlist);

    void reSet() {
        for (int i = 0; i < this->num_vertices(); i++) {
            vl[i].partition = -1;
            vl[i].visited = false;
        }
    }

    void dumpGraph() {
        for (int i = 0; i < num_vertices(); i++) {
            std::cout << "(" << vl[i].node << " " << vl[i].start << "): ";
            for (auto t : out_edges(i)) {
                std::cout << "(" << vl[t].node << " " << vl[t].start << "), ";
            }
            std::cout << "\n";
        }
    }

#ifdef GRAIL
    // void readGraph(istream &);
    void printGraph();
    // void clear();
    void strTrimRight(string &str);
    const double actualgap(const int);
    const double tcs(const int);
    bool incrementalContains(int src, int trg, int cur);
    bool contains(int src, int trg, int dim);
    void printMap(unordered_map<int, vector<int>> &inlist, unordered_map<int, vector<int>> &outlist);
#endif
};

class subGraph : public Graph {
   private:
    int gid;

   public:
    int get_gid();

    void set_gid(int gid);
};
}  // namespace ptree
#endif  //_GRAPH_H
