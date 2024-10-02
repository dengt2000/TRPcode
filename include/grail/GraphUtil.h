/* Copyright (c) Hilmi Yildirim 2010,2011.

The software is provided on an as is basis for research purposes.
There is no additional support offered, nor are the author(s)
or their institutions liable under any circumstances.
*/
#ifndef Graph_UTIL_H_
#define Graph_UTIL_H_

#include <sys/time.h>

#include <unordered_map>

#include "graph.h"

namespace grail {
class GraphUtil {
   public:
    static void dfs(ptree::Graph &g, int vid, vector<int> &preorder, vector<int> &postorder, vector<bool> &visited);
    static void topo_leveler(ptree::Graph &g);
    static int topo_level(ptree::Graph &g, int vid);
    static void topological_sort(ptree::Graph g, vector<int> &ts);
    static void tarjan(ptree::Graph &g, int vid, int &index, unordered_map<int, pair<int, int> > &order,
                       vector<int> &sn, multimap<int, int> &sccmap, int &scc);
    static void mergeSCC(ptree::Graph &g, int *on, vector<int> &ts);
    static void traverse(ptree::Graph &tree, int vid, int &pre_post, vector<bool> &visited);
    static void pre_post_labeling(ptree::Graph &tree);

    static void genRandomGraph(int n, double c, char *filename);
};
}  // namespace grail
#endif