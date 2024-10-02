# Clique-Center Recusive Label Indexing
---
## 1. Build and Test
use `test.sh` to build the code automatically and test the efficiency of our algorithme：
``` bash
# build with flag O3
./test.sh build
# test all datasets
./test.sh test
# test given dataset: filepath->path of dataset
# data/gfile/*->mode g
# data/rfile/*->mode r
./test.sh test filepath mode 
``` 

## 2. Some details of implementation
### 2.1 Transformation from TBG to DAG
As mentioned in our paper, all cliques associated with only one upper layer vertex are redundant. During the clique construction process, we can calculate the exact number of upper layer vertices each clique is associated with based on the movement of the upper vertex, more specifically, startTime and endTime correspond to the activation and deactivation of the link between the upper vertex and the lower vertex, respectively. 
Therefore, when constructing cliques for the lower layer vertex $v$, we can iterate over $E(v)$ and for $k$ timestamps, $startTime(i)$ and $endTime(i)$ record the number of edges activated and deactivated at time $i \in [1,...,k]$.Then we have the following recurrence relation:
$$clique(i+1).pnum = \begin{cases}
 startTime(i+1) & i = 0 \\ 
 clique(i).pnum+startTime(i+1)-endTime(i+1) & i \in [1,...,k-2] \\
 \end{cases}$$

Note that here we are only considering cliques that belong to the same partition. With $k$ timestamps, we have $k-1$ cliques. After that, we can remove the clique with pnum=1 and add the corresponding edges to the transformed graph.

### 2.2 Mark Red and Black nodes
