#include <gflags/gflags.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
#include "basebfs/bfs.h"
#include "compact.h"
#include "graph.h"
#include "nlohmann/json.hpp"
#include "omp.h"
#include "partition.h"
#include "utils.h"
using namespace std;
using json = nlohmann::json;

double getMemoryUsageGB() {
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.substr(0, 6) == "VmSize") {
            size_t pos = line.find_first_of("0123456789");
            if (pos != std::string::npos) {
                // 将读取的内存大小从 kB 转换为 GB
                size_t memSizeKB = std::stoull(line.substr(pos));
                return memSizeKB / 1048576.0;  // 1 GB = 1048576 kB
            }
        }
    }
    return 0;  // 如果没有找到 VmSize，返回0
}

double maxMem = 0;
void monitorMemory() {
    while (true) {
        double memoryUsageGB = getMemoryUsageGB();  // 读取内存使用量并转换为 GB
        // std::cout << "Current memory usage: " << memoryUsageGB << " GB" << std::endl;
        // LOG("Current memory usage = {} GB", memoryUsageGB);
        maxMem = max(maxMem, memoryUsageGB);
        sleep(10);  // 每10秒执行一次
    }
}

DEFINE_string(inputPath, "", "输入文件路径,绝对路径");
DEFINE_string(queryFilePath, "", "查询测试文件路径,绝对路径");
DEFINE_string(outputPath, "", "输出文件路径,绝对路径");
DEFINE_string(algorithm, "no func", "试图使用的算法");

DEFINE_int32(maxBlock, 64, "最大分区数");
DEFINE_int32(minBlock, 4, "最小分区数");
DEFINE_int32(numStep, 100, "循环次数");
DEFINE_int32(delta, std::numeric_limits<int32_t>::max(), "delta 截断");
DEFINE_int32(numThread, 4, "线程数");

string querypairfilename = "../result/query";

int main(int argc, char *argv[]) {
    if (argc == 1) return 1;
    spdlog::cfg::load_env_levels();
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    auto input_path = FLAGS_inputPath;
    auto output_path = FLAGS_outputPath;
    auto algo = FLAGS_algorithm;
    auto name = CCR::split(input_path, "/.", 100);
    auto fn = CCR::extractFileName(input_path);
    std::thread backThread(monitorMemory);
    backThread.detach();
    // cout << "开始读取 " << fn << endl;
    LOG("开始读取 {}", fn);
    LOG("output path = {}", output_path);
    querypairfilename += "_" + fn;
    BiGraph bg(input_path, "r");

    if (bg.adj_matrix_l.empty()) {
        exit(-1);
    }

    auto pairs = CCR::readQuery(FLAGS_queryFilePath);

    CCR::Timer timer;
    vector<int> queryRes(pairs.size(), 0);
    unsigned bgEdgeNum = 0;
    for (const auto &Neighbors : bg.adj_matrix_l) bgEdgeNum += Neighbors.size();
    timer.start();
    ptree::Graph g = CCR::newTransformation(bg, FLAGS_delta);
    timer.ticker();
    LOG("transformation time = {}(ms)", timer.get_last_consuming());
    LOG("g.num_vertices = {}, g.num_edges = {}", g.num_vertices(), g.num_edges());
    auto transTime = timer.get_last_consuming();
    srand48(time(nullptr));

    Partitioner partitioner(FLAGS_minBlock);
    omp_set_num_threads(FLAGS_numThread);
    // 分区
    partitioner.runPartition(g, algo);
    timer.ticker();
    LOG("run partition time =  {} (ms)", timer.get_last_consuming());

    timer.ticker();
    partitioner.ComputeBorder();  //
    timer.ticker();
    LOG("ComputeBorder = {} (ms)", timer.get_last_consuming());

    timer.ticker();
    partitioner.runLocalReachability();  //
    timer.ticker();
    LOG("local algorithm {}, time = {} (ms)", algo, timer.get_last_consuming());

    timer.ticker();
    partitioner.runBorderReachability();  //
    timer.ticker();
    LOG("border Reachability = {} (ms)", timer.get_last_consuming());

    timer.ticker();
    partitioner.runQueryWithBfs(bg, g, pairs, queryRes);
    timer.end();
    int reach = 0;

    json j;
    j["分区数"] = FLAGS_minBlock;
    j["转换时间ms"] = transTime;
    j["方法"] = algo;
    j["index构建时间"] = timer.get_total_consuming() - timer.get_last_consuming();
    j["并行构建时间"] = "Nan";
    j["边界图大小"] = std::to_string(partitioner.borderGraph.num_vertices()) + " , " +
                      std::to_string(partitioner.borderGraph.num_edges());
    j["内存大小G"] = maxMem;
    j["查询时间ms"] = timer.get_last_consuming();
    j["平均查询时间"] = (double)timer.get_last_consuming() / (double)pairs.size();
    j["index size"] = partitioner.getIndexSize();
    j["graph nodes"] = g.num_vertices();
    j["graph edge"] = g.num_edges();

    double sum = 0;
    for (int i = 0; i < g.num_vertices(); i++) {
        for (auto edge : g.out_edges(i)) {
            if (g[i].partition != g[edge].partition) {
                sum++;
            }
        }
    }
    LOG("edge cut = {}", sum / g.num_edges());
    j["edeg_cut"] = sum / g.num_edges();
    std::ofstream ofile(FLAGS_outputPath);
    ofile << j << std::endl;
    LOG("run {} query time = {}(ms), total time = {}, build index time = {}, max memory usage = {}", pairs.size(),
        timer.get_last_consuming(), timer.get_total_consuming(),
        timer.get_total_consuming() - timer.get_last_consuming(), maxMem);
}
