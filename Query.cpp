#include <string.h>
#include <fstream>
#include <iostream>
#include "utils.h"

using namespace std;

int main(int argv, char *argc[]) {
    if (argv != 10) {
        cout << "请给出一下参数：" << endl;
        cout << "\t--from：出发节点" << endl;
        cout << "\t--to：到达节点" << endl;
        cout << "\t--start：开始时间" << endl;
        cout << "\t--end：开始时间" << endl;
        cout << "\t--graphfile：链分解文件" << endl;
        cout << "\t--chainfile：链分解文件" << endl;
        cout << "\t--timefile：时间戳文件" << endl;
        cout << "\t--linfile：入表文件" << endl;
        cout << "\t--loutfile：出表文件" << endl;
        cout << "示例：./Query.exe --from=1 --to=3 --start=2011/01/01:12:15:14 "
                "--end=2011/02/01:12:15:14 --graphfile= "
                "--chainfile=result/Chains.txt --timefile=result/timeTable.txt "
                "--linfile=result/Lin.txt --loutfile=result/Lout.txt"
             << endl;
        return 0;
    }

    map<string, string> cmd;
    for (int i = 1; i < argv; i++) {
        auto slices = split(argc[i], "-=", 50);
        cout << slices[0] << ":" << slices[1] << endl;
        cmd[slices[0]] = slices[1];
    }

    time_t start, end;
    start = date_to_unixtime(cmd["start"]);
    end = date_to_unixtime(cmd["end"]);
    if (start >= end) {
        cout << "时间区间错误" << endl;
        return 0;
    }

    fstream timeFile(cmd["timefile"], ios::in);
    fstream LinFile(cmd["linfile"], ios::in);
    fstream LoutFile(cmd["loutfile"], ios::in);
    fstream ChainFile(cmd["chainfile"], ios::in);

    BiGraph g(cmd["graphfile"]);

    string buffer;
    vector<time_t> timeTable;
    while (getline(timeFile, buffer)) {
        auto slices = split(buffer, "\t", 30);
        timeTable.push_back(atoi(slices[1].c_str()));
    }

    io_index io;
    while (getline(LoutFile, buffer)) {
        auto slices = split(buffer, "\t", 1000);
        auto slice_0 = split(slices[0], ",", 50);
        int from = atoi(slice_0[0].c_str()) * g.adj_matrix_l.size() + atoi(slice_0[1].c_str());
        for (auto iter = ++slices.begin(); iter != slices.end(); iter++) {
            auto slice = split(*iter, ",", 50);
            int to = atoi(slice[0].c_str()) * g.adj_matrix_l.size() + atoi(slice[1].c_str());
            io.first[from].insert(to);
        }
    }
    while (getline(LinFile, buffer)) {
        auto slices = split(buffer, "\t", 1000);
        auto slice_0 = split(slices[0], ",", 50);
        int to = atoi(slice_0[0].c_str()) * g.adj_matrix_l.size() + atoi(slice_0[1].c_str());
        for (auto iter = ++slices.begin(); iter != slices.end(); iter++) {
            auto slice = split(*iter, ",", 50);
            int from = atoi(slice[0].c_str()) * g.adj_matrix_l.size() + atoi(slice[1].c_str());
            io.second[to].insert(from);
        }
    }

    Chains chains;
    getline(ChainFile, buffer);
    auto slices = split(buffer, " ", 40);
    int timeN = atoi(slices[1].c_str());
    int posN = atoi(slices[2].c_str());
    chains.second = vector<pair<int, int>>(timeN * posN, pair<int, int>(0, 0));
    int i = 0;
    while (getline(ChainFile, buffer)) {
        vector<int> chain;
        auto slices = split(buffer, "\t", 1000);
        int j = 0;
        for (auto iter = slices.begin(); iter != slices.end(); iter++) {
            auto coord = split(*iter, ",", 10);
            int node = atoi((coord[0]).c_str()) * posN + atoi((coord[1]).c_str());
            chains.second[node].first = j;
            chains.second[node].second = i;
            chain.push_back(node);
            j++;
        }
        i++;
        chains.first.push_back(chain);
    }

    cout << "读取完成" << endl;

    int from, to;
    from = atoi(cmd["from"].c_str());
    to = atoi(cmd["to"].c_str());

    cout << "开始查询" << endl;

    if (query(from, to, start, end, timeTable, io, chains, g))
        cout << "可达" << endl;
    else
        cout << "不可达" << endl;
}