#include "InputOutput.h"
#include <fstream>
#include "utility.h"

using namespace xxf_utility;
using namespace std;

namespace min_diff_dp {

UMatrix::UMatrix(const String &path) {
    ifstream ifs;
    ifs.open(path);
    if (!ifs.is_open()) {
        mylog << "无法打开文件" << path <<= logsw_error;
    }
    String line;
    getline(ifs, line);
    istringstream stream(line);
    stream >> nb_set_ele >> nb_subset_ele;
    matrix.resize(nb_set_ele);
    for (int i = 0; i < nb_set_ele; ++i)matrix[i].resize(nb_set_ele, 0.0);
    int node1, node2;
    Distance dis;
    while (getline(ifs, line)) {
        istringstream sst(line);
        sst >> node1 >> node2 >> dis;
        matrix[node1][node2] = dis;
        matrix[node2][node1] = dis;
    }
    mylog << "加载文件" << path <<= logsw_info;
    mylog << "矩阵节点数目为:" << nb_set_ele <<
        "  要选的节点数目为:" << nb_subset_ele <<= logsw_info;
    ifs.close();
}

void Solution::randomInit() {
    int num_select = 0;
    int i = 0;
    while (num_select < nb_subset_ele) {
        while (true) {
            i = myrand.gen(nb_set_ele - 1);
            if (ele_value[i] == 0)break;
        }
        ele_value[i] = 1;
        //select_nodes.insert(i);
        num_select++;
    }
}

void Solution::print() {
    mylog << "目标函数值如下：" << object <<= logsw_info;
    mylog << "选出的节点如下： " <<= logsw_info;
    for (int i = 0; i < nb_set_ele; ++i) {
        if (ele_value[i] == 1)mylog << i << " ";
    }
    mylog <<= logsw_info;
}

}