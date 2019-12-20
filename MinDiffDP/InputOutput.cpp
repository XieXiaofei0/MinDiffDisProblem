#include "InputOutput.h"
#include <fstream>
#include "utility.h"

using namespace xxf_utility;
using namespace std;

namespace min_diff_dp {

//List<List<Distance> > _matrix;         //attention:ȫ�ֱ����Ķ��壺ֻ�ܶ���һ�Σ�

bool judgeDistanceEqual(const Distance &a, const Distance &b) {    //���a��b�󣬷���true
    return fabs(a - b) < EPS;
    //return a > b + EPS;
}

UMatrix::UMatrix(const String &path) {
    ifstream ifs;
    ifs.open(path);
    if (!ifs.is_open()) {
        mylog << "�޷����ļ�" << path <<= logsw_error;
    }
    String line;
    getline(ifs, line);
    istringstream stream(line);
    stream >> nb_set_ele >> nb_subset_ele;
    //(nb_set_ele, List<Distance>(nb_set_ele, 0.0));
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
    mylog << "�����ļ�" << path <<= logsw_info;
    mylog << "����ڵ���ĿΪ:" << nb_set_ele <<
        "  Ҫѡ�Ľڵ���ĿΪ:" << nb_subset_ele <<= logsw_info;
    ifs.close();
}

void Solution::randomInit() {
    int num_select = 0;
    int i = 0;
    srand(myrand.getSeed());
    while (num_select < nb_subset_ele) {
        while (true) {
            i = rand() % nb_set_ele;
            //i = myrand.gen(nb_set_ele - 1);
            if (ele_value[i] == 0)break;
        }
        ele_value[i] = 1;
        num_select++;
    }
}

bool Solution::check(const UMatrix &_matrix) {
    List<int> select_nodes;
    select_nodes.reserve(nb_subset_ele);
    for (int i = 0; i < nb_set_ele; ++i) {
        if (ele_value[i])select_nodes.push_back(i);
    }
    Distance min = DISTANCE_MAX, max = 0;
    for (int i = 0; i < nb_subset_ele; ++i) {
        Distance cur = 0;
        for (int j = 0; j < nb_subset_ele; ++j) {
            cur += _matrix.dis_nodes(select_nodes[i], select_nodes[j]);
        }
        if (cur > max)max = cur;
        if (cur < min)min = cur;
    }
    Distance obj = max - min;
    mylog << "����Ŀ�꺯��ֵΪ��" << obj <<= logsw_info;
    mylog << "��ǰ��Ŀ�꺯��ֵΪ��" << object <<= logsw_info;
    return judgeDistanceEqual(obj, object);
}

void Solution::print() {
    mylog << "Ŀ�꺯��ֵ���£�" << object <<= logsw_info;
    mylog << "ѡ���Ľڵ����£� " <<= logsw_info;
    for (int i = 0; i < nb_set_ele; ++i) {
        if (ele_value[i] == 1)mylog << i << " ";
    }
    mylog <<= logsw_info;
}

}