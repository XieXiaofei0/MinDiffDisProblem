#include<iostream>
#include<iomanip>
#include<sstream>
#include "utility.h"
#include "InputOutput.h"
#include "LocalSearch.h"

using namespace std;
using namespace xxf_utility;
using namespace min_diff_dp;

void test_localSearch(String &filename) {
    //TODO����������
    double _param = 0.3;           //��ʾ�����С�Ĳ���
    UMatrix matrix(filename);
    Solution sl(matrix.setele_num(), matrix.subsetele_num());
    sl.randomInit();
    mylog << "�ܽڵ�����" << sl.get_nb_nodes()
        <<"��ѡ�Ľڵ�����"<< sl.get_sub_nb_nodes() <<= logsw_info;
    time_t t = time(0);
    myrand.setSeed(t);
    mylog << "������ӣ�" << t <<= logsw_info;
    LocalSearch local(matrix, _param, sl);
    Solution sol = local.solve();
    sol.print();
}

void benchmark(void test(string&)) {
    constexpr int max_num_calculations = 3;          //���ļ������
    List<string> instances_id{
        "a_1_n500_m50","a_2_n500_m50",
        "a_3_n500_m50","a_4_n500_m50"
        "a_21_n2000_m200","a_22_n2000_m200",
        "a_23_n2000_m200","a_24_n2000_m200"
    };
    for (int i = 0; i < instances_id.size(); ++i) {
        int count = 0;
        for (int c_time = 0; c_time <= max_num_calculations; ++c_time) {
            //string file_name = "Instances/MDG-" + instances_id[i] + ".txt";
            //�Լ�����-ʵ���ҵ������� why?
            string file_name = "../Deploy/Instances/MDG-" + instances_id[i] + ".txt";
            ++count;
            mylog << "��" << count << "�β�������" << instances_id[i] << "" <<= LogSwitch(1, 1, "BenchMark");
            test(file_name);
        }
    }
}

int main(int argc, char* argv[]) {
    benchmark(test_localSearch);
    system("pause");
    return 0;
}