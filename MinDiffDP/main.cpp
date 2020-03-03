#include<iostream>
#include<iomanip>
#include<time.h>
#include<sstream>
#include<Windows.h>
#include "utility.h"
#include "InputOutput.h"
#include "LocalSearch.h"

using namespace std;
using namespace xxf_utility;
using namespace min_diff_dp;

void test_localSearch(String &filename) {
    //TODO������������θ�����(_param��tabuStep�Ǹ��������޸ĵ�)
    const double _param = 0.3;           //��ʾ�����С�Ĳ���
    const int tabuStep = 35;
    const int sizeTabu = 100000000;
    const double param1 = 1.8;
    const double param2 = 1.9;
    const double param3 = 2.0;
    UMatrix matrix(filename);
    Solution sl(matrix.setele_num(), matrix.subsetele_num());
    time_t t = time(0);
    myrand.setSeed(t);
    //test
    //time_t t = 1581402801;
    //myrand.setSeed(1581402801);
    //test end
    sl.randomInit();
    mylog << "�������Ϊ��" << t <<= logsw_info;
    LocalSearch local(matrix, _param, sl, tabuStep, sizeTabu, param1, param2, param3);
    //Solution sol = local.solve_tabu();
    Solution sol = local.solve();
    if (!sol.check(matrix))mylog << "Ŀ�꺯��ֵ��ͻ" <<= logsw_error;
    sol.print();
}

void benchmark(void test(string&)) {
    constexpr int max_num_calculations = 5;         //���ļ������
    List<string> instances_id{

        //MDG-b-20������
         "MDG-b_10_n500_m50",
         // "MDG-b_11_n500_m50",
         //"MDG-b_12_n500_m50",
         //"MDG-b_15_n500_m50",
         //"MDG-b_16_n500_m50",
         //"MDG-b_17_n500_m50",
         "MDG-b_20_n500_m50",
         //"MDG-b_3_n500_m50",
         ////"MDG-b_5_n500_m50",
         ////"MDG-b_6_n500_m50",
         //"MDG-b_1_n500_m50",
         //"MDG-b_9_n500_m50",
         //////MDG-a-20����
         "MDG-a_1_n500_m50",
         //"MDG-a_5_n500_m50",
         //"MDG-a_8_n500_m50",
         //"MDG-a_9_n500_m50",
         //"MDG-a_14_n500_m50",
         "MDG-a_17_n500_m50",
         //        //MDG-A-40����
                 "MDG-a_22_n2000_m200",
                 "MDG-a_24_n2000_m200",
                 //"MDG-b_21_n2000_m200",
                 //"MDG-b_24_n2000_m200",
                 "10Type1_52.10_n500m200",
                 "13Type1_52.13_n500m200",
                 //"MDG-a_27_n2000_m200",
                 //"MDG-a_32_n2000_m200",
                 "MDG-a_35_n2000_m200",
                 //"MDG-a_38_n2000_m200",
                 //"MDG-a_40_n2000_m200",
         //        //MDG-b-40����
                 "MDG-b_21_n2000_m200",
                 //"MDG-b_24_n2000_m200",
                 "MDG-b_27_n2000_m200",
                 //"MDG-b_33_n2000_m200",
                 "MDG-b_38_n2000_m200",
                 //         //      //DM1A����
         //"10Type1_52.10_n500m200",
         //"13Type1_52.13_n500m200",
         //"15Type1_52.15_n500m200",
         //"16Type1_52.16_n500m200",
         //"19Type1_52.19_n500m200",
         //"04Type1_52.4_n500m200",
         //"07Type1_52.7_n500m200",
         //"01Type1_52.1_n500m200",
         //        ////MDG-C����
                 //"MDG-c_1_n3000_m300",
                 //"MDG-c_5_n3000_m300",
                 // "MDG-c_14_n3000_m500",
                  "MDG-c_17_n3000_m600",
                  //"MDG-c_20_n3000_m600",
                  "MDG-c_11_n3000_m500",
    };
    for (int i = 0; i < instances_id.size(); ++i) {
        int count = 0;
        for (int c_time = 0; c_time < max_num_calculations; ++c_time) {
            //string file_name = "Instances/MDG-" + instances_id[i] + ".txt";
            //�Լ�����-ʵ���ҵ��� ����
            //�����csv�ļ���
            Date cur;
            string time = cur.shortDateStr();
            ofstream outFile;                                //д�ļ�
            outFile.open("../Deploy/Logs/log.csv", ios::app);
            outFile << time << ',' << instances_id[i] << ',' << count << ',';
            outFile.close();
            //�������
            string file_name = "../Deploy/Instances/" + instances_id[i] + ".txt";
            ++count;
            mylog << "��" << count << "�β�������" << instances_id[i] << "" <<= LogSwitch(1, 1, "BenchMark");
            test(file_name);
        }
    }
}

int main(int argc, char* argv[]) {
    //����ϸ���ݼ�¼��csv�ļ���
    ifstream infile("../Deploy/Logs/log.csv");                           //�ж��ļ��Ƿ���ڣ�д��־�ļ�
    if (!infile.good())                                         //�ļ�������
    {
        infile.close();
        ofstream outFile;                                //д�ļ�
        outFile.open("../Deploy/Logs/log.csv", ios::out);
        outFile << "CurrentTime" << ',' << "Instance" << ',' << "Number" << ',' << "Spend_time" << ',' << "nb_nodes" << ',' << "nb_sub_nodes" << ',' << "Seed" << ',' << "Iter" << ',' << "Local_Best_Obj" << endl;
        outFile.close();
    }
    benchmark(test_localSearch);
    system("pause");
    return 0;
}