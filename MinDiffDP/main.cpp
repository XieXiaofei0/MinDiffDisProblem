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
    //TODO：各个参数如何更方便(_param和tabuStep是根据算例修改的)
    const double _param = 0.3;           //表示邻域大小的参数
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
    mylog << "随机种子为：" << t <<= logsw_info;
    LocalSearch local(matrix, _param, sl, tabuStep, sizeTabu, param1, param2, param3);
    //Solution sol = local.solve_tabu();
    Solution sol = local.solve();
    if (!sol.check(matrix))mylog << "目标函数值冲突" <<= logsw_error;
    sol.print();
}

void benchmark(void test(string&)) {
    constexpr int max_num_calculations = 5;         //最大的计算次数
    List<string> instances_id{

        "MDG-a_1_n500_m50",
        //"MDG-a_7_n500_m50",
        //"MDG-a_13_n500_m50",
        //"MDG-a_3_n500_m50",
        //"MDG-a_5_n500_m50",
        //"MDG-a_8_n500_m50",
        //"MDG-a_10_n500_m50",
        //"MDG-a_12_n500_m50",
        //"MDG-a_15_n500_m50",
        //"MDG-a_17_n500_m50",
        //"MDG-a_18_n500_m50",
        //"MDG-a_23_n2000_m200",
   //     "MDG-a_24_n2000_m200",
   // "MDG-a_32_n2000_m200",
   //"MDG-a_35_n2000_m200",
   //     //MDG-b算例
   //     "MDG-b_5_n500_m50",
   //     "MDG-b_11_n500_m50",
   //     "MDG-b_14_n500_m50",
   //     "MDG-b_15_n500_m50",
   // "MDG-b_32_n2000_m200",
   //     //MDG-C算例
   //"MDG-c_3_n3000_m300",
   //      //"MDG-c_7_n3000_m400",
   //      //"MDG-c_11_n3000_m500",
   // "MDG-c_19_n3000_m600",

   //      //DM1A算例
         //"01Type1_52.1_n500m200",
         //"09Type1_52.9_n500m200",
         //"15Type1_52.15_n500m200",
         //SOM-b算例  ---重新设置长度
//"SOM-b_2_n100_m20",
//"SOM-b_4_n100_m40",
//"SOM-b_5_n200_m20",
//"SOM-b_8_n200_m80",
//"SOM-b_9_n300_m30",
//"SOM-b_12_n300_m120",
//"SOM-b_13_n400_m40",
//"SOM-b_16_n400_m160",

         //需要修改参数
         //APOM算例
        // "01a050m10",
        //"02a050m20",
        //"09a250m50",
        //"10a250m100",
        //"13b100m20",
        //"14b100m40",
        //"17b200m40",
        //"18b200m80",
        //"25c150m30",
        //"26c150m60",
        // GKD-b
        // "GKD-b_21_n100_m10",
        // "GKD-b_30_n100_m30",
        // "GKD-b_41_n150_m15",
        // "GKD-b_50_n150_m45",
        // //GKD-c
         //"GKD-c_6_n500_m50",
         //"GKD-c_2_n500_m50",
         //"GKD-c_15_n500_m50",
    };
    for (int i = 0; i < instances_id.size(); ++i) {
        int count = 0;
        for (int c_time = 0; c_time < max_num_calculations; ++c_time) {
            //string file_name = "Instances/MDG-" + instances_id[i] + ".txt";
            //自己电脑-实验室电脑 区别
            //输出到csv文件中
            Date cur;
            string time = cur.shortDateStr();
            ofstream outFile;                                //写文件
            outFile.open("../Deploy/Logs/log.csv", ios::app);
            outFile << time << ',' << instances_id[i] << ',' << count << ',';
            outFile.close();
            //输出结束
            string file_name = "../Deploy/Instances/" + instances_id[i] + ".txt";
            ++count;
            mylog << "第" << count << "次测试算例" << instances_id[i] << "" <<= LogSwitch(1, 1, "BenchMark");
            test(file_name);
        }
    }
}

int main(int argc, char* argv[]) {
    //将详细数据记录到csv文件中
    ifstream infile("../Deploy/Logs/log.csv");                           //判断文件是否存在，写日志文件
    if (!infile.good())                                         //文件不存在
    {
        infile.close();
        ofstream outFile;                                //写文件
        outFile.open("../Deploy/Logs/log.csv", ios::out);
        outFile << "CurrentTime" << ',' << "Instance" << ',' << "Number" << ',' << "Spend_time" << ',' << "nb_nodes" << ',' << "nb_sub_nodes" << ',' << "Seed" << ',' << "Iter" << ',' << "Local_Best_Obj" << endl;
        outFile.close();
    }
    benchmark(test_localSearch);
    system("pause");
    return 0;
}