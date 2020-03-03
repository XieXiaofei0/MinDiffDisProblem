#include"LocalSearch.h"
#include"utility.h"
#include"InputOutput.h"
#include <queue>

using namespace std;
using namespace xxf_utility;

namespace min_diff_dp {

LogSwitch logsw_local(1, 1, "PINCOL_PT");

bool compareByAscend(const pair<int, Distance> &a, const  pair<int, Distance> &b) {
    return a.second < b.second;
}
bool compareByDescend(const pair<int, Distance> &a, const  pair<int, Distance> &b) {
    return a.second > b.second;
}

LocalSearch::LocalSearch(const UMatrix &_matrix, const double param, const Solution &_init_sol, const int _tabu_step,      //xxf:done right--12.10
    const int _size_of_tabu, const double _param1, const double _param2, const double _param3) :
    ins(_matrix), max_select_node(-1), min_select_node(DISTANCE_MAX), 
    best_max_select_node(-1), best_min_select_node(DISTANCE_MAX),
    best_hashfun_one(-1),best_hashfun_two(-1),best_hashfun_three(0),
    nb_nodes(_matrix.setele_num()),nb_sub_nodes(_matrix.subsetele_num()), 
    iter(0),                   //xxf:修改算例运行最长时间ms
    max_time(_matrix.setele_num() *1.5* 1000),
    rate_of_sele_nodes(param), tabu_step(_tabu_step), size_of_tabu_list(_size_of_tabu),
    hashFun_one_param(_param1), hashFun_two_param(_param2), hashFun_three_param(_param3)
{
    node_value.resize(nb_nodes);
    node_value = _init_sol.node_values();
    local_best.resize(nb_nodes);
    local_best = _init_sol.node_values();
    tabu_list_one.resize(size_of_tabu_list, 0);
    tabu_list_two.resize(size_of_tabu_list, 0);
    tabu_list_three.resize(size_of_tabu_list, 0);
    hash_key_temp_one.resize(nb_nodes, 0);
    hash_key_temp_two.resize(nb_nodes, 0);
    hash_key_temp_three.resize(nb_nodes, 0);
    init();
}

void LocalSearch::init() {              //xxf:done,right-12.10;临时哈希函数中间键值无错-12.18
    size_neighbor_struc = (int)((nb_nodes - nb_sub_nodes)*rate_of_sele_nodes);   //初始化List大小和参数
    node_dis_sum.reserve(nb_nodes);
    best_solu_dis_sum.reserve(nb_nodes);
    select_nodes.reserve(nb_sub_nodes);
    no_select_nodes.reserve(nb_nodes - nb_sub_nodes);
    long long sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
    for (int i = 0; i < nb_nodes; ++i) {                         //初始化node_dis_sum和目标函数
        int one = (int)(floor(pow(i, hashFun_one_param)));
        int two = (int)(floor(pow(i, hashFun_two_param)));
        int three = (int)(floor(pow(i, hashFun_three_param)));
        hash_key_temp_one[i] = one;
        hash_key_temp_two[i] = two;
        hash_key_temp_three[i] = three;
        double dis = 0.0;
        for (int j = 0; j < nb_nodes; ++j) {
            if (node_value[j]) {
                dis += ins.dis_nodes(i, j);
            } 
        }
        node_dis_sum.push_back(dis);
        best_solu_dis_sum.push_back(dis);
        if (node_value[i]) {
            if (max_select_node < dis)max_select_node = dis;
            if (min_select_node > dis)min_select_node = dis;
            sum1 += hash_key_temp_one[i];                      //计算当前历史最优解三个哈希函数的键值
            sum2 += hash_key_temp_two[i];
            sum3 += hash_key_temp_three[i];
        }
    }
    best_hashfun_one = sum1 % size_of_tabu_list;
    best_hashfun_two = sum2 % size_of_tabu_list;
    best_hashfun_three = sum3 % size_of_tabu_list;
    cur_obj = max_select_node - min_select_node;
    local_best_obj = cur_obj;
    Distance temp_sum = (max_select_node + min_select_node) / 2;
    for (int i = 0; i < nb_nodes; ++i) {                                     //初始化两个辅助结构no_select_nodes和select_nodes
        Distance temp = fabs(node_dis_sum[i] - temp_sum);
        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
        else no_select_nodes.push_back(make_pair(i, temp));
    }
    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //未选中的升序排列
    length_disturbance = (int)(nb_nodes / 2 * nb_sub_nodes);
}

Solution LocalSearch::solve() {
    Timer time(max_time);
    int fixed_value = 4000;   //MDG_a.b_1-20;GKD_c;
    //int fixed_value = 2000;    //MDG_a.b_21-40;DM1A
    //int fixed_value = 500;     //MDG_c
    int step_length = 0;           //记录多少步之内不能改进历史最优解
    bool tabu_flag = false;          //判断是否邻域解都在禁忌中
    bool length_flag = false;
    clock_t start_time = clock();
    while (!time.isTimeOut()) {           //当时间未超时时，进行局部搜索 
        int _hashfun_one = best_hashfun_one;
        int _hashfun_two = best_hashfun_two;
        int _hashfun_three = best_hashfun_three;
        //随机扰动
        if (length_flag) {
            for (int i = 0; i < length_disturbance; ++i)
                    stochastic_perturbation(_hashfun_one, _hashfun_two, _hashfun_three, max_select_node, min_select_node);
            length_flag = false;
            step_length = 0;
        }
        else {
            node_value = local_best;               //强化搜索策略：用历史最优解更新当前解
            cur_obj = local_best_obj;
            //若历史最优解!=当前解,则要更新当前解的部分数据结构node_dis_sum，max_select_node， min_select_node，select_nodes和no_select_nodes
            max_select_node = best_max_select_node;                   //xxf:解决大bug1：之前替换当前解时，未更新相关数据结构node_dis_sum，select_nodes和no_select_nodes
            min_select_node = best_min_select_node;
            no_select_nodes.clear();            //更新辅助结构select_nodes和no_select_nodes
            select_nodes.clear();
            Distance temp_sum = (best_max_select_node + best_min_select_node) / 2;
            for (int i = 0; i < nb_nodes; ++i) {
                node_dis_sum[i] = best_solu_dis_sum[i];
                Distance temp = fabs(best_solu_dis_sum[i] - temp_sum);
                if (node_value[i])select_nodes.push_back(make_pair(i, temp));
                else no_select_nodes.push_back(make_pair(i, temp));
            }
        }
        sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //未选中的升序排列

        int count = 0;

        while (count <= tabu_step) {
            pair<int, int> swap_pair(-1, -1);              //邻域结构中找最好的非禁忌解:保存非禁忌的最好的交换对;第一个I1-->I0，第二个I0-->I1
            pair<Distance, Distance> new_obj(DISTANCE_MAX, 0);         //保存对应的目标函数的最大距离和最小距离
            find_best_move(swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three);   //_hashfun_one已经是更新解的哈希函数值
            //TODO:判断是否所有邻域解都在禁忌中
            if (swap_pair.first == -1) {       //判断是否邻域解都在禁忌中
                mylog << "邻域解都在禁忌中；当前iter：" << iter <<= logsw_local;
                tabu_flag = true;
                break;
            } 
            if (update_solu(swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three, step_length))count = 0;   //更新当前解、历史最优解、count
            else count++; 
            iter++;
            if (step_length == fixed_value) {
                length_flag = true;
                step_length = 0;
                break;
            }
        }
        if (tabu_flag)break;
    }
    mylog << "\n总迭代步数：" << iter <<= logsw_local;
    clock_t end_time = clock();
    //关键信息输出到文件中
    ofstream outFile;                                
    outFile.open("../Deploy/Logs/log.csv", ios::app);
    outFile << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << ',' << nb_nodes << ',' << nb_sub_nodes << ',' << myrand.getSeed() << ',' << iter << ',' << local_best_obj << endl;
    outFile.close();
    //文件输出结束
    return Solution(nb_nodes, nb_sub_nodes, local_best, local_best_obj);
}

void LocalSearch::find_best_move(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three){   //xxf:done,right--12.10
    int new_hashone = _hash_one;     //保存当前解的哈希函数值
    int new_hashtwo = _hash_two;
    int new_hashthree = _hash_three;
    int num = 1;                             //保存邻域解相等的个数
    for (int i = 0; i < nb_sub_nodes; ++i)
    {
        int one_toZero_node = select_nodes[i].first;
        for (int j = 0; j < size_neighbor_struc; ++j) {
            int zero_toOne_node = no_select_nodes[j].first; 
            int _new_hash_one = new_hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //计算交换后新解的哈希函数值
            int _new_hash_two = new_hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
            int _new_hash_three = new_hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
            _new_hash_one = (_new_hash_one+ size_of_tabu_list) % size_of_tabu_list;                //xxf：解决bug2：防止出现负数和>size_of_tabu_list的数
            _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
            _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
            if (tabu_list_three[_new_hash_three]) {       //xxf：两个或三个哈希函数为了减小冲突碰撞；一个哈希函数容易不同解映射到同一个value值，且很容易出现邻域解都被禁忌情况，因为是从强化搜索策略开始
                if (tabu_list_two[_new_hash_two])
                    if (tabu_list_one[_new_hash_one]) {
                        continue;
                    }
            }
            Distance temp_min = node_dis_sum[zero_toOne_node] - ins.dis_nodes(zero_toOne_node, one_toZero_node);      //记录当前一次动作的最大值和最小值;初始化为交换之后新的选中节点的距离之和
            Distance temp_max = temp_min;
            for (int k = 0; k < nb_sub_nodes; ++k) {
                if (k == i)continue;
                int node = select_nodes[k].first;
                Distance update_dis = node_dis_sum[node] - ins.dis_nodes(node, one_toZero_node) + ins.dis_nodes(node, zero_toOne_node);
                if (update_dis > temp_max)temp_max = update_dis;
                else if (temp_min > update_dis)temp_min = update_dis;
                else;
            }
            if ((_new_obj.first - _new_obj.second) > (temp_max - temp_min)) {   //更新邻域动作
                _pair.first = one_toZero_node;
                _pair.second = zero_toOne_node;
                _new_obj.first = temp_max;
                _new_obj.second = temp_min;
                _hash_one = _new_hash_one;
                _hash_two = _new_hash_two;
                _hash_three = _new_hash_three;
            }
        }
    }
}

bool LocalSearch::update_solu(const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, int &step) {  //xxf:done,right-12.10
    bool flag = false;         //表示是否更新历史最优解
    node_value[_pair.first] = 0;                       //更新当前解
    node_value[_pair.second] = 1;
    cur_obj = _new_obj.first - _new_obj.second;
    max_select_node = _new_obj.first;
    min_select_node = _new_obj.second;
    if (local_best_obj > cur_obj)       //更新历史最优解，更新历史最优解的相关结构
    {
        local_best = node_value;  //如果能改进历史最优解，则更新历史最优解,返回true
        local_best_obj = cur_obj;
        best_hashfun_one = _hash_one;
        best_hashfun_two = _hash_two;
        best_hashfun_three = _hash_three;
        flag = true;
        step = 0;
    }
    step++;
    tabu_list_one[_hash_one] = 1;         //更新三个禁忌列表
    tabu_list_two[_hash_two] = 1;
    tabu_list_three[_hash_three] = 1;
    no_select_nodes.clear();           //更新辅助结构select_nodes和no_select_nodes
    select_nodes.clear();
    Distance temp_sum = (max_select_node + min_select_node) / 2.0;
    if (flag) {           //若更新了历史最优解，则保存历史最优解的相关数据结构best_solu_dis_sum
        best_max_select_node = max_select_node;
        best_min_select_node = min_select_node;
        for (int i = 0; i < nb_nodes; ++i) {
            node_dis_sum[i] = node_dis_sum[i] - ins.dis_nodes(i, _pair.first) + ins.dis_nodes(i, _pair.second);
            best_solu_dis_sum[i] = node_dis_sum[i];
            Distance temp = fabs(node_dis_sum[i] - temp_sum);
            if (node_value[i])select_nodes.push_back(make_pair(i, temp));
            else no_select_nodes.push_back(make_pair(i, temp));
        }
    }
    else {
        for (int i = 0; i < nb_nodes; ++i) {
            node_dis_sum[i] = node_dis_sum[i] - ins.dis_nodes(i, _pair.first) + ins.dis_nodes(i, _pair.second);
            Distance temp = fabs(node_dis_sum[i] - temp_sum);
            if (node_value[i])select_nodes.push_back(make_pair(i, temp));
            else no_select_nodes.push_back(make_pair(i, temp));
        }
    }
    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //未选中的升序排列
    return flag;
}

//另一种随机扰动方法:从当前历史最优解的候选解里随机挑
bool LocalSearch::stochastic_perturbation(int &hashone, int &hashtwo, int &hashthree, Distance &max, Distance &min) {
    int one_toZero_node = 0;
    int zero_toOne_node = 0;
    int i = 0;
    int _new_hash_one = 0;
    int _new_hash_two = 0;
    int _new_hash_three = 0;
    while (true) {
        srand((unsigned)time(NULL));
        i = rand() % nb_sub_nodes;
        int j = rand() % size_neighbor_struc;
        one_toZero_node = select_nodes[i].first;
        zero_toOne_node = no_select_nodes[j].first;
        int _new_hash_one = hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //计算交换后新解的哈希函数值
        int _new_hash_two = hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
        int _new_hash_three = hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
        _new_hash_one = (_new_hash_one + size_of_tabu_list) % size_of_tabu_list;                //xxf：解决bug2：防止出现负数和>size_of_tabu_list的数
        _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
        _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
        if (tabu_list_three[_new_hash_three]) {       //xxf：两个或三个哈希函数为了减小冲突碰撞；一个哈希函数容易不同解映射到同一个value值，且很容易出现邻域解都被禁忌情况，因为是从强化搜索策略开始
            if (tabu_list_two[_new_hash_two]) {
                if (tabu_list_one[_new_hash_one]);
                else break;
            }
            else break;
        }
        else break;
    }
    Distance temp_min = node_dis_sum[zero_toOne_node] - ins.dis_nodes(zero_toOne_node, one_toZero_node);      //记录当前一次动作的最大值和最小值;初始化为交换之后新的选中节点的距离之和
    Distance temp_max = temp_min;
    for (int k = 0; k < nb_sub_nodes; ++k) {
        if (k == i)continue;
        int node = select_nodes[k].first;
        Distance update_dis = node_dis_sum[node] - ins.dis_nodes(node, one_toZero_node) + ins.dis_nodes(node, zero_toOne_node);
        if (update_dis > temp_max)temp_max = update_dis;
        else if (temp_min > update_dis)temp_min = update_dis;
        else;
    }
    node_value[one_toZero_node] = 0;                       //更新当前解
    node_value[zero_toOne_node] = 1;
    cur_obj = temp_max - temp_min;
    max_select_node = temp_max;
    min_select_node = temp_min;

    tabu_list_one[_new_hash_one] = 1;         //更新三个禁忌列表
    tabu_list_two[_new_hash_two] = 1;
    tabu_list_three[_new_hash_three] = 1;
    hashone = _new_hash_one;
    hashtwo = _new_hash_two;
    hashthree = _new_hash_three;
        for (int i = 0; i < nb_nodes; ++i) {
            node_dis_sum[i] = node_dis_sum[i] - ins.dis_nodes(i, one_toZero_node) + ins.dis_nodes(i, zero_toOne_node);
        }
        no_select_nodes.clear();            //更新辅助结构select_nodes和no_select_nodes
        select_nodes.clear();
        Distance cur_temp_sum = (max_select_node + min_select_node) / 2;
        for (int i = 0; i < nb_nodes; ++i) {
            Distance temp = fabs(node_dis_sum[i] - cur_temp_sum);
            if (node_value[i])select_nodes.push_back(make_pair(i, temp));
            else no_select_nodes.push_back(make_pair(i, temp));
        }
    //sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //未选中的升序排列
    return true;
}

}