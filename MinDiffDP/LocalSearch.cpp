#include"LocalSearch.h"
#include"utility.h"

using namespace std;
using namespace xxf_utility;

namespace min_diff_dp {

LogSwitch logsw_local(1, 1, "PINCOL_PT");

bool judgeDistanceGreat(const Distance &a, const Distance &b) {    //如果a比b大，返回true
    return a > b + EPS;
}

bool compareByAscend(const pair<int, Distance> &a, const  pair<int, Distance> &b) {
    //return a.second + EPS < b.second;
    return judgeDistanceGreat(b.second, a.second);
}
bool compareByDescend(const pair<int, Distance> &a, const  pair<int, Distance> &b) {
    //return a.second > b.second + EPS;
    return judgeDistanceGreat(a.second, b.second);
}

LocalSearch::LocalSearch(const UMatrix &_matrix, double param, const Solution &_init_sol):
    ins(_matrix),max_select_node(-1),min_select_node(DISTANCE_MAX),nb_nodes(_matrix.setele_num()),
    nb_sub_nodes(_matrix.subsetele_num()),max_time(_matrix.setele_num()*1000),
    rate_of_sele_nodes(param),iter(0)
{
    node_value.resize(nb_nodes);
    node_value = _init_sol.node_values();
    local_best.resize(nb_nodes);
    local_best = _init_sol.node_values();
    init();
}

void LocalSearch::init() {              //xxf:done,right
    size_neighbor_struc = (int)((nb_nodes - nb_sub_nodes)*rate_of_sele_nodes);   //初始化List大小和参数
    node_dis_sum.resize(nb_nodes);
    select_nodes.reserve(nb_sub_nodes);
    no_select_nodes.reserve(nb_nodes - nb_sub_nodes);
    for (int i = 0; i < nb_nodes; ++i) {                     //初始化node_dis_sum和目标函数
        double dis = 0.0;
        for (int j = 0; j < nb_nodes; ++j) {
            if (node_value[j]) {
                dis += ins.dis_nodes(i, j);
            } 
        }
        node_dis_sum[i] = dis;
        if (node_value[i]) {
            if (max_select_node < dis)max_select_node= dis;
            if (min_select_node > dis)min_select_node = dis;
        }
    }
    cur_obj = max_select_node - min_select_node;
    local_best_obj = cur_obj;
    Distance temp_sum = (max_select_node + min_select_node) / 2;
    for (int i = 0; i < nb_nodes; ++i) {                     //初始化两个辅助结构no_select_nodes和select_nodes
        Distance temp = fabs(node_dis_sum[i] - temp_sum);
        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
        else no_select_nodes.push_back(make_pair(i, temp));
    }
    sort(select_nodes.begin(), select_nodes.end(), compareByDescend);     //选中的集合降序排列
    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //未选中的升序排列
}

Solution LocalSearch::solve() {
    Timer time(max_time);
    while (!time.isTimeOut()) {          //当时间未超时时，进行局部搜索
        //在邻域结构中找最好的动作 Xmm
        pair<int, int> swap_pair(-1, -1);             //保存最好的交换对;第一个I1-->I0，第二个I0-->I1
        pair<Distance, Distance> new_obj(DISTANCE_MAX, 0);         //保存对应的目标函数的最大距离和最小距离
        find_best_move(swap_pair, new_obj);
        //更新当前解和局部最优解，及目标函数 和 max min
        update_solu(swap_pair, new_obj);
        //更新三个数据结构并对两个辅助结构进行排序
        update_auxiliary_structure(swap_pair);
        iter++;
    }
    mylog << "\n总迭代步数：" << iter <<= logsw_local;
    return Solution(nb_nodes, nb_sub_nodes, local_best, local_best_obj);
}

void LocalSearch::find_best_move(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj) {   //xxf:done,right
    for(int i=0;i<nb_sub_nodes;++i){
        int one_toZero_node = select_nodes[i].first;
        for (int j = 0; j < size_neighbor_struc; ++j) {
            int zero_toOne_node = no_select_nodes[j].first;
            Distance temp_min = node_dis_sum[zero_toOne_node] - ins.dis_nodes(zero_toOne_node, one_toZero_node);      //记录当前一次动作的最大值和最小值;初始化为交换之后新的选中节点的距离之和
            Distance temp_max = temp_min;
            for (int k = 0; k < nb_sub_nodes; ++k) {
                if (k == i)continue;
                Distance update_dis = node_dis_sum[select_nodes[k].first] - ins.dis_nodes(select_nodes[k].first, one_toZero_node) + ins.dis_nodes(select_nodes[k].first, zero_toOne_node);  
                if (judgeDistanceGreat(update_dis, temp_max))temp_max = update_dis;
                if (judgeDistanceGreat(temp_min, update_dis))temp_min = update_dis;
            }
            if (judgeDistanceGreat(_new_obj.first - _new_obj.second, temp_max - temp_min)) {
                _new_obj.first = temp_max;
                _new_obj.second = temp_min;
                _pair.first = one_toZero_node;
                _pair.second = zero_toOne_node;
            }
        }
    }
}

void LocalSearch::update_solu(const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj) {  //xxf:done,right
    node_value[_pair.first] = 0;                       //更新当前解
    node_value[_pair.second] = 1;
    cur_obj = _new_obj.first - _new_obj.second;
    max_select_node = _new_obj.first;
    min_select_node = _new_obj.second;
    if (judgeDistanceGreat(local_best_obj, cur_obj)) {    //如果能改进历史最优解，则更新历史最优解
        local_best.assign(node_value.begin(), node_value.end());
        local_best_obj = cur_obj;
    }
}

void LocalSearch::update_auxiliary_structure(const pair<int, int> &_pair) {      //xxf:done,right
    no_select_nodes.clear();
    select_nodes.clear();
    select_nodes.reserve(nb_sub_nodes);
    no_select_nodes.reserve(nb_nodes - nb_sub_nodes);
    Distance temp_sum = (max_select_node + min_select_node) / 2;
    for (int i = 0; i < nb_nodes; ++i) {
        node_dis_sum[i] = node_dis_sum[i] - ins.dis_nodes(i, _pair.first) + ins.dis_nodes(i, _pair.second);
        Distance temp = fabs(node_dis_sum[i] - temp_sum);
        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
        else no_select_nodes.push_back(make_pair(i, temp));
    }
    sort(select_nodes.begin(), select_nodes.end(), compareByDescend);     //选中的集合降序排列
    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //未选中的升序排列
}

}