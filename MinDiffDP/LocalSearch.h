#pragma once
#ifndef MIN_DIFF_DP_LOCALSEARCH_H
#define MIN_DIFF_DP_LOACLSEARCH_H
#include "InputOutput.h"
//#include <utility>
constexpr auto EPS = 1e-6;

using namespace std;

namespace min_diff_dp{

class LocalSearch {

public:
    LocalSearch(const UMatrix &_matrix, double param, const Solution &_init_sol);
    Solution solve();

private:
    void init();         //初始化数据结构

private:
    void find_best_move(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj);   //局部搜索中在邻域结构中找到最好的交换动作
    void update_solu(const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj);    //做完邻域动作后更新当前解和历史最优解
    void update_auxiliary_structure(const pair<int, int> &_pair);        //更新邻域结构及两个复制结构并排序
    
private:
    const UMatrix &ins;      //算例矩阵

    List<int> node_value;    //当前解        ？（是用List保存，还是Set保存；更快）
    Distance cur_obj;        //当前目标
    List<int> local_best;     //局部最优解
    Distance local_best_obj;        //局部最优目标函数

    List<Distance> node_dis_sum;     //每个节点与选中集合中节点的距离和--邻域结构
    Distance max_select_node;                  //记录node_dis_sum中距离最大值的选中节点
    Distance min_select_node;                  //记录node_dis_sum中距离最小值的选中节点
    List<Pair<int, Distance>> no_select_nodes;         //未选中的节点排序--辅助结构,Dis-(max+min)/2
    List<Pair<int, Distance>> select_nodes;             //选中的节点排序--辅助结构
    //TODO：三个禁忌列表

    int nb_nodes;                   //图中节点数目
    int nb_sub_nodes;               //所选的节点数目

    int iter;                      //迭代次数
    int size_neighbor_struc;          //参与邻域动作的I0元素个数

    //参数相关
    long long max_time;             //迭代运行的最长时间=算例的节点数目
    double rate_of_sele_nodes;        //从I0中选择邻域动作的大小比例
    //TODO:禁忌步数参数；L大小；三个禁忌函数参数
};

}

#endif