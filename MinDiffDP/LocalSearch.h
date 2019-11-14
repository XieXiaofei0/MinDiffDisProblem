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
    void init();         //��ʼ�����ݽṹ

private:
    void find_best_move(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj);   //�ֲ�������������ṹ���ҵ���õĽ�������
    void update_solu(const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj);    //��������������µ�ǰ�����ʷ���Ž�
    void update_auxiliary_structure(const pair<int, int> &_pair);        //��������ṹ���������ƽṹ������
    
private:
    const UMatrix &ins;      //��������

    List<int> node_value;    //��ǰ��        ��������List���棬����Set���棻���죩
    Distance cur_obj;        //��ǰĿ��
    List<int> local_best;     //�ֲ����Ž�
    Distance local_best_obj;        //�ֲ�����Ŀ�꺯��

    List<Distance> node_dis_sum;     //ÿ���ڵ���ѡ�м����нڵ�ľ����--����ṹ
    Distance max_select_node;                  //��¼node_dis_sum�о������ֵ��ѡ�нڵ�
    Distance min_select_node;                  //��¼node_dis_sum�о�����Сֵ��ѡ�нڵ�
    List<Pair<int, Distance>> no_select_nodes;         //δѡ�еĽڵ�����--�����ṹ,Dis-(max+min)/2
    List<Pair<int, Distance>> select_nodes;             //ѡ�еĽڵ�����--�����ṹ
    //TODO�����������б�

    int nb_nodes;                   //ͼ�нڵ���Ŀ
    int nb_sub_nodes;               //��ѡ�Ľڵ���Ŀ

    int iter;                      //��������
    int size_neighbor_struc;          //������������I0Ԫ�ظ���

    //�������
    long long max_time;             //�������е��ʱ��=�����Ľڵ���Ŀ
    double rate_of_sele_nodes;        //��I0��ѡ���������Ĵ�С����
    //TODO:���ɲ���������L��С���������ɺ�������
};

}

#endif