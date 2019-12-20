#pragma once
#ifndef MIN_DIFF_DP_LOCALSEARCH_H
#define MIN_DIFF_DP_LOACLSEARCH_H
#include "InputOutput.h"

using namespace std;

namespace min_diff_dp{

class LocalSearch {

public:
    LocalSearch(const UMatrix &_matrix, double param, const Solution &_init_sol, const int _tabu_step,
        const int _size_of_tabu, const double _param1, const double _param2, const double _param3);
    Solution solve();

private:
    void init();         //��ʼ�����ݽṹ

private:
    void find_best_move(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three);   //�ֲ�������������ṹ���ҵ���õĽ�������
    bool update_solu(const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three);   //��������������µ�ǰ�����ʷ���Ž�

//private:
//    int hash_function_one();
//    int hash_function_two();
//    int hash_function_three();
//    int hash_function_temp_one(const vector<int>& temp);
//    int hash_function_temp_two(const vector<int>& temp);
//    int hash_function_temp_three(const vector<int>& temp);

private:
    const UMatrix &ins;      //��������

    List<int> node_value;    //��ǰ��        ��������List���棬����Set���棻���죩
    Distance cur_obj;        //��ǰĿ��
    List<int> local_best;     //�ֲ����Ž�
    Distance local_best_obj;        //�ֲ�����Ŀ�꺯��

    List<Distance> node_dis_sum;     //ÿ���ڵ���ѡ�м����нڵ�ľ����--����ṹ
    Distance max_select_node;                  //��¼node_dis_sum�о������ֵ��ѡ�нڵ�
    Distance min_select_node;                  //��¼node_dis_sum�о�����Сֵ��ѡ�нڵ�
                                 //������ʷ���Ž��������ݣ�ǿ���������ԣ�S<-��ʷ���Ž⣩�з���ʹ��
    List<Distance> best_solu_dis_sum;
    Distance best_max_select_node;
    Distance best_min_select_node;

    List<Pair<int, Distance>> no_select_nodes;         //δѡ�еĽڵ�����--�����ṹ,Dis-(max+min)/2
    List<Pair<int, Distance>> select_nodes;             //ѡ�еĽڵ�����--�����ṹ
    List<int> tabu_list_one;                   //���������б�
    List<int> tabu_list_two;
    List<int> tabu_list_three;
    List<int> hash_key_temp_one;          //�������нڵ��hash�м��ֵ(int)(floor(pow(i, hashFun_one_param)))
    List<int> hash_key_temp_two;
    List<int> hash_key_temp_three;

    int best_hashfun_one;      //�м�ֵ:��ʷ���Ž��������ϣ����ֵ
    int best_hashfun_two;
    int best_hashfun_three;             

    int nb_nodes;                   //ͼ�нڵ���Ŀ
    int nb_sub_nodes;               //��ѡ�Ľڵ���Ŀ

    int iter;                      //�ܵ�������
    int size_neighbor_struc;          //������������I0Ԫ�ظ���

    //�������
    long long max_time;             //�������е��ʱ��=�����Ľڵ���Ŀ
    double rate_of_sele_nodes;        //��I0��ѡ���������Ĵ�С����
    int tabu_step;                 //���ɲ���
    int size_of_tabu_list;        //L��С
    double hashFun_one_param;     //������ϣ�����Ĳ���
    double hashFun_two_param;
    double hashFun_three_param;
};

}

#endif