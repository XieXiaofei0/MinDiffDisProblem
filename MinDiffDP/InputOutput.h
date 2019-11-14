#pragma once
#ifndef MIN_DIFF_DP_INPUTOUTPUT_H
#define MIN_DIFF_DP_INPUTOUTPUT_H
#include "common.h"

namespace min_diff_dp{

class UMatrix {
public:
    UMatrix() :nb_set_ele(0), nb_subset_ele(0) {}
    UMatrix(const String &path);

    size_t setele_num() const { return nb_set_ele; }
    size_t subsetele_num() const { return nb_subset_ele; }
    const Distance& dis_nodes(int i, int j) const { return matrix[i][j]; }


private:
    size_t nb_set_ele;              //�ܽڵ���
    size_t nb_subset_ele;           //��ѡ�Ľڵ�����
    List<List<Distance>> matrix;     //�������
};

class Solution {
public:
    Solution(size_t nb_nodes = 0, size_t nb_sel_nodes = 0) :
        nb_set_ele(nb_nodes), nb_subset_ele(nb_sel_nodes), object(DISTANCE_MAX) {
        ele_value.resize(nb_nodes, 0);
    }
    Solution(Solution &other) :nb_set_ele(other.nb_set_ele), nb_subset_ele(other.nb_subset_ele),
        ele_value(other.ele_value), object(other.object) {
    }
    Solution(Solution &&other) :nb_set_ele(other.nb_set_ele), nb_subset_ele(other.nb_subset_ele),
        ele_value(std::move(other.ele_value)), object(other.object) {
    }
    Solution& operator=(Solution &rhs) {
        nb_set_ele = rhs.nb_set_ele;
        nb_subset_ele = rhs.nb_subset_ele;
        ele_value = rhs.ele_value;
        object = rhs.object;
        return *this;
    }
    Solution& operator=(Solution &&rhs) {
        nb_set_ele = rhs.nb_set_ele;
        nb_subset_ele = rhs.nb_subset_ele;
        ele_value = std::move(rhs.ele_value);
        object = rhs.object;
        return *this;
    }
    Solution(int nb_nodes, int nb_sub_nodes, const List<int> &_sol, Distance _best_obj) {
        nb_set_ele = nb_nodes;
        nb_subset_ele = nb_sub_nodes;
        ele_value.resize(nb_nodes, 0);
        ele_value.assign(_sol.begin(), _sol.end());
        object = _best_obj;
    }

    int& operator[] (int node) { return ele_value[node]; }
    const int get_nb_nodes() const { return nb_set_ele; }
    const int get_sub_nb_nodes() const { return nb_subset_ele; }
    const List<int>& node_values() const { return ele_value; }
    const Distance get_object() const { return object; }

    void randomInit();
    void print();

private:
    int nb_set_ele;              //�ܽڵ���Ŀ
    int nb_subset_ele;           //��ѡ�Ľڵ���Ŀ
    List<int> ele_value;         //һάvector��ʶ�ڵ��Ƿ�ѡ��
    //Set<int> select_nodes;        //ѡ�еĽڵ㼯�ϡ�����������List���滹��Set���棩
    Distance object;               //Ŀ�꺯��ֵ
};

}

#endif