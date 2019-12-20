#include"LocalSearch.h"
#include"utility.h"
#include"InputOutput.h"

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
    iter(0),//max_time(100000*1000), 
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

void LocalSearch::init() {              //xxf:done,right-12.10;��ʱ��ϣ�����м��ֵ�޴�-12.18
    size_neighbor_struc = (int)((nb_nodes - nb_sub_nodes)*rate_of_sele_nodes);   //��ʼ��List��С�Ͳ���
    node_dis_sum.reserve(nb_nodes);
    best_solu_dis_sum.reserve(nb_nodes);
    select_nodes.reserve(nb_sub_nodes);
    no_select_nodes.reserve(nb_nodes - nb_sub_nodes);
    long long sum1 = 0, sum2 = 0, sum3 = 0;
    for (int i = 0; i < nb_nodes; ++i) {                         //��ʼ��node_dis_sum��Ŀ�꺯��
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
            sum1 += hash_key_temp_one[i];                      //���㵱ǰ��ʷ���Ž�������ϣ�����ļ�ֵ
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
    for (int i = 0; i < nb_nodes; ++i) {                                     //��ʼ�����������ṹno_select_nodes��select_nodes
        Distance temp = fabs(node_dis_sum[i] - temp_sum);
        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
        else no_select_nodes.push_back(make_pair(i, temp));
    }
    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //δѡ�е���������
}

Solution LocalSearch::solve() {
    Timer time(max_time);
    while (!time.isTimeOut()) {           //��ʱ��δ��ʱʱ�����оֲ�����
        node_value = local_best;               //ǿ���������ԣ�����ʷ���Ž���µ�ǰ��
        cur_obj = local_best_obj;
        //����ʷ���Ž�!=��ǰ��,��Ҫ���µ�ǰ��Ĳ������ݽṹnode_dis_sum��max_select_node�� min_select_node��select_nodes��no_select_nodes
        max_select_node = best_max_select_node;
        min_select_node = best_min_select_node;
        no_select_nodes.clear();            //���¸����ṹselect_nodes��no_select_nodes
        select_nodes.clear();
        Distance temp_sum = (best_max_select_node + best_min_select_node) / 2;
        for (int i = 0; i < nb_nodes; ++i) {
            node_dis_sum[i] = best_solu_dis_sum[i];
            Distance temp = fabs(best_solu_dis_sum[i] - temp_sum);
            if (node_value[i])select_nodes.push_back(make_pair(i, temp));
            else no_select_nodes.push_back(make_pair(i, temp));
        }
        sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //δѡ�е���������

        int count = 0;
        int _hashfun_one = best_hashfun_one;
        int _hashfun_two = best_hashfun_two;
        int _hashfun_three = best_hashfun_three;
        bool tabu_flag = false;                     //����Ƿ�����ⶼ�ڽ�����
        while (count <= tabu_step) {
            pair<int, int> swap_pair(-1, -1);              //����ṹ������õķǽ��ɽ�:����ǽ��ɵ���õĽ�����;��һ��I1-->I0���ڶ���I0-->I1
            pair<Distance, Distance> new_obj(DISTANCE_MAX, 0);         //�����Ӧ��Ŀ�꺯�������������С����
            //test:�ж�_hashfun_one�Ƿ��ǵ�ǰ���hashֵ
            //if (hash_function_one() != _hashfun_one) mylog << "whileѭ���У�hashֵ����\n" <<= logsw_local;
            //if (hash_function_two() != _hashfun_two) mylog << "whileѭ���У�hashֵ����\n" <<= logsw_local;
            //test end
            find_best_move(swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three);   //_hashfun_one�Ѿ��Ǹ��½�Ĺ�ϣ����ֵ
            //TODO:�ж��Ƿ���������ⶼ�ڽ�����
            if (swap_pair.first == -1) {       //�ж��Ƿ�����ⶼ�ڽ�����
                mylog << "����ⶼ�ڽ����У���ǰiter��" << iter <<= logsw_local;
                tabu_flag = true;
                break;
            } 
            if (update_solu(swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three))count = 0;    //���µ�ǰ�⡢��ʷ���Ž⡢count
            else count++; 
            iter++;
        }
        if (tabu_flag)break;
    }
    mylog << "\n�ܵ���������" << iter <<= logsw_local;
    return Solution(nb_nodes, nb_sub_nodes, local_best, local_best_obj);
}

void LocalSearch::find_best_move(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three){   //xxf:done,right--12.10
    int new_hashone = _hash_one;     //���浱ǰ��Ĺ�ϣ����ֵ
    int new_hashtwo = _hash_two;
    int new_hashthree = _hash_three;
    for (int i = 0; i < nb_sub_nodes; ++i)
    {
        int one_toZero_node = select_nodes[i].first;
        for (int j = 0; j < size_neighbor_struc; ++j) {
            int zero_toOne_node = no_select_nodes[j].first; 
            int _new_hash_one = new_hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //���㽻�����½�Ĺ�ϣ����ֵ
            int _new_hash_two = new_hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
            int _new_hash_three = new_hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
            _new_hash_one = _new_hash_one % size_of_tabu_list;
            _new_hash_two = _new_hash_two % size_of_tabu_list;
            _new_hash_three = _new_hash_three % size_of_tabu_list;
            //test:�ж������Ľ����Ƿ����
            //1.�жϵ�ǰ���hashֵ�Ƿ����
            //if (iter > 75) {
            //    if (hash_function_three() != new_hashthree) mylog << "find_best_move�У�hashֵ����\n" <<= logsw_local;
            //    list<int> node_temp(node_value.begin(), node_value.end());
            //    if (node_temp[zero_toone_node] == 0)node_temp[zero_toone_node] = 1;
            //    else mylog << "find_best_move�У��ҽ����ڵ����" <<= logsw_info;
            //    if (node_temp[one_tozero_node] == 1)node_temp[one_tozero_node] = 0;
            //    else mylog << "find_best_move�У��ҽ����ڵ����" <<= logsw_info;
            //    if (hash_function_temp_one(node_temp) != _new_hash_one) mylog << "find_best_move�У�������hashֵ����\n" <<= logsw_local;
            //    if (hash_function_temp_two(node_temp) != _new_hash_two) mylog << "find_best_move�У�������hashֵ����\n" <<= logsw_local;
            //    if (hash_function_temp_three(node_temp) != _new_hash_three) mylog << "find_best_move�У�������hashֵ����\n   " << hash_function_temp_three(node_temp) <<= logsw_local;
            //    //test end
            //}
            if (tabu_list_one[_new_hash_one]) {       //xxf��������˵�������ϣ�����ܹ������жϽ���״̬��ʱ�临�Ӷ�
                if (tabu_list_two[_new_hash_two])
                    if (tabu_list_three[_new_hash_three]) {
                        //test
                        //mylog << "iter��" << iter << ";  ���ɽ⣺" << _new_hash_one << " " << _new_hash_two << " " << _new_hash_three <<= logsw_local;
                        //test end
                        continue;
                    }
            }
            Distance temp_min = node_dis_sum[zero_toOne_node] - ins.dis_nodes(zero_toOne_node, one_toZero_node);      //��¼��ǰһ�ζ��������ֵ����Сֵ;��ʼ��Ϊ����֮���µ�ѡ�нڵ�ľ���֮��
            Distance temp_max = temp_min;
            for (int k = 0; k < nb_sub_nodes; ++k) {
                if (k == i)continue;
                int node = select_nodes[k].first;
                Distance update_dis = node_dis_sum[node] - ins.dis_nodes(node, one_toZero_node) + ins.dis_nodes(node, zero_toOne_node);
                if (update_dis > temp_max)temp_max = update_dis;
                else if (temp_min > update_dis)temp_min = update_dis;
                else;
            }
            if ((_new_obj.first - _new_obj.second) > (temp_max - temp_min)){   //����������
                _new_obj.first = temp_max;
                _new_obj.second = temp_min;
                _pair.first = one_toZero_node;
                _pair.second = zero_toOne_node;
                _hash_one = _new_hash_one;
                _hash_two = _new_hash_two;
                _hash_three = _new_hash_three;
            }
        }
    }
    //test
    //if (iter > 34822) {
        //mylog << "find_end: pair: " << _pair.first << " " << _pair.second <<= logsw_local;
        //if (node_value[_pair.first] != 1 || node_value[_pair.second] != 0)mylog << "�����ڵ����bug" <<= logsw_local;
        //mylog << "find_end: _new_obj: " << _new_obj.first << " " << _new_obj.second <<= logsw_local;
        //mylog << "find_end: hash: " << _hash_one << " " << _hash_two << " " << _hash_three <<= logsw_local;
        //if (_pair.first == -1)mylog << "find_move:����ⶼ�ڽ�����" <<= logsw_info;
    //}
}

bool LocalSearch::update_solu(const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three) {  //xxf:done,right-12.10
    bool flag = false;         //��ʾ�Ƿ������ʷ���Ž�
    node_value[_pair.first] = 0;                       //���µ�ǰ��
    node_value[_pair.second] = 1;
    cur_obj = _new_obj.first - _new_obj.second;
    max_select_node = _new_obj.first;
    min_select_node = _new_obj.second;
    //test 
    //if (hash_function_one() != _hash_one) mylog << "update:  hashֵ����" <<= logsw_local;
    //test end
    if (local_best_obj > cur_obj)       //������ʷ���Ž⣬������ʷ���Ž����ؽṹ
    {    
        local_best = node_value;  //����ܸĽ���ʷ���Ž⣬�������ʷ���Ž�,����true
        local_best_obj = cur_obj;
        best_hashfun_one = _hash_one;
        best_hashfun_two = _hash_two;
        best_hashfun_three = _hash_three;
        //test:TODO:���ٲ�֮���ܵ�������
        mylog << "\n��ǰΪ��" << local_best_obj << "     ������" << iter <<= logsw_local;
        flag = true;
    }
    //test:�ж�hashֵ�Ƿ��ǵ�ǰ���hashֵ
    //if (hash_function_three() != _hash_three) mylog << "update�У�hashֵ����\n" <<= logsw_local;
    //mylog << "��Ľ���ֵ�� " << tabu_list_one[_hash_one] << " " << tabu_list_two[_hash_two] << " " << tabu_list_three[_hash_three] <<= logsw_info;
    //test end
    tabu_list_one[_hash_one] = 1;         //�������������б�
    tabu_list_two[_hash_two] = 1;
    tabu_list_three[_hash_three] = 1;
    //mylog << "��Ľ���ֵ -���º� " << tabu_list_one[_hash_one] << " " << tabu_list_two[_hash_two] << " " << tabu_list_three[_hash_three] <<= logsw_info;
    //test
    //mylog << "��ǰ���ɵĽ�Ϊ��" << _hash_one << "  " << _hash_two << "  " << _hash_three <<= logsw_info;
    //test end
    no_select_nodes.clear();           //���¸����ṹselect_nodes��no_select_nodes
    select_nodes.clear();
    Distance temp_sum = (max_select_node + min_select_node) / 2.0;
    if (flag) {
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
    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //δѡ�е���������
    return flag;
}

//void LocalSearch::update_auxiliary_structure(const pair<int, int> &_pair, const int &_hash_one, const  int &_hash_two, const int &_hash_three) {      //xxf:done,right
//    //
//    tabu_list_one[_hash_one] = 1;         //�������������б�
//    tabu_list_two[_hash_two] = 1;
//    tabu_list_three[_hash_three] = 1;
//    no_select_nodes.clear();           //���¸����ṹselect_nodes��no_select_nodes
//    select_nodes.clear();
//    select_nodes.reserve(nb_sub_nodes);
//    no_select_nodes.reserve(nb_nodes - nb_sub_nodes);
//    Distance temp_sum = (max_select_node + min_select_node) / 2;
//    for (int i = 0; i < nb_nodes; ++i) {
//        //node_dis_sum[i] = node_dis_sum[i] - _matrix[i][_pair.first] + _matrix[i][_pair.second];
//        node_dis_sum[i] = node_dis_sum[i] - ins.dis_nodes(i, _pair.first) + ins.dis_nodes(i, _pair.second);
//        Distance temp = fabs(node_dis_sum[i] - temp_sum);
//        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
//        else no_select_nodes.push_back(make_pair(i, temp));
//    }
//    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //δѡ�е���������
//}

//int LocalSearch::hash_function_one() {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (node_value[i]) {
//            sum += (int)(floor(pow(i, hashFun_one_param)));
//        }
//    }
//    return sum % size_of_tabu_list;
//}
//
//int LocalSearch::hash_function_two() {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (node_value[i]) {
//            sum += hash_key_temp_two[i];
//        }
//    }
//    return sum % size_of_tabu_list;
//}
//
//int LocalSearch::hash_function_three() {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (node_value[i]) {
//            sum += hash_key_temp_three[i];
//        }
//    }
//    return sum % size_of_tabu_list;
//}
//
//int LocalSearch::hash_function_temp_one(const vector<int>& temp) {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (temp[i]) {
//            sum += hash_key_temp_one[i];
//        }
//    }
//    return sum % size_of_tabu_list;
//}
//
//int LocalSearch::hash_function_temp_two(const vector<int>& temp) {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (temp[i]) {
//            sum += hash_key_temp_two[i];
//        }
//    }
//    return sum % size_of_tabu_list;
//}
//
//int LocalSearch::hash_function_temp_three(const vector<int>& temp) {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (temp[i]) {
//            sum += hash_key_temp_three[i];
//        }
//    }
//    return sum % size_of_tabu_list;
//}

}