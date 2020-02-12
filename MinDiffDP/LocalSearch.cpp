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
    iter(0),//max_time(1000000*1000),                    //xxf:�޸����������ʱ��ms
    max_time(_matrix.setele_num() *2* 1000),
    //max_time(_matrix.setele_num()*1000),
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
    //��Ľ���
    tabu_nodes.resize(nb_nodes, 0);
    init();
}

void LocalSearch::init() {              //xxf:done,right-12.10;��ʱ��ϣ�����м��ֵ�޴�-12.18
    size_neighbor_struc = (int)((nb_nodes - nb_sub_nodes)*rate_of_sele_nodes);   //��ʼ��List��С�Ͳ���
    node_dis_sum.reserve(nb_nodes);
    best_solu_dis_sum.reserve(nb_nodes);
    select_nodes.reserve(nb_sub_nodes);
    no_select_nodes.reserve(nb_nodes - nb_sub_nodes);
    long long sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
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
    //���ý��ɲ���
    int tabu_length = 200;
    //end
    Timer time(max_time);
    int step_length = 0;           //��¼���ٲ�֮�ڲ��ܸĽ���ʷ���Ž�
    bool tabu_flag = false;          //�ж��Ƿ�����ⶼ�ڽ�����
    //test
    //double k = 0.07;              //��ֵ̬������ṹ�ĳ���
    //test end
    while (!time.isTimeOut()) {           //��ʱ��δ��ʱʱ�����оֲ����� 
        int _hashfun_one = best_hashfun_one;
        int _hashfun_two = best_hashfun_two;
        int _hashfun_three = best_hashfun_three;

        //����Ŷ�
        //if (tabu_flag) {                 //�������ⶼ�ڽ����л���һ�������ڸĽ�������ʷ���Ž⣬��ʱ��Ϊ5000
        //    //test
        //    mylog << "��ǰ��ʷ���Ž⣺ " << local_best_obj <<= logsw_info;
        //    mylog << "��ǰ�⣺ " << cur_obj <<= logsw_info;
        //    //test end
        //    stochastic_perturbation(_hashfun_one, _hashfun_two, _hashfun_three, max_select_node, min_select_node);
        //    tabu_flag = false;          //��־λ��λ
        //    step_length = 0;
        //    //test
        //    mylog << "����⣺ " << cur_obj <<= logsw_info;
        ////mylog << "�����Ϊ�� ";
        ////for (int i = 0; i < nb_nodes; ++i) {
        ////    if (node_value[i])mylog << i << "  ";
        ////}
        //}
        //else {
        //    node_value = local_best;               //ǿ���������ԣ�����ʷ���Ž���µ�ǰ��
        //    cur_obj = local_best_obj;
        //    //����ʷ���Ž�!=��ǰ��,��Ҫ���µ�ǰ��Ĳ������ݽṹnode_dis_sum��max_select_node�� min_select_node��select_nodes��no_select_nodes
        //    max_select_node = best_max_select_node;                   //xxf:�����bug1��֮ǰ�滻��ǰ��ʱ��δ����������ݽṹnode_dis_sum��select_nodes��no_select_nodes
        //    min_select_node = best_min_select_node;
        //    no_select_nodes.clear();            //���¸����ṹselect_nodes��no_select_nodes
        //    select_nodes.clear();
        //    Distance temp_sum = (best_max_select_node + best_min_select_node) / 2;
        //    for (int i = 0; i < nb_nodes; ++i) {
        //        node_dis_sum[i] = best_solu_dis_sum[i];
        //        Distance temp = fabs(best_solu_dis_sum[i] - temp_sum);
        //        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
        //        else no_select_nodes.push_back(make_pair(i, temp));
        //    }
        //}
                //�����ѡ���������Ŷ�
        //if (tabu_flag) {
        //    for (int i = 0; i < 3; ++i) {
        //        stochastic_perturbation(_hashfun_one, _hashfun_two, _hashfun_three, max_select_node, min_select_node);
        //    }
        //    no_select_nodes.clear();           //���¸����ṹselect_nodes��no_select_nodes
        //    select_nodes.clear();
        //    Distance temp_sum = (max_select_node + min_select_node) / 2.0;
        //    for (int i = 0; i < nb_nodes; ++i) {
        //        Distance temp = fabs(node_dis_sum[i] - temp_sum);
        //        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
        //        else no_select_nodes.push_back(make_pair(i, temp));
        //    }
        //    sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //δѡ�е���������
        //    tabu_flag = false;
        //    mylog << "��ǰ��ʷ���Ž⣺ " << local_best_obj <<= logsw_info;
        //    mylog << "��ǰ�⣺ " << cur_obj <<= logsw_info;
        //}

        node_value = local_best;               //ǿ���������ԣ�����ʷ���Ž���µ�ǰ��
        cur_obj = local_best_obj;
        //����ʷ���Ž�!=��ǰ��,��Ҫ���µ�ǰ��Ĳ������ݽṹnode_dis_sum��max_select_node�� min_select_node��select_nodes��no_select_nodes
        max_select_node = best_max_select_node;                   //xxf:�����bug1��֮ǰ�滻��ǰ��ʱ��δ����������ݽṹnode_dis_sum��select_nodes��no_select_nodes
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

        while (count <= tabu_step) {
            pair<int, int> swap_pair(-1, -1);              //����ṹ������õķǽ��ɽ�:����ǽ��ɵ���õĽ�����;��һ��I1-->I0���ڶ���I0-->I1
            pair<Distance, Distance> new_obj(DISTANCE_MAX, 0);         //�����Ӧ��Ŀ�꺯�������������С����
            find_best_move(swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three);   //_hashfun_one�Ѿ��Ǹ��½�Ĺ�ϣ����ֵ
            //TODO:�ж��Ƿ���������ⶼ�ڽ�����
            if (swap_pair.first == -1) {       //�ж��Ƿ�����ⶼ�ڽ�����
                mylog << "����ⶼ�ڽ����У���ǰiter��" << iter <<= logsw_local;
                tabu_flag = true;
                break;
            } 
            if (update_solu(swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three, step_length))count = 0;   //���µ�ǰ�⡢��ʷ���Ž⡢count
            else count++; 
            iter++;
            //if (step_length == 5000) {
            //    //tabu_flag = true;
            //    size_neighbor_struc += (int)((nb_nodes - nb_sub_nodes)*k);
            //    k -= 0.01;
            //    if (k < 0)k = 0;
            //    //break;
            //}
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
    int num = 1;                             //�����������ȵĸ���
    for (int i = 0; i < nb_sub_nodes; ++i)
    {
        int one_toZero_node = select_nodes[i].first;
        for (int j = 0; j < size_neighbor_struc; ++j) {
            int zero_toOne_node = no_select_nodes[j].first; 
            int _new_hash_one = new_hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //���㽻�����½�Ĺ�ϣ����ֵ
            int _new_hash_two = new_hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
            int _new_hash_three = new_hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
            _new_hash_one = (_new_hash_one+ size_of_tabu_list) % size_of_tabu_list;                //xxf�����bug2����ֹ���ָ�����>size_of_tabu_list����
            _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
            _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
            if (tabu_list_three[_new_hash_three]) {       //xxf��������������ϣ����Ϊ�˼�С��ͻ��ײ��һ����ϣ�������ײ�ͬ��ӳ�䵽ͬһ��valueֵ���Һ����׳�������ⶼ�������������Ϊ�Ǵ�ǿ���������Կ�ʼ
                if (tabu_list_two[_new_hash_two])
                    if (tabu_list_one[_new_hash_one]) {
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
            if ((_new_obj.first - _new_obj.second) > (temp_max - temp_min)) {   //����������
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
    bool flag = false;         //��ʾ�Ƿ������ʷ���Ž�
    node_value[_pair.first] = 0;                       //���µ�ǰ��
    node_value[_pair.second] = 1;
    cur_obj = _new_obj.first - _new_obj.second;
    max_select_node = _new_obj.first;
    min_select_node = _new_obj.second;
    if (local_best_obj > cur_obj)       //������ʷ���Ž⣬������ʷ���Ž����ؽṹ
    {
        local_best = node_value;  //����ܸĽ���ʷ���Ž⣬�������ʷ���Ž�,����true
        local_best_obj = cur_obj;
        best_hashfun_one = _hash_one;
        best_hashfun_two = _hash_two;
        best_hashfun_three = _hash_three;
        //test:TODO:���ٲ�֮���ܵ�������
        mylog << "\n��ǰΪ��" << local_best_obj << "     ������" << iter <<= logsw_local;
        //test end
        flag = true;
        step = 0;
    }
    step++;
    tabu_list_one[_hash_one] = 1;         //�������������б�
    tabu_list_two[_hash_two] = 1;
    tabu_list_three[_hash_three] = 1;
    no_select_nodes.clear();           //���¸����ṹselect_nodes��no_select_nodes
    select_nodes.clear();
    Distance temp_sum = (max_select_node + min_select_node) / 2.0;
    if (flag) {           //����������ʷ���Ž⣬�򱣴���ʷ���Ž��������ݽṹbest_solu_dis_sum
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

Solution LocalSearch::solve_tabu() {
    //���ý��ɲ���
    int tabu_length = 40;
    int step_length = 0;

    Timer time(max_time);
    bool tabu_flag = false;          //�ж��Ƿ�����ⶼ�ڽ�����
    while (!time.isTimeOut()) {           //��ʱ��δ��ʱʱ�����оֲ����� 
        int _hashfun_one = best_hashfun_one;
        int _hashfun_two = best_hashfun_two;
        int _hashfun_three = best_hashfun_three;

        node_value = local_best;               //ǿ���������ԣ�����ʷ���Ž���µ�ǰ��
        cur_obj = local_best_obj;
        //����ʷ���Ž�!=��ǰ��,��Ҫ���µ�ǰ��Ĳ������ݽṹnode_dis_sum��max_select_node�� min_select_node��select_nodes��no_select_nodes
        max_select_node = best_max_select_node;                   //xxf:�����bug1��֮ǰ�滻��ǰ��ʱ��δ����������ݽṹnode_dis_sum��select_nodes��no_select_nodes
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

        while (count <= tabu_step) {
            pair<int, int> swap_pair(-1, -1);              //����ṹ������õķǽ��ɽ�:����ǽ��ɵ���õĽ�����;��һ��I1-->I0���ڶ���I0-->I1
            pair<Distance, Distance> new_obj(DISTANCE_MAX, 0);         //�����Ӧ��Ŀ�꺯�������������С����
            bool flag = find_best_move_tabu(swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three, iter);   //_hashfun_one�Ѿ��Ǹ��½�Ĺ�ϣ����ֵ
            //TODO:�ж��Ƿ���������ⶼ�ڽ�����
            if (swap_pair.first == -1) {       //�ж��Ƿ�����ⶼ�ڽ�����
                mylog << "����ⶼ�ڽ����У���ǰiter��" << iter <<= logsw_local;
                tabu_flag = true;
                break;
            }
            if (update_solu_tabu(flag, swap_pair, new_obj, _hashfun_one, _hashfun_two, _hashfun_three, step_length, tabu_length))count = 0;   //���µ�ǰ�⡢��ʷ���Ž⡢count
            else count++;
            iter++;
            //if (local_best_obj == 39) {
            //    tabu_length = 30;
            //}
        }
        if (tabu_flag)break;
    }
    mylog << "\n�ܵ���������" << iter <<= logsw_local;
    return Solution(nb_nodes, nb_sub_nodes, local_best, local_best_obj);
}

bool LocalSearch::find_best_move_tabu(pair<int, int> &_pair, pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, const int &iter) {   //xxf:done,right--12.10
    int new_hashone = _hash_one;     //���浱ǰ��Ĺ�ϣ����ֵ
    int new_hashtwo = _hash_two;
    int new_hashthree = _hash_three;

    pair<int, int> pair_tabu(-1, -1);
    pair<int, int> pair_notabu(-1,-1);
    pair<Distance, Distance> _new_obj_tabu(DISTANCE_MAX, 0);
    pair<Distance, Distance> _new_obj_notabu(DISTANCE_MAX, 0);
    int hash_one_tabu;
    int hash_two_tabu;
    int hash_three_tabu;
    int hash_one_notabu;
    int hash_two_notabu;
    int hash_three_notabu;

    for (int i = 0; i < nb_sub_nodes; ++i)
    {
        int one_toZero_node = select_nodes[i].first;
        for (int j = 0; j < size_neighbor_struc; ++j) {
            int zero_toOne_node = no_select_nodes[j].first;
            int _new_hash_one = new_hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //���㽻�����½�Ĺ�ϣ����ֵ
            int _new_hash_two = new_hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
            int _new_hash_three = new_hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
            _new_hash_one = (_new_hash_one + size_of_tabu_list) % size_of_tabu_list;                //xxf�����bug2����ֹ���ָ�����>size_of_tabu_list����
            _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
            _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
            if (tabu_list_three[_new_hash_three]) {       //xxf��������������ϣ����Ϊ�˼�С��ͻ��ײ��һ����ϣ�������ײ�ͬ��ӳ�䵽ͬһ��valueֵ���Һ����׳�������ⶼ�������������Ϊ�Ǵ�ǿ���������Կ�ʼ
                if (tabu_list_two[_new_hash_two])
                    if (tabu_list_one[_new_hash_one]) {
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
            //�ֱ��ҳ������е����Ž�ͷǽ����е����Ž�
            if (tabu_nodes[zero_toOne_node] > iter && tabu_nodes[one_toZero_node] > iter) {   //����״̬
                if ((_new_obj_tabu.first - _new_obj_tabu.second) > (temp_max - temp_min)) {   //����������
                    pair_tabu.first = one_toZero_node;
                    pair_tabu.second = zero_toOne_node;
                    _new_obj_tabu.first = temp_max;
                    _new_obj_tabu.second = temp_min;
                    hash_one_tabu = _new_hash_one;
                    hash_two_tabu = _new_hash_two;
                    hash_three_tabu = _new_hash_three;
                }
            }
            else {
                if ((_new_obj_notabu.first - _new_obj_notabu.second) > (temp_max - temp_min)) {   //����������
                    pair_notabu.first = one_toZero_node;
                    pair_notabu.second = zero_toOne_node;
                    _new_obj_notabu.first = temp_max;
                    _new_obj_notabu.second = temp_min;
                    hash_one_notabu = _new_hash_one;
                    hash_two_notabu = _new_hash_two;
                    hash_three_notabu = _new_hash_three;
                }
            }
        }

    }
    ////test
    //mylog << "test   " <<= logsw_info;
    //mylog << "������ɽ⣺  " << pair_tabu.first << "  " << pair_tabu.second << "  " << _new_obj_tabu.first - _new_obj_tabu.second <<= logsw_info;
    //mylog << "����ǽ��ɽ�:    " << pair_notabu.first << "  " << pair_notabu.second << "  " << _new_obj_notabu.first - _new_obj_notabu.second <<= logsw_info;
    ////test end
    if (pair_tabu.first == -1) {
        _pair = pair_notabu;
        _new_obj = _new_obj_notabu;
        _hash_one = hash_one_notabu;
        _hash_two = hash_two_notabu;
        _hash_three = hash_three_notabu;
        return false;
    }
    else if (pair_notabu.first == -1) {
        _pair = pair_tabu;
        _new_obj = _new_obj_tabu;
        _hash_one = hash_one_tabu;
        _hash_two = hash_two_tabu;
        _hash_three = hash_three_tabu;
        return true;
    }
    else {
        Distance tabu_obj = _new_obj_tabu.first - _new_obj_tabu.second;
        if ((tabu_obj < _new_obj_notabu.first - _new_obj_notabu.second) && (tabu_obj < local_best_obj)) {  //������ɽ����ڷǽ��ɽ����ܹ��Ľ���ʷ���Ž�
            mylog << "                   ���ɽ�Ľ���ʷ���Ž�  " <<= logsw_info;
            _pair = pair_tabu;
            _new_obj = _new_obj_tabu;
            _hash_one = hash_one_tabu;
            _hash_two = hash_two_tabu;
            _hash_three = hash_three_tabu;
            return true;
        }
        else
        {
            _pair = pair_notabu;
            _new_obj = _new_obj_notabu;
            _hash_one = hash_one_notabu;
            _hash_two = hash_two_notabu;
            _hash_three = hash_three_notabu;
            return false;
        }
    }

}

bool LocalSearch::update_solu_tabu(bool tabu_flag, const pair<int, int> &_pair, const pair<Distance, Distance> &_new_obj, int &_hash_one, int &_hash_two, int &_hash_three, int &step, int tabu_length) {  //xxf:done,right-12.10
    bool flag = false;         //��ʾ�Ƿ������ʷ���Ž�
    node_value[_pair.first] = 0;                       //���µ�ǰ��
    node_value[_pair.second] = 1;
    cur_obj = _new_obj.first - _new_obj.second;
    max_select_node = _new_obj.first;
    min_select_node = _new_obj.second;
    if (local_best_obj > cur_obj)       //������ʷ���Ž⣬������ʷ���Ž����ؽṹ
    {
        local_best = node_value;  //����ܸĽ���ʷ���Ž⣬�������ʷ���Ž�,����true
        local_best_obj = cur_obj;
        best_hashfun_one = _hash_one;
        best_hashfun_two = _hash_two;
        best_hashfun_three = _hash_three;
        //test:TODO:���ٲ�֮���ܵ�������
        mylog << "\n��ǰΪ��" << local_best_obj << "     ������" << iter <<= logsw_local;
        //test end
        flag = true;
        step = 0;
    }
    step++;
    tabu_list_one[_hash_one] = 1;         //�������������б�
    tabu_list_two[_hash_two] = 1;
    tabu_list_three[_hash_three] = 1;
    //if (tabu_flag) {
        //test
        //mylog << "             ���             " <<= logsw_info;
        //test end
        //���ɽ��ܹ��Ľ���ʷ���Ž⣬ִ�н������
        //if (tabu_nodes[_pair.first] > iter)tabu_nodes[_pair.first] = iter;
        //else tabu_nodes[_pair.first] = iter + tabu_length + rand() % 10;
        //if (tabu_nodes[_pair.second] > iter)tabu_nodes[_pair.second] = iter;
        //else tabu_nodes[_pair.second] = iter + tabu_length + rand() % 10;
    //}
    //else {
    tabu_nodes[_pair.first] = iter + tabu_length + rand() % 5;
    tabu_nodes[_pair.second] = iter + tabu_length + rand() % 5;
    //}
    no_select_nodes.clear();           //���¸����ṹselect_nodes��no_select_nodes
    select_nodes.clear();
    Distance temp_sum = (max_select_node + min_select_node) / 2.0;
    if (flag) {           //����������ʷ���Ž⣬�򱣴���ʷ���Ž��������ݽṹbest_solu_dis_sum
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

//��һ������Ŷ�����:�ӵ�ǰ��ʷ���Ž�ĺ�ѡ���������
//bool LocalSearch::stochastic_perturbation(int &hashone, int &hashtwo, int &hashthree, Distance &max, Distance &min) {
//    int one_toZero_node = 0;
//    int zero_toOne_node = 0;
//    int i = 0;
//    int _new_hash_one = 0;
//    int _new_hash_two = 0;
//    int _new_hash_three = 0;
//    while (true) {
//        srand((unsigned)time(NULL));
//        i = rand() % nb_sub_nodes;
//        int j = rand() % size_neighbor_struc;
//        one_toZero_node = select_nodes[i].first;
//        zero_toOne_node = no_select_nodes[j].first;
//        int _new_hash_one = hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //���㽻�����½�Ĺ�ϣ����ֵ
//        int _new_hash_two = hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
//        int _new_hash_three = hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
//        _new_hash_one = (_new_hash_one + size_of_tabu_list) % size_of_tabu_list;                //xxf�����bug2����ֹ���ָ�����>size_of_tabu_list����
//        _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
//        _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
//        if (tabu_list_three[_new_hash_three]) {       //xxf��������������ϣ����Ϊ�˼�С��ͻ��ײ��һ����ϣ�������ײ�ͬ��ӳ�䵽ͬһ��valueֵ���Һ����׳�������ⶼ�������������Ϊ�Ǵ�ǿ���������Կ�ʼ
//            if (tabu_list_two[_new_hash_two]) {
//                if (tabu_list_one[_new_hash_one]);
//                else break;
//            }
//            else break;
//        }
//        else break;
//
//    }
//    Distance temp_min = node_dis_sum[zero_toOne_node] - ins.dis_nodes(zero_toOne_node, one_toZero_node);      //��¼��ǰһ�ζ��������ֵ����Сֵ;��ʼ��Ϊ����֮���µ�ѡ�нڵ�ľ���֮��
//    Distance temp_max = temp_min;
//    for (int k = 0; k < nb_sub_nodes; ++k) {
//        if (k == i)continue;
//        int node = select_nodes[k].first;
//        Distance update_dis = node_dis_sum[node] - ins.dis_nodes(node, one_toZero_node) + ins.dis_nodes(node, zero_toOne_node);
//        if (update_dis > temp_max)temp_max = update_dis;
//        else if (temp_min > update_dis)temp_min = update_dis;
//        else;
//    }
//    node_value[one_toZero_node] = 0;                       //���µ�ǰ��
//    node_value[zero_toOne_node] = 1;
//    cur_obj = temp_max - temp_min;
//    max_select_node = temp_max;
//    min_select_node = temp_min;
//
//
//    tabu_list_one[_new_hash_one] = 1;         //�������������б�
//    tabu_list_two[_new_hash_two] = 1;
//    tabu_list_three[_new_hash_three] = 1;
//    hashone = _new_hash_one;
//    hashtwo = _new_hash_two;
//    hashthree = _new_hash_three;
//        for (int i = 0; i < nb_nodes; ++i) {
//            node_dis_sum[i] = node_dis_sum[i] - ins.dis_nodes(i, one_toZero_node) + ins.dis_nodes(i, zero_toOne_node);
//        }
//    //sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //δѡ�е���������
//    return true;
//}

//����ʷ���Ž��������Ŷ�
//bool LocalSearch::stochastic_perturbation(int &hashone, int &hashtwo, int &hashthree, Distance &max, Distance &min) {
//    hashone = 0;
//    hashtwo = 0;
//    hashthree = 0;
//    int num_of_random_select = nb_sub_nodes / 5;       //����Ŷ��������ѡ���Ľڵ�ĸ���
//    vector<int> node_temp(nb_nodes, 0);
//    vector<int> selectNodes;           //��������Ŷ�ѡ���Ľڵ�
//    selectNodes.reserve(nb_sub_nodes);
//    int num_select = 0, i = 0;      //�ڵ�������
//    srand((unsigned)time(NULL));
//    while (num_select < num_of_random_select) {        //�����ѡ��һ��ڵ�
//        while (true) {
//            i = rand() % nb_nodes;
//            if (local_best[i] == 1 && node_temp[i] == 0)break;
//        }
//        node_temp[i] = 1;
//        selectNodes.push_back(i);
//        hashone += hash_key_temp_one[i];      //ͬʱ�����µĽ��hashֵ
//        hashtwo += hash_key_temp_two[i];
//        hashthree += hash_key_temp_three[i];
//        num_select++;
//    }
//
//    //���ѡ��һ���ֽڵ�
//    //int num = 0;
//    //srand((unsigned)time(NULL));
//    //while (num < nb_sub_nodes - num_of_random_select) {        //�����ѡ��һ��ڵ�
//    //    while (true) {
//    //        i = rand() % nb_nodes;
//    //        if (local_best[i] == 0 && node_temp[i] == 0)break;
//    //    }
//    //    node_temp[i] = 1;
//    //    selectNodes.push_back(i);
//    //    hashone += hash_key_temp_one[i];      //ͬʱ�����µĽ��hashֵ
//    //    hashtwo += hash_key_temp_two[i];
//    //    hashthree += hash_key_temp_three[i];
//    //    num++;
//    //    //test
//    //    //new_nodes.push_back(i);
//    //    //test end
//    //}
//
//    //̰����ѡ��һ���ֽڵ�
//    Distance temp_sum = (best_max_select_node + best_min_select_node) / 2;     //̰����ѡ��һ��ڵ�
//    priority_queue<pair<Distance, int>, vector<pair<Distance, int>>, greater<pair<Distance, int>> > q;
//    for (int m = 0; m < nb_nodes; ++m) {
//        if (local_best[m] == 0) {
//            q.push(make_pair(fabs(best_solu_dis_sum[m] - temp_sum), m));
//        }
//    }
//    for (int m = 0; m < nb_sub_nodes - num_of_random_select; ++m) {
//        pair <Distance, int> first = q.top();
//        node_temp[first.second] = 1;
//        selectNodes.push_back(first.second);
//        hashone += hash_key_temp_one[first.second];      //ͬʱ�����µĽ��hashֵ
//        hashtwo += hash_key_temp_two[first.second];
//        hashthree += hash_key_temp_three[first.second];
//        q.pop();
//    }
//    hashone = hashone % size_of_tabu_list;
//    hashtwo = hashtwo % size_of_tabu_list;
//    hashthree = hashthree % size_of_tabu_list;
//    tabu_list_one[hashone] = 1;
//    tabu_list_two[hashtwo] = 1;
//    tabu_list_three[hashthree] = 1;
//    node_value = node_temp;
//    max = 0, min = DISTANCE_MAX;
//    for (int l = 0; l < nb_nodes; ++l) {       //����node_dis_num
//        node_dis_sum[l] = 0;
//        for (int m = 0; m < nb_sub_nodes; ++m) {
//            node_dis_sum[l] += ins.dis_nodes(l, selectNodes[m]);
//        }
//        if (node_value[l]) {
//            if (max < node_dis_sum[l])max = node_dis_sum[l];
//            if (min > node_dis_sum[l])min = node_dis_sum[l];
//        }
//    }
//    cur_obj = max - min;
//    no_select_nodes.clear();            //���¸����ṹselect_nodes��no_select_nodes
//    select_nodes.clear();
//    Distance cur_temp_sum = (max + min) / 2;
//    for (int i = 0; i < nb_nodes; ++i) {
//        Distance temp = fabs(node_dis_sum[i] - cur_temp_sum);
//        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
//        else no_select_nodes.push_back(make_pair(i, temp));
//    }
//    return true;
//}

//�ӵ�ǰ���������Ŷ�
//bool LocalSearch::stochastic_perturbation(int &hashone, int &hashtwo, int &hashthree, Distance &max, Distance &min) {
//    hashone = 0;
//    hashtwo = 0;
//    hashthree = 0;
//    int num_of_random_select = nb_sub_nodes / 4 * 3;       //����Ŷ��������ѡ���Ľڵ�ĸ���
//    vector<int> node_temp(nb_nodes, 0);
//    vector<int> selectNodes;           //��������Ŷ�ѡ���Ľڵ�
//    selectNodes.reserve(nb_sub_nodes);
//    int num_select = 0, i = 0;      //�ڵ�������
//    srand((unsigned)time(NULL));
//    while (num_select < num_of_random_select) {        //�����ѡ��һ��ڵ�
//        while (true) {
//            i = rand() % nb_nodes;
//            if (node_value[i] == 1 && node_temp[i] == 0)break;
//        }
//        node_temp[i] = 1;
//        selectNodes.push_back(i);
//        hashone += hash_key_temp_one[i];      //ͬʱ�����µĽ��hashֵ
//        hashtwo += hash_key_temp_two[i];
//        hashthree += hash_key_temp_three[i];
//        num_select++;
//    }
//
//    //���ѡ��һ���ֽڵ�
//    //int num = 0;
//    //srand((unsigned)time(NULL));
//    //while (num < nb_sub_nodes - num_of_random_select) {        //�����ѡ��һ��ڵ�
//    //    while (true) {
//    //        i = rand() % nb_nodes;
//    //        if (local_best[i] == 0 && node_temp[i] == 0)break;
//    //    }
//    //    node_temp[i] = 1;
//    //    selectNodes.push_back(i);
//    //    hashone += hash_key_temp_one[i];      //ͬʱ�����µĽ��hashֵ
//    //    hashtwo += hash_key_temp_two[i];
//    //    hashthree += hash_key_temp_three[i];
//    //    num++;
//    //    test
//    //    new_nodes.push_back(i);
//    //    test end
//    //}
//
//    //̰����ѡ��һ���ֽڵ�
//    Distance temp_sum = (max_select_node + min_select_node) / 2;     //̰����ѡ��һ��ڵ�
//    priority_queue<pair<Distance, int>, vector<pair<Distance, int>>, greater<pair<Distance, int>> > q;
//    for (int m = 0; m < nb_nodes; ++m) {
//        if (node_value[m] == 0) {
//            q.push(make_pair(fabs(node_dis_sum[m] - temp_sum), m));
//        }
//    }
//    for (int m = 0; m < nb_sub_nodes - num_of_random_select; ++m) {
//        pair <Distance, int> first = q.top();
//        node_temp[first.second] = 1;
//        selectNodes.push_back(first.second);
//        hashone += hash_key_temp_one[first.second];      //ͬʱ�����µĽ��hashֵ
//        hashtwo += hash_key_temp_two[first.second];
//        hashthree += hash_key_temp_three[first.second];
//        q.pop();
//    }
//    hashone = hashone % size_of_tabu_list;
//    hashtwo = hashtwo % size_of_tabu_list;
//    hashthree = hashthree % size_of_tabu_list;
//    node_value = node_temp;
//    max = 0, min = DISTANCE_MAX;
//    for (int l = 0; l < nb_nodes; ++l) {       //����node_dis_num
//        node_dis_sum[l] = 0;
//        for (int m = 0; m < nb_sub_nodes; ++m) {
//            node_dis_sum[l] += ins.dis_nodes(l, selectNodes[m]);
//        }
//        if (node_value[l]) {
//            if (max < node_dis_sum[l])max = node_dis_sum[l];
//            if (min > node_dis_sum[l])min = node_dis_sum[l];
//        }
//    }
//    cur_obj = max - min;
//    no_select_nodes.clear();            //���¸����ṹselect_nodes��no_select_nodes
//    select_nodes.clear();
//    Distance cur_temp_sum = (max + min) / 2;
//    for (int i = 0; i < nb_nodes; ++i) {
//        Distance temp = fabs(node_dis_sum[i] - cur_temp_sum);
//        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
//        else no_select_nodes.push_back(make_pair(i, temp));
//    }
//    return true;
//}

//ѡ����ǰ�����ʷ���Ž�Ĺ�ͬ�ڵ㣬Ȼ�������ӽڵ�
//bool LocalSearch::stochastic_perturbation(int &hashone, int &hashtwo, int &hashthree, Distance &max, Distance &min) {
//    vector<int> node_temp(nb_nodes, 0);
//    hashone = 0;
//    hashtwo = 0;
//    hashthree = 0;
//    set<int> best_nodes;
//    set<int> cur_nodes;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (node_value[i])cur_nodes.insert(i);
//        if (local_best[i])best_nodes.insert(i);
//    }
//    vector<int> ivec(nb_sub_nodes); 
//    auto iter = set_intersection(best_nodes.begin(), best_nodes.end(), cur_nodes.begin(), cur_nodes.end(), ivec.begin());
//    ivec.resize(iter - ivec.begin());//����ȷ��ivec��С
//    for (int i = 0; i < ivec.size(); ++i) {
//        node_temp[ivec[i]] = 1;
//        hashone += hash_key_temp_one[ivec[i]];      //ͬʱ�����µĽ��hashֵ
//        hashtwo += hash_key_temp_two[ivec[i]];
//        hashthree += hash_key_temp_three[ivec[i]];
//    }
//    int remain = nb_sub_nodes - ivec.size();
//    int num = 0, i = 0;
//    ivec.resize(nb_nodes);
//    srand((unsigned)time(NULL));
//    while (num < remain) {        //�����ѡ��һ��ڵ�
//        while (true) {
//            i = rand() % nb_nodes;
//            if (local_best[i] == 0 && node_temp[i] == 0)break;
//        }
//        node_temp[i] = 1;
//        ivec.push_back(i);
//        hashone += hash_key_temp_one[i];      //ͬʱ�����µĽ��hashֵ
//        hashtwo += hash_key_temp_two[i];
//        hashthree += hash_key_temp_three[i];
//        num++;
//    }
//    
//    hashone = hashone % size_of_tabu_list;
//    hashtwo = hashtwo % size_of_tabu_list;
//    hashthree = hashthree % size_of_tabu_list;
//
//    node_value = node_temp;
//
//    max = 0, min = DISTANCE_MAX;
//    for (int l = 0; l < nb_nodes; ++l) {       //����node_dis_num
//        node_dis_sum[l] = 0;
//        for (int m = 0; m < nb_sub_nodes; ++m) {
//            node_dis_sum[l] += ins.dis_nodes(l, ivec[m]);
//        }
//        if (node_value[l]) {
//            if (max < node_dis_sum[l])max = node_dis_sum[l];
//            if (min > node_dis_sum[l])min = node_dis_sum[l];
//        }
//    }
//    cur_obj = max - min;
//
//    no_select_nodes.clear();            //���¸����ṹselect_nodes��no_select_nodes
//    select_nodes.clear();
//    Distance cur_temp_sum = (max + min) / 2;
//    for (int i = 0; i < nb_nodes; ++i) {
//        Distance temp = fabs(node_dis_sum[i] - cur_temp_sum);
//        if (node_value[i])select_nodes.push_back(make_pair(i, temp));
//        else no_select_nodes.push_back(make_pair(i, temp));
//    }
//    return true;
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
//int LocalSearch::hash_function_two() {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (node_value[i]) {
//            sum += (int)(floor(pow(i, hashFun_two_param)));
//        }
//    }
//    return sum % size_of_tabu_list;
//}
//int LocalSearch::hash_function_three() {
//    long long sum = 0;
//    for (int i = 0; i < nb_nodes; ++i) {
//        if (node_value[i]) {
//            sum += hash_key_temp_three[i];
//        }
//    }
//    return sum % size_of_tabu_list;
//}

}