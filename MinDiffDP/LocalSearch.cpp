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
    iter(0),                   //xxf:�޸����������ʱ��ms
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
    length_disturbance = (int)(nb_nodes / 2 * nb_sub_nodes);
}

Solution LocalSearch::solve() {
    Timer time(max_time);
    int fixed_value = 4000;   //MDG_a.b_1-20;GKD_c;
    //int fixed_value = 2000;    //MDG_a.b_21-40;DM1A
    //int fixed_value = 500;     //MDG_c
    int step_length = 0;           //��¼���ٲ�֮�ڲ��ܸĽ���ʷ���Ž�
    bool tabu_flag = false;          //�ж��Ƿ�����ⶼ�ڽ�����
    bool length_flag = false;
    clock_t start_time = clock();
    while (!time.isTimeOut()) {           //��ʱ��δ��ʱʱ�����оֲ����� 
        int _hashfun_one = best_hashfun_one;
        int _hashfun_two = best_hashfun_two;
        int _hashfun_three = best_hashfun_three;
        //����Ŷ�
        if (length_flag) {
            for (int i = 0; i < length_disturbance; ++i)
                    stochastic_perturbation(_hashfun_one, _hashfun_two, _hashfun_three, max_select_node, min_select_node);
            length_flag = false;
            step_length = 0;
        }
        else {
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
            if (step_length == fixed_value) {
                length_flag = true;
                step_length = 0;
                break;
            }
        }
        if (tabu_flag)break;
    }
    mylog << "\n�ܵ���������" << iter <<= logsw_local;
    clock_t end_time = clock();
    //�ؼ���Ϣ������ļ���
    ofstream outFile;                                
    outFile.open("../Deploy/Logs/log.csv", ios::app);
    outFile << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << ',' << nb_nodes << ',' << nb_sub_nodes << ',' << myrand.getSeed() << ',' << iter << ',' << local_best_obj << endl;
    outFile.close();
    //�ļ��������
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

//��һ������Ŷ�����:�ӵ�ǰ��ʷ���Ž�ĺ�ѡ���������
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
        int _new_hash_one = hashone + hash_key_temp_one[zero_toOne_node] - hash_key_temp_one[one_toZero_node];  //���㽻�����½�Ĺ�ϣ����ֵ
        int _new_hash_two = hashtwo + hash_key_temp_two[zero_toOne_node] - hash_key_temp_two[one_toZero_node];
        int _new_hash_three = hashthree + hash_key_temp_three[zero_toOne_node] - hash_key_temp_three[one_toZero_node];
        _new_hash_one = (_new_hash_one + size_of_tabu_list) % size_of_tabu_list;                //xxf�����bug2����ֹ���ָ�����>size_of_tabu_list����
        _new_hash_two = (_new_hash_two + size_of_tabu_list) % size_of_tabu_list;
        _new_hash_three = (_new_hash_three + size_of_tabu_list) % size_of_tabu_list;
        if (tabu_list_three[_new_hash_three]) {       //xxf��������������ϣ����Ϊ�˼�С��ͻ��ײ��һ����ϣ�������ײ�ͬ��ӳ�䵽ͬһ��valueֵ���Һ����׳�������ⶼ�������������Ϊ�Ǵ�ǿ���������Կ�ʼ
            if (tabu_list_two[_new_hash_two]) {
                if (tabu_list_one[_new_hash_one]);
                else break;
            }
            else break;
        }
        else break;
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
    node_value[one_toZero_node] = 0;                       //���µ�ǰ��
    node_value[zero_toOne_node] = 1;
    cur_obj = temp_max - temp_min;
    max_select_node = temp_max;
    min_select_node = temp_min;

    tabu_list_one[_new_hash_one] = 1;         //�������������б�
    tabu_list_two[_new_hash_two] = 1;
    tabu_list_three[_new_hash_three] = 1;
    hashone = _new_hash_one;
    hashtwo = _new_hash_two;
    hashthree = _new_hash_three;
        for (int i = 0; i < nb_nodes; ++i) {
            node_dis_sum[i] = node_dis_sum[i] - ins.dis_nodes(i, one_toZero_node) + ins.dis_nodes(i, zero_toOne_node);
        }
        no_select_nodes.clear();            //���¸����ṹselect_nodes��no_select_nodes
        select_nodes.clear();
        Distance cur_temp_sum = (max_select_node + min_select_node) / 2;
        for (int i = 0; i < nb_nodes; ++i) {
            Distance temp = fabs(node_dis_sum[i] - cur_temp_sum);
            if (node_value[i])select_nodes.push_back(make_pair(i, temp));
            else no_select_nodes.push_back(make_pair(i, temp));
        }
    //sort(no_select_nodes.begin(), no_select_nodes.end(), compareByAscend);    //δѡ�е���������
    return true;
}

}