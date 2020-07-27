//Contributors: Sibo Wang, Renchi Yang
#ifndef __ALGO_H__
#define __ALGO_H__

#include "graph.h"
#include "heap.h"
#include "config.h"
#include "rng.h"
#include <tuple>
#include <boost/random.hpp>
// #include "sfmt/SFMT.h"


struct PredResult{
    double topk_avg_relative_err;
    double topk_avg_abs_err;
    double topk_recall;
    double topk_precision;
    int real_topk_source_count;
    PredResult(double mae=0, double mre=0, double rec=0, double pre=0, int count=0):
        topk_avg_relative_err(mae),
        topk_avg_abs_err(mre),
        topk_recall(rec),
        topk_precision(pre),
        real_topk_source_count(count){}
};

unordered_map<int, PredResult> pred_results;

Fwdidx fwd_idx;
IFwdidx ifwd_idx;
Bwdidx bwd_idx;
iMap<double> ppr;
iMap<int> topk_filter;

// vector< boost::atomic<double> > vec_ppr;
iMap<int> rw_counter;
// RwIdx rw_idx;
atomic<unsigned long long> num_hit_idx;
atomic<unsigned long long> num_total_rw;
long num_iter_topk;
vector<int> rw_idx;
vector< pair<unsigned long long, unsigned long> > rw_idx_info;

map< int, vector< pair<int ,double> > > exact_topk_pprs;
vector< pair<int ,double> > topk_pprs;

iMap<double> upper_bounds;
iMap<double> lower_bounds;

vector<pair<int, double>> map_lower_bounds;
vector<double> index_use_ratio;


// for hubppr
vector<int> hub_fwd_idx;
//pointers to compressed fwd_idx, nodeid:{ start-pointer, start-pointer, start-pointer,...,end-pointer }
map<int, vector<unsigned long long> > hub_fwd_idx_cp_pointers;
vector<vector<unsigned long long>> hub_fwd_idx_ptrs;
vector<int> hub_fwd_idx_size;
iMap<int> hub_fwd_idx_size_k;
vector<int> hub_sample_number;
iMap<int> hub_counter;
PendingQueue pq;
map<int, vector<HubBwdidxWithResidual>> hub_bwd_idx;
vector<double> pg_values;

unsigned concurrency;

vector<int> ks;

vector<unordered_map<int, double>> residual_maps;
vector<map<int, double>> reserve_maps;

inline uint32_t xor128(void){
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;
    t = x ^ (x << 11);   
    x = y; y = z; z = w;   
    return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}

inline static unsigned long new_xshift_lrand(){
    static rng::rng128 rng;
    static uint64_t no_opt=0;
    no_opt |=rng();
    return no_opt;
}

inline static bool new_xshift_drand(){
        return ((double)new_xshift_lrand()/(double)UINT_MAX)< config.alpha;
}

inline static unsigned long xshift_lrand(){
    return (unsigned long)xor128();
}

inline static bool xshift_drand(){
    return ((double)xshift_lrand()/(double)UINT_MAX)<config.alpha;
}

inline static unsigned long lrand() {

    static boost::taus88 rngG(time(0));
    return rngG();

    //return rand();
    // return sfmt_genrand_uint32(&sfmtSeed);
}

inline static bool drand(){
    static boost::bernoulli_distribution <> bernoulli(config.alpha);
    static boost::lagged_fibonacci607 rngG(time(0));
    static boost::variate_generator<boost::lagged_fibonacci607&, boost::bernoulli_distribution<> > bernoulliRngG(rngG, bernoulli);

    return bernoulliRngG();
    //return rand()*1.0f/RAND_MAX;
    // return sfmt_genrand_real1(&sfmtSeed);
}

inline int random_walk(int start, const Graph& graph){
    int cur = start;
    unsigned long k;
    if(graph.g[start].size()==0){
        return start;
    }
    while (true) {
        if (drand()) {
            return cur;
        }
        if (graph.g[cur].size()){
            k = lrand()%graph.g[cur].size();
            cur = graph.g[cur][ k ];
        }
        else{
            cur = start;
        }
    }
}

inline int random_walk_no_zero_hop(int start, const Graph& graph){
    int cur = start;
    unsigned long k;
    if(graph.g[start].size()==0){
        return start;
    }

    k = lrand()%graph.g[cur].size();
    cur = graph.g[cur][ k ];

    while (true) {
        if (drand()) {
            return cur;
        }
        if (graph.g[cur].size()){
            k = lrand()%graph.g[cur].size();
            cur = graph.g[cur][ k ];
        }
        else{
            cur = start;
        }
    }
}


unsigned int SEED=1;
inline static unsigned long lrand_thd(int core_id) {
    //static thread_local std::mt19937 gen(core_id+1);
    //static std::uniform_int_distribution<> dis(0, INT_MAX);
    //return dis(gen);
    return rand_r(&SEED);
}

inline static double drand_thd(int core_id){
	return ((double)lrand_thd(core_id)/(double)INT_MAX);
}

inline int random_walk_thd(int start, const Graph& graph, int core_id){
    int cur = start;
    unsigned long k;
    if(graph.g[start].size()==0){
        return start;
    }
    while (true) {
        if (drand_thd(core_id) < config.alpha) {
            return cur;
        }
        if (graph.g[cur].size()){
            k = lrand_thd(core_id)%graph.g[cur].size();
            cur = graph.g[cur][ k ];
        }
        else{
            cur = start;
        }
    }
}

void count_hub_dest(){
    // {   Timer tm(101);
    int remaining;
    unsigned long long last_beg_ptr;
    unsigned long long end_ptr;
    int hub_node;
    int blocked_num;
    int bit_pos;
    for(int i=0; i<hub_counter.occur.m_num; i++){
        hub_node = hub_counter.occur[i];
        last_beg_ptr = hub_fwd_idx_ptrs[hub_node][hub_fwd_idx_ptrs[hub_node].size()-2];
        end_ptr = hub_fwd_idx_ptrs[hub_node][hub_fwd_idx_ptrs[hub_node].size()-1];

        remaining = hub_counter[hub_node];

        if(remaining > hub_fwd_idx_size_k[hub_node]){
            for(unsigned long long ptr=last_beg_ptr; ptr<end_ptr; ptr+=2){
                if(rw_counter.notexist(hub_fwd_idx[ptr])){
                    rw_counter.insert(hub_fwd_idx[ptr], hub_fwd_idx[ptr+1]);
                }
                else{
                    rw_counter[hub_fwd_idx[ptr]]+=hub_fwd_idx[ptr+1];
                }
                remaining-=hub_fwd_idx[ptr+1];
            }
        }


        for(int j=0; j< hub_fwd_idx_ptrs[hub_node].size()-2; j++){
            bit_pos = 1<<j;
            if(bit_pos & remaining){
                for(unsigned long long ptr=hub_fwd_idx_ptrs[hub_node][j]; ptr<hub_fwd_idx_ptrs[hub_node][j+1]; ptr+=2){
                    if(rw_counter.notexist(hub_fwd_idx[ptr])){
                        rw_counter.insert(hub_fwd_idx[ptr], hub_fwd_idx[ptr+1]);
                    }
                    else{
                        rw_counter[hub_fwd_idx[ptr]]+=hub_fwd_idx[ptr+1];
                    }
                }
            }
        }
    }
    // }
}

inline int random_walk_with_compressed_forward_oracle(int start, const Graph &graph) {
    int cur = start;
    unsigned long k;
    if(graph.g[start].size()==0){
        return start;
    }
    while (true) {
        if ( hub_fwd_idx_size[cur]!=0 && (hub_counter.notexist(cur) || hub_counter[cur] <  hub_fwd_idx_size[cur]) ){
            if(hub_counter.notexist(cur))
                hub_counter.insert(cur, 1);
            else
                hub_counter[cur]+=1;
            return -1;
        }

        if (drand() < config.alpha) {
            return cur;
        }

        if (graph.g[cur].size()){
            k = lrand()%graph.g[cur].size();
            cur = graph.g[cur][k];
        }
        else
            cur  = start;
    }
}

inline void generate_accumulated_fwd_randwalk(int s, const Graph& graph, unsigned long long num_rw){
    if(graph.g[s].size()==0){
        if(rw_counter.notexist(s)){
            rw_counter.insert(s, num_rw);
        }
        else{
            rw_counter[s] += num_rw;
        }
        return;
    }

    if(config.with_rw_idx){
        if(hub_fwd_idx_size[s]!=0){
            if(num_rw <= hub_fwd_idx_size[s]){
                hub_counter.insert(s,  num_rw);
                num_rw=0;
            }
            else{
                hub_counter.insert(s,  hub_fwd_idx_size[s]);
                num_rw-=hub_fwd_idx_size[s];
            }
            for(unsigned long long i=0; i<num_rw; i++){
                int t = random_walk_with_compressed_forward_oracle(s, graph);
                if(rw_counter.notexist(t)){
                    rw_counter.insert(t, 1);
                }
                else{
                    rw_counter[t]+=1;
                }
            }
        }
        else{
            rw_counter.insert(s, num_rw*config.alpha);
            // rw_counter[s] += num_rw*config.alpha;
            num_rw = num_rw*(1-config.alpha);
            for(unsigned long long i=0; i<num_rw; i++){
                int v = graph.g[s][lrand()%graph.g[s].size()];
                int t = random_walk_with_compressed_forward_oracle(v, graph);
                if(rw_counter.notexist(t)){
                    rw_counter.insert(t, 1);
                }
                else{
                    rw_counter[t]+=1;
                }
            }
        }

        count_hub_dest();
        hub_counter.clean();
    }
    else{
        for(unsigned long long i=0; i<num_rw; i++){
            int t = random_walk(s, graph);
            if(rw_counter.notexist(t)){
                rw_counter.insert(t, 1);
            }
            else{
                rw_counter[t]+=1;
            }
        }
    }
}

inline void split_line(){
    INFO("-----------------------------");
}

inline void display_setting(){
    INFO(config.delta);
    INFO(config.pfail);
    INFO(config.rmax);
    INFO(config.omega);
}

inline void display_fwdidx(){
    for(int i=0; i<fwd_idx.first.occur.m_num;i++){
        int nodeid = fwd_idx.first.occur[i];
        cout << "k:" << nodeid << " v:" << fwd_idx.first[nodeid] << endl;
    }

    cout << "=======================" << endl;

    for(int i=0; i<fwd_idx.second.occur.m_num;i++){
        int nodeid = fwd_idx.second.occur[i];
        cout << "k:" << nodeid << " v:" << fwd_idx.second[nodeid] << endl;
    }
}

inline void display_ppr(){
    for(int i=0; i< ppr.occur.m_num; i++){
        cout << ppr.occur[i] << "->" << ppr[ ppr.occur[i] ] << endl;
    }
}

static void display_time_usage(int used_counter, int query_size){
    if(config.algo == FORA){
        cout << "Total cost (s): " << Timer::used(used_counter) << endl;
        cout <<  Timer::used(RONDOM_WALK)*100.0/Timer::used(used_counter) << "%" << " for random walk cost" <<endl;
        cout <<  Timer::used(FWD_LU)*100.0/Timer::used(used_counter) << "%" << " for forward push cost" << endl;
        // if(config.action == TOPK)
            // cout <<  Timer::used(SORT_MAP)*100.0/Timer::used(used_counter) << "%" << " for sorting top k cost" << endl;
        split_line();
    }
    else if(config.algo == BIPPR){
        cout << "Total cost (s): " << Timer::used(used_counter) << endl;
        cout <<  Timer::used(RONDOM_WALK)*100.0/Timer::used(used_counter) << "%" << " for random walk cost" <<endl;
        cout <<  Timer::used(BWD_LU)*100.0/Timer::used(used_counter) << "%" << " for backward push cost" << endl;
    }
    else if(config.algo == MC){
        cout << "Total cost (s): " << Timer::used(used_counter) << endl;
        cout <<  Timer::used(RONDOM_WALK)*100.0/Timer::used(used_counter) << "%" << " for random walk cost" <<endl;
    }
    else if(config.algo == FWDPUSH){
        cout << "Total cost (s): " << Timer::used(used_counter) << endl;
        cout <<  Timer::used(FWD_LU)*100.0/Timer::used(used_counter) << "%" << " for forward push cost" << endl;
    }

    if(config.with_rw_idx)
        cout << "Average rand-walk idx hit ratio: " << num_hit_idx*100.0/num_total_rw << "%" << endl;
    
    if(config.action == TOPK){
        assert(result.real_topk_source_count>0);
        cout << "Average top-K Precision: " << result.topk_precision/result.real_topk_source_count << endl;
        cout << "Average top-K Recall: " << result.topk_recall/result.real_topk_source_count << endl;
    }
    
    cout << "Average query time (s):"<<Timer::used(used_counter)/query_size<<endl;
    cout << "Memory usage (MB):" << get_proc_memory()/1000.0 << endl << endl; 
}

static void set_result(const Graph& graph, int used_counter, int query_size){
    config.query_size = query_size;

    result.m = graph.m;
    result.n = graph.n;
    result.avg_query_time = Timer::used(used_counter)/query_size;

    result.total_mem_usage = get_proc_memory()/1000.0;
    result.total_time_usage = Timer::used(used_counter);

    result.num_randwalk = num_total_rw;
    
    if(config.with_rw_idx){
        result.num_rw_idx_use = num_hit_idx;
        result.hit_idx_ratio = num_hit_idx*1.0/num_total_rw;
    }

    result.randwalk_time = Timer::used(RONDOM_WALK);
    result.randwalk_time_ratio = Timer::used(RONDOM_WALK)*100/Timer::used(used_counter);

    if(config.algo == FORA){
        result.propagation_time = Timer::used(FWD_LU);
        result.propagation_time_ratio = Timer::used(FWD_LU)*100/Timer::used(used_counter);
        // result.source_dist_time = Timer::used(SOURCE_DIST);
        // result.source_dist_time_ratio = Timer::used(SOURCE_DIST)*100/Timer::used(used_counter);
    }
    else if(config.algo == BIPPR){
        result.propagation_time = Timer::used(BWD_LU);
        result.propagation_time_ratio = Timer::used(BWD_LU)*100/Timer::used(used_counter);
    }

    if(config.action == TOPK){
        result.topk_sort_time = Timer::used(SORT_MAP);
        // result.topk_precision = avg_topk_precision;
        // result.topk_sort_time_ratio = Timer::used(SORT_MAP)*100/Timer::used(used_counter);
    }
}

inline void bippr_setting(int n, long long m){
    config.rmax = config.epsilon*sqrt(m*1.0*config.delta/3.0/log(2.0/config.pfail));
    // config.omega = m/config.rmax;
    config.rmax *= config.rmax_scale;
    config.omega = config.rmax*3*log(2.0/config.pfail)/config.delta/config.epsilon/config.epsilon;
}

inline void hubppr_topk_setting(int n, long long m){
    config.rmax = config.epsilon*sqrt(m*1.0*config.delta/3.0/log(2.0/config.pfail));
    config.rmax *= config.rmax_scale;
    config.omega = config.rmax*3*log(2.0/config.pfail)/config.delta/config.epsilon/config.epsilon;
}

inline void fora_setting(int n, long long m){
    config.rmax = config.epsilon*sqrt(config.delta/3/m/log(2/config.pfail));
    if(config.opt)
        config.rmax *= config.rmax_scale/(1-config.alpha);
    else
        config.rmax *=config.rmax_scale;
    // config.rmax *= config.multithread_param;
    config.omega = (2+config.epsilon)*log(2/config.pfail)/config.delta/config.epsilon/config.epsilon;
}


inline void fora_topk_setting(int n, long long m){
    config.rmax = config.epsilon*sqrt(config.delta/3/m/log(2/config.pfail));
    if(config.opt)
        config.rmax *= sqrt( 1.0*m*config.rmax)*config.rmax_scale*3;
    else
        config.rmax *= sqrt( 1.0*m*config.rmax)*config.rmax_scale*3;
    // config.rmax *= config.multithread_param;
    config.omega = (2+config.epsilon)*log(2/config.pfail)/config.delta/config.epsilon/config.epsilon;
}


inline void montecarlo_setting(){
    double fwd_rw_count = 3*log(2/config.pfail)/config.epsilon/config.epsilon/config.delta;
    config.omega = fwd_rw_count;
    // INFO(config.delta);
    // INFO(config.pfail);
    // INFO(config.omega);
}

inline void fwdpush_setting(int n, long long m){
    // below is just a estimate value, has no accuracy guarantee
    // since for undirected graph, error |ppr(s, t)-approx_ppr(s, t)| = sum( r(s, v)*ppr(v, t)) <= d(t)*rmax
    // |ppr(s, t)-approx_ppr(s, t)| <= epsilon*ppr(s, t)
    // d(t)*rmax <= epsilon*ppr(s, t)
    // rmax <= epsilon*ppr(s, t)/d(t)
    // d(t): use average degree d=m/n
    // ppr(s, t): use minimum ppr value, delta, i.e., 1/n
    // thus, rmax <= epsilon*delta*n/m = epsilon/m
    // use config.rmax_scale to tune rmax manually
    config.rmax = config.rmax_scale*config.delta*config.epsilon*n/m;
}

inline void generate_ss_query(int n){
    string filename = config.graph_location + "ssquery.txt";
    if(exists_test(filename)){
        INFO("ss query set exists");
        return;
    }
    ofstream queryfile(filename);
    for(int i=0; i<config.query_size; i++){
        int v = rand()%n;
        queryfile<<v<<endl;
    }
}

void load_ss_query(vector<int>& queries){
    string filename = config.graph_location+"ssquery.txt";
     if(!file_exists_test(filename)){
        cerr<<"query file does not exist, please generate ss query files first"<<endl;
        exit(0);
    }
    ifstream queryfile(filename);
    int v;
    while(queryfile>>v){
        queries.push_back(v);
    }
}

void compute_precision(int v){
    double precision=0.0;
    double recall=0.0;
    //INFO(topk_pprs.size());
    if( exact_topk_pprs.size()>=1 && exact_topk_pprs.find(v)!=exact_topk_pprs.end() ){

        unordered_map<int, double> topk_map;
        for(auto &p: topk_pprs){
            if(p.second>0){
                topk_map.insert(p);
            }
        }

        unordered_map<int, double> exact_map;
        int size_e = min( config.k, (unsigned int)exact_topk_pprs[v].size() );

        for(int i=0; i<size_e; i++){
            pair<int ,double>& p = exact_topk_pprs[v][i];
            if(p.second>0){
                exact_map.insert(p);
                if(topk_map.find(p.first)!=topk_map.end())
                    recall++;
            }
        }

        for(auto &p: topk_map){
            if(exact_map.find(p.first)!=exact_map.end()){
                precision++;
            }
        }

        // for(int i=0; i<config.k; i++){
        //     cout << "NO." << i << " pred:" << topk_pprs[i].first << ", " << topk_pprs[i].second << "\t exact:" << exact_topk_pprs[v][i].first << ", " << exact_topk_pprs[v][i].second << endl;
        // }

        assert(exact_map.size() > 0);
        assert(topk_map.size() > 0);
        //if(exact_map.size()<=1)
        //    return;

        recall = recall*1.0/exact_map.size();
        precision = precision*1.0/exact_map.size();
        INFO(exact_map.size(), recall, precision);
        result.topk_recall += recall;
        result.topk_precision += precision;

        result.real_topk_source_count++;
    }
}

inline bool cmp(double x, double y){
    return x>y;
}
// obtain the top-k ppr values from ppr map
double kth_ppr(){
    Timer tm(SORT_MAP);

    static vector<double> temp_ppr;
    temp_ppr.clear();
    temp_ppr.resize(ppr.occur.m_num);
    for(int i; i<ppr.occur.m_num; i++){
        temp_ppr[i] = ppr[ ppr.occur[i] ];
    }

    nth_element(temp_ppr.begin(), temp_ppr.begin()+config.k-1, temp_ppr.end(), cmp);
    return temp_ppr[config.k-1];
}

double topk_ppr(){
    topk_pprs.clear();
    topk_pprs.resize(config.k);

    static unordered_map< int, double > temp_ppr;
    temp_ppr.clear();
    // temp_ppr.resize(ppr.occur.m_num);
    int nodeid;
    for(long i=0; i<ppr.occur.m_num; i++){
        nodeid = ppr.occur[i];
        // INFO(nodeid);
        temp_ppr[nodeid] = ppr[ nodeid ];
    }

    partial_sort_copy(temp_ppr.begin(), temp_ppr.end(), topk_pprs.begin(), topk_pprs.end(), 
            [](pair<int, double> const& l, pair<int, double> const& r){return l.second > r.second;});
    
    return topk_pprs[config.k-1].second;
}

double topk_of(vector< pair<int, double> >& top_k){
    top_k.clear();
    top_k.resize(config.k);
    static vector< pair<int, double> > temp_ppr;
    temp_ppr.clear();
    temp_ppr.resize(ppr.occur.m_num);
    for(int i=0; i<ppr.occur.m_num; i++){
        temp_ppr[i] = MP(ppr.occur[i], ppr[ ppr.occur[i] ]);
    }

    partial_sort_copy(temp_ppr.begin(), temp_ppr.end(), top_k.begin(), top_k.end(), 
        [](pair<int, double> const& l, pair<int, double> const& r){return l.second > r.second;});
    
    return top_k[config.k-1].second;
}

void compute_precision_for_dif_k(int v){
    if( exact_topk_pprs.size()>=1 && exact_topk_pprs.find(v)!=exact_topk_pprs.end() ){
        for(auto k: ks){

            int j=0;
            unordered_map<int, double> topk_map;
            for(auto &p: topk_pprs){
                if(p.second>0){
                    topk_map.insert(p);
                }
                j++;
                if(j==k){ // only pick topk
                    break;
                }
            }

            double recall=0.0;
            unordered_map<int, double> exact_map;
            int size_e = min( k, (int)exact_topk_pprs[v].size() );
            for(int i=0; i<size_e; i++){
                pair<int ,double>& p = exact_topk_pprs[v][i];
                if(p.second>0){
                    exact_map.insert(p);
                    if(topk_map.find(p.first)!=topk_map.end())
                        recall++;
                }
            }

            double precision=0.0;
            for(auto &p: topk_map){
                if(exact_map.find(p.first)!=exact_map.end()){
                    precision++;
                }
            }

            //if(exact_map.size()<=1)
            //   continue;

            precision = precision*1.0/exact_map.size();
            recall = recall*1.0/exact_map.size();

            pred_results[k].topk_precision += precision;
            pred_results[k].topk_recall += recall;
            pred_results[k].real_topk_source_count++;
        }
    }
}

inline void display_precision_for_dif_k(){
    split_line();
    cout << config.algo << endl;
    for(auto k: ks){
        cout << k << "\t";
    }
    cout << endl << "Precision:" << endl;
    //assert(pred_results[k].real_topk_source_count>0);
    for(auto k: ks){
        cout << pred_results[k].topk_precision/pred_results[k].real_topk_source_count << "\t";
    }
    cout << endl << "Recall:" << endl;
    for(auto k: ks){
        cout << pred_results[k].topk_recall/pred_results[k].real_topk_source_count << "\t";
    }
    cout << endl;
}

inline void init_multi_setting(int n){
    INFO("multithreading mode...");
    concurrency = std::thread::hardware_concurrency()+1;
    INFO(concurrency);
    assert(concurrency >= 2);
    config.rmax_scale = sqrt( concurrency*1.0 );
    INFO(config.rmax_scale);
}

static void reverse_local_update_linear(int t, const Graph &graph, double init_residual = 1) {
    bwd_idx.first.clean();
    bwd_idx.second.clean();

    static unordered_map<int, bool> idx;
    idx.clear();

    vector<int> q;
    q.reserve(graph.n);
    q.push_back(-1);
    unsigned long left = 1;

    double myeps = config.rmax;

    q.push_back(t);
    bwd_idx.second.insert(t, init_residual);

    idx[t] = true;
    while (left < q.size()) {
        int v = q[left];
        idx[v] = false;
        left++;
        if (bwd_idx.second[v] < myeps)
            break;

        if(bwd_idx.first.notexist(v))
            bwd_idx.first.insert(v, bwd_idx.second[v]*config.alpha);
        else
            bwd_idx.first[v] += bwd_idx.second[v]*config.alpha;

        double residual = (1 - config.alpha) * bwd_idx.second[v];
        bwd_idx.second[v] = 0;
        if(graph.gr[v].size()>0){
            for (int next : graph.gr[v]) {
                int cnt = graph.g[next].size();
                if(bwd_idx.second.notexist(next))
                    bwd_idx.second.insert(next, residual/cnt);
                else
                    bwd_idx.second[next] += residual/cnt;

                if (bwd_idx.second[next] > myeps && idx[next] != true) {
                    // put next into q if next is not in q
                    idx[next] = true;//(int) q.size();
                    q.push_back(next);
                }
            }
        }
    }
}

static void reverse_local_update_heap(int t, const Graph &graph, double init_residual = 1) {
    static BinaryHeap<double, greater<double> > heap(graph.n, greater<double>());

    bwd_idx.first.clean();
    bwd_idx.second.clean();

    double myeps = config.rmax;

    heap.clear();
    heap.insert(t, init_residual);

    while (heap.size()) {
        auto top = heap.extract_top();
        double residual = top.first;
        int v = top.second;
        if (residual < myeps)
            break;

        heap.delete_top();
        if(bwd_idx.first.notexist(v)){
            bwd_idx.first.insert(v, residual * config.alpha);
        }
        else{
            bwd_idx.first[v]+=residual * config.alpha;
        }
        double resi = (1 - config.alpha) * residual;
        for (int next : graph.gr[v]) {
            int cnt = graph.g[next].size();
            double delta = resi / cnt;
            if (heap.has_idx(next))
                heap.modify(next, heap.get_value(next) + delta);
            else
                heap.insert(next, delta);
        }
    }

    for(auto item: heap.get_elements()){
        bwd_idx.second.insert(item.second, item.first);
    }
}

void fwd_with_hub_oracle(const Graph& graph, int start){
    num_total_rw += config.omega;

    rw_counter.clean();
    hub_counter.clean();

    unsigned long long num_tries = int64(config.omega/graph.n*hub_sample_number[start]);
    // INFO(num_tries, config.omega);
    
    if(graph.g[start].size()==0) {
        if(rw_counter.notexist(start))
            rw_counter.insert(start, num_tries);
        else
            rw_counter[start] +=num_tries;
        return;
    }

    if(hub_fwd_idx_size[start]!=0){
        if(num_tries <= hub_fwd_idx_size[start]){
            hub_counter.insert(start, num_tries);
            num_tries=0;
        }
        else{
            hub_counter.insert(start, hub_fwd_idx_size[start]);
            num_tries-=hub_fwd_idx_size[start];
        }
    }else{
        if(rw_counter.notexist(start))
            rw_counter.insert(start, num_tries*config.alpha);
        else
            rw_counter[start] += num_tries*config.alpha;
        num_tries = num_tries*(1-config.alpha);
    }

    for (int64 i = 0; i < num_tries; i++) {
        if(hub_fwd_idx_size[start] !=0){
            int l = random_walk_with_compressed_forward_oracle(start, graph);
            if(l >=0){
                if(rw_counter.notexist(l)){
                    rw_counter.insert(l, 1);
                }
                else{
                    rw_counter[l]+=1;
                }
            }
        }else{
            int random_out_neighbor = drand()*(graph.g[start].size()-1);
            random_out_neighbor = graph.g[start][random_out_neighbor];
            int l = random_walk_with_compressed_forward_oracle(random_out_neighbor, graph);
            if(l >=0){
                if(rw_counter.notexist(l)){
                    rw_counter.insert(l, 1);
                }
                else{
                    rw_counter[l]+=1;
                }
            }
        }
    }
}

void bwd_with_hub_oracle(const Graph& graph, int t){
    bwd_idx.first.clean();
    bwd_idx.second.clean();

    static unordered_map<int, bool> idx;
    idx.clear();

    // static vector<bool> idx(graph.n);
    // std::fill(idx.begin(), idx.end(), false);

    vector<int> q;
    q.reserve(graph.n);
    q.push_back(-1);
    unsigned long left = 1;

    double myeps = config.rmax;
    // residual.clear();
    // exist.clear();

    q.push_back(t);
    // residual[t] = init_residual;
    bwd_idx.second.insert(t, 1);

    idx[t] = true;
    while (left < q.size()) {
        int v = q[left];
        idx[v] = false;
        left++;
        if (bwd_idx.second[v] < myeps)
            break;
        
        if (hub_bwd_idx.find(v) != hub_bwd_idx.end()) {
            vector<HubBwdidxWithResidual> &idxv = hub_bwd_idx[v];
            for (int i = idxv.size()-1; i >= 0; i--) {
                HubBwdidxWithResidual &x = idxv[i];
                if (x.first >= bwd_idx.second[v]) {
                    HubBwdidx &useidx = x.second;
                    for (auto &residualkv:useidx.first) {
                        int next = residualkv.first;
                        double delta = residualkv.second * bwd_idx.second[v] / x.first;
                        if(bwd_idx.first.notexist(next)){
                            bwd_idx.first.insert(next, delta);
                        }
                        else{
                            bwd_idx.first[next]+=delta;
                        }
                    }
                    if (useidx.second.size() <= 1){
                        bwd_idx.second[v] = 0;
                        break;
                    }
                    for (auto &residualkv:useidx.second) {
                        int next = residualkv.first;
                        double delta = residualkv.second * bwd_idx.second[v] / x.first;

                        if(bwd_idx.second.notexist(next))
                            bwd_idx.second.insert(next, delta);
                        else
                            bwd_idx.second[next] +=delta;

                        if (bwd_idx.second[next] >= myeps && idx[next] != true) {
                            // put next into q if next is not in q
                            idx[next] = true;//(int) q.size();
                            q.push_back(next);
                        }
                    }
                    bwd_idx.second[v] = 0;
                    break;
                }
            }
        }
        else{
            if(bwd_idx.first.notexist(v))
                bwd_idx.first.insert(v, bwd_idx.second[v]*config.alpha);
            else
                bwd_idx.first[v] += bwd_idx.second[v]*config.alpha;

            double residual = (1 - config.alpha) * bwd_idx.second[v];
            
            for (int next : graph.gr[v]) {
                int cnt = graph.g[next].size();
                // residual[next] += ((1 - config.alpha) * residual[v]) / cnt;
                if(bwd_idx.second.notexist(next))
                    bwd_idx.second.insert(next, residual/cnt);
                else
                    bwd_idx.second[next] += residual/cnt;

                if (bwd_idx.second[next] > myeps && idx[next] != true) {
                    // put next into q if next is not in q
                    idx[next] = true;//(int) q.size();
                    q.push_back(next);
                }
            }
            bwd_idx.second[v] = 0;
        }

    }
}

void forward_local_update_linear(int s, const Graph &graph, double& rsum, double rmax, double init_residual = 1.0){
    fwd_idx.first.clean();
    fwd_idx.second.clean();

    static vector<bool> idx(graph.n);
    std::fill(idx.begin(), idx.end(), false);

    if(graph.g[s].size()==0){
        fwd_idx.first.insert( s, 1);
        rsum =0;
        return; 
    }

    double myeps = rmax;//config.rmax;

    vector<int> q;  //nodes that can still propagate forward
    q.reserve(graph.n);
    q.push_back(-1);
    unsigned long left = 1;
    q.push_back(s);

    // residual[s] = init_residual;
    fwd_idx.second.insert(s, init_residual);
    
    idx[s] = true;
    
    while (left < (int) q.size()) {
        int v = q[left];
        idx[v] = false;
        left++;
        double v_residue = fwd_idx.second[v];
        fwd_idx.second[v] = 0;
        if(!fwd_idx.first.exist(v))
            fwd_idx.first.insert( v, v_residue * config.alpha);
        else
            fwd_idx.first[v] += v_residue * config.alpha;

        int out_neighbor = graph.g[v].size();
        rsum -=v_residue*config.alpha;
        if(out_neighbor == 0){
            fwd_idx.second[s] += v_residue * (1-config.alpha);
            if(graph.g[s].size()>0 && fwd_idx.second[s]/graph.g[s].size() >= myeps && idx[s] != true){
                idx[s] = true;
                q.push_back(s);   
            }
            continue;
        }

        double avg_push_residual = ((1.0 - config.alpha) * v_residue) / out_neighbor;
        for (int next : graph.g[v]) {
            // total_push++;
            if( !fwd_idx.second.exist(next) )
                fwd_idx.second.insert( next,  avg_push_residual);
            else
                fwd_idx.second[next] += avg_push_residual;

            //if a node's' current residual is small, but next time it got a laerge residual, it can still be added into forward list
            //so this is correct
            if ( fwd_idx.second[next]/graph.g[next].size() >= myeps && idx[next] != true) {  
                idx[next] = true;//(int) q.size();
                q.push_back(next);    
            }
        }
    }
}

void forward_local_update_linear_topk(int s, const Graph &graph, double& rsum, double rmax, double lowest_rmax, vector<int>& forward_from){
    double myeps = rmax;

    static vector<bool> in_forward(graph.n);
    static vector<bool> in_next_forward(graph.n);

    std::fill(in_forward.begin(), in_forward.end(), false);
    std::fill(in_next_forward.begin(), in_next_forward.end(), false);

    vector<int> next_forward_from;
    next_forward_from.reserve(graph.n);
    for(auto &v: forward_from)
        in_forward[v] = true;

    unsigned long i=0;
    while( i < forward_from.size() ){
        int v = forward_from[i];
        i++;
        in_forward[v] = false;
        if( fwd_idx.second[v]/graph.g[v].size() >= myeps ){
            int out_neighbor = graph.g[v].size();
            double v_residue = fwd_idx.second[v];
            fwd_idx.second[v] = 0;
            if(!fwd_idx.first.exist(v)){
                fwd_idx.first.insert( v, v_residue * config.alpha );
            }
            else{
                fwd_idx.first[v] += v_residue * config.alpha;
            }        

            rsum -= v_residue*config.alpha;
            if(out_neighbor == 0){
                fwd_idx.second[s] += v_residue * (1-config.alpha);
                if(graph.g[s].size()>0 && in_forward[s]!=true && fwd_idx.second[s]/graph.g[s].size() >= myeps){
                    forward_from.push_back(s);
                    in_forward[s] = true;
                }
                else{
                    if(graph.g[s].size()>=0 && in_next_forward[s]!=true && fwd_idx.second[s]/graph.g[s].size() >= lowest_rmax){
                        next_forward_from.push_back(s);
                        in_next_forward[s] = true;
                    }
                }
                continue;
            }            
            double avg_push_residual = ((1 - config.alpha) * v_residue) / out_neighbor;
            for( int next: graph.g[v] ){
                if(!fwd_idx.second.exist(next))
                    fwd_idx.second.insert(next, avg_push_residual);
                else
                    fwd_idx.second[next] += avg_push_residual;

                if( in_forward[next]!=true && fwd_idx.second[next]/graph.g[next].size() >= myeps ){
                    forward_from.push_back(next);
                    in_forward[next] = true;
                }
                else{
                    if(in_next_forward[next]!=true && fwd_idx.second[next]/graph.g[next].size() >= lowest_rmax){
                        next_forward_from.push_back(next);
                        in_next_forward[next] = true;
                    }
                }
            }
        }
        else{
            if(in_next_forward[v]!=true &&  fwd_idx.second[v]/graph.g[v].size() >= lowest_rmax){
                next_forward_from.push_back(v);
                in_next_forward[v] = true;
            }
        }
    }
    
    forward_from = next_forward_from;
}

extern double threshold;
bool if_stop(){
    // Timer tm(SORT_MAP);

    if(kth_ppr()>=2.0*config.delta)
        return true;

    if(config.delta>=threshold) return false;
    
    const static double error = 1.0+config.epsilon;
    const static double error_2 = 1.0+config.epsilon;

    topk_pprs.clear();
    topk_pprs.resize(config.k);
    topk_filter.clean();

    static vector< pair<int, double> > temp_bounds;
    temp_bounds.clear();
    temp_bounds.resize(lower_bounds.occur.m_num);
    int nodeid;
    for(int i=0; i<lower_bounds.occur.m_num; i++){
        nodeid = lower_bounds.occur[i];
        temp_bounds[i] = MP( nodeid, lower_bounds[nodeid] );
    }

    //sort topk nodes by lower bound
    partial_sort_copy(temp_bounds.begin(), temp_bounds.end(), topk_pprs.begin(), topk_pprs.end(), 
        [](pair<int, double> const& l, pair<int, double> const& r){return l.second > r.second;});
    
    //for topk nodes, upper-bound/low-bound <= 1+epsilon
    double ratio=0.0;
    double largest_ratio=0.0;
    for(auto &node: topk_pprs){
        topk_filter.insert(node.first, 1);
        ratio = upper_bounds[node.first]/lower_bounds[node.first];
        if(ratio>largest_ratio)
            largest_ratio = ratio;
        if(ratio>error_2){
            return false;
        }
    }

    // INFO("ratio checking passed--------------------------------------------------------------");
    //for remaining NO. k+1 to NO. n nodes, low-bound of k > the max upper-bound of remaining nodes 
    /*int actual_exist_ppr_num = lower_bounds.occur.m_num;
    if(actual_exist_ppr_num == 0) return true;
    int actual_k = min(actual_exist_ppr_num-1, (int) config.k-1);
    double low_bound_k = topk_pprs[actual_k].second;*/
    double low_bound_k = topk_pprs[config.k-1].second;
    if(low_bound_k<= config.delta){
        return false;
    }
    for(int i=0; i<upper_bounds.occur.m_num; i++){
        nodeid = upper_bounds.occur[i];
        if(topk_filter.exist(nodeid) || ppr[nodeid]<=0)
            continue;

        double upper_temp = upper_bounds[nodeid];
        double lower_temp = lower_bounds[nodeid];
        if(upper_temp> low_bound_k*error){
             if(upper_temp > (1+config.epsilon)/(1-config.epsilon)*lower_temp)
                 continue;
             else {
                 return false;
             }
        }else{
            continue;
        }
    }

    return true;
}


inline double calculate_lambda(double rsum, double pfail, double upper_bound, long total_rw_num){
    return 1.0/3*log(2/pfail)*rsum/total_rw_num + 
    sqrt(4.0/9.0*log(2.0/pfail)*log(2.0/pfail)*rsum*rsum +
        8*total_rw_num*log(2.0/pfail)*rsum*upper_bound)
    /2.0/total_rw_num;
}

double zero_ppr_upper_bound = 1.0;
double threshold = 0.0;
void set_ppr_bounds(const Graph& graph, double rsum, long total_rw_num){
    Timer tm(100);

    const static double min_ppr = 1.0/graph.n;
    const static double sqrt_min_ppr = sqrt(1.0/graph.n);


    double epsilon_v_div = sqrt(2.67*rsum*log(2.0/config.pfail)/total_rw_num);
    double default_epsilon_v = epsilon_v_div/sqrt_min_ppr;

    int nodeid;
    double ub_eps_a;
    double lb_eps_a;
    double ub_eps_v;
    double lb_eps_v;
    double up_bound;
    double low_bound;
    // INFO(total_rw_num);
    // INFO(zero_ppr_upper_bound);
    //INFO(rsum, 1.0/config.pfail, log(2/config.pfail), zero_ppr_upper_bound, total_rw_num);
    zero_ppr_upper_bound = calculate_lambda( rsum,  config.pfail,  zero_ppr_upper_bound,  total_rw_num );
    for(long i=0; i<ppr.occur.m_num; i++){
        nodeid = ppr.occur[i];
        if(ppr[nodeid]<=0)
            continue;
        double reserve=0.0;
        if(fwd_idx.first.exist(nodeid))
            reserve = fwd_idx.first[nodeid];
        double epsilon_a = 1.0;
        if( upper_bounds.exist(nodeid)  ){
            assert(upper_bounds[nodeid]>0.0);
            if(upper_bounds[nodeid] > reserve)
                //epsilon_a = calculate_lambda( rsum, config.pfail, upper_bounds[nodeid] - reserve, total_rw_num);
                epsilon_a = calculate_lambda( rsum, config.pfail, upper_bounds[nodeid] - reserve, total_rw_num);
            else
                epsilon_a = calculate_lambda( rsum, config.pfail, 1 - reserve, total_rw_num);
        }else{
            /*if(zero_ppr_upper_bound > reserve)
                epsilon_a = calculate_lambda( rsum, config.pfail, zero_ppr_upper_bound-reserve, total_rw_num);
            else
                epsilon_a = calculate_lambda( rsum, config.pfail, 1.0-reserve, total_rw_num);*/
            epsilon_a = calculate_lambda( rsum, config.pfail, 1.0-reserve, total_rw_num);
        }

        ub_eps_a = ppr[nodeid]+epsilon_a;
        lb_eps_a = ppr[nodeid]-epsilon_a;
        if(!(lb_eps_a>0))
            lb_eps_a = 0;

        double epsilon_v=default_epsilon_v;
        if(fwd_idx.first.exist(nodeid) && fwd_idx.first[nodeid]>min_ppr){
            if(lower_bounds.exist(nodeid))
                reserve = max(reserve, lower_bounds[nodeid]);
            epsilon_v = epsilon_v_div/sqrt(reserve);
        }else{
            if(lower_bounds[nodeid]>0)
                epsilon_v = epsilon_v_div/sqrt(lower_bounds[nodeid]);
        }


        ub_eps_v = 1.0;
        lb_eps_v = 0.0;
        if(1.0 - epsilon_v > 0){
            ub_eps_v = ppr[nodeid] / (1.0 - epsilon_v);
            lb_eps_v = ppr[nodeid] / (1.0 + epsilon_v);
        }

        up_bound = min( min(ub_eps_a, ub_eps_v), 1.0 );
        low_bound = max( max(lb_eps_a, lb_eps_v), reserve );
        if(up_bound>0){
            if(!upper_bounds.exist(nodeid))
                upper_bounds.insert(nodeid, up_bound);
            else
                upper_bounds[nodeid] = up_bound;
        }
        
        if(low_bound>=0){
            if(!lower_bounds.exist(nodeid))
                lower_bounds.insert(nodeid, low_bound);
            else
                lower_bounds[nodeid] = low_bound;
        }
    }
}

void set_martingale_bound(double lambda, unsigned long long total_num_rw, int t, double reserve, double rmax, 
                          double pfail_star, double min_ppr, double m_omega){
    double upper;
    double lower;
    double ey = total_num_rw*max( min_ppr-reserve, map_lower_bounds[t].second-reserve );
    double epsilon;
    if(m_omega>0 && ey>0){
        epsilon = (sqrt(pow(rmax*pfail_star/3, 2) - 2*ey*pfail_star) - rmax*pfail_star/3)/ey;
        upper = m_omega+lambda;
        if( 1-epsilon>0 && upper > m_omega/(1-epsilon) )
            upper = m_omega/(1-epsilon);
        upper /= total_num_rw;
        upper += reserve;
        lower = reserve + max( m_omega/(1+epsilon), m_omega-lambda )/total_num_rw;
    }
    else{
        upper = reserve + (m_omega + lambda)/total_num_rw;
        lower = reserve + (m_omega - lambda)/total_num_rw;
    }

    // INFO(m_omega, ey, epsilon, lower, reserve);

    if(upper > 0 && upper < upper_bounds[t])
        upper_bounds[t] = upper;
    if(lower < 1 && lower > map_lower_bounds[t].second)
        map_lower_bounds[t].second = lower;
}

int reverse_local_update_topk(int s, int t, map<int, double>& map_reserve, double rmax, unordered_map<int, double>& map_residual, const Graph &graph) {
    static vector<bool> bwd_flag(graph.n);

    vector<int> bwd_list;
    bwd_list.reserve(graph.n);
    unsigned long ptr = 0;

    for(auto &item: map_residual){
        if( item.second >= rmax ){
            bwd_list.push_back(item.first);
            bwd_flag[item.first] = true;
        }
    }

    int push_count=0;
    int v;
    while (ptr < bwd_list.size()) {
        v = bwd_list[ptr];
        bwd_flag[v] = false;
        ptr++;

        map_reserve[v] += map_residual[v] * config.alpha;

        double resi = (1 - config.alpha) * map_residual[v];
        map_residual.erase(v);
        push_count+=graph.gr[v].size();
        for (int next: graph.gr[v]) {
            long cnt = graph.g[next].size();
            map_residual[next] += resi / cnt;
            if( map_residual[next]>=rmax && bwd_flag[next]!=true ){
                bwd_flag[next] = true;
                bwd_list.push_back(next);
            }
        }
    }
    return push_count;
}

#endif
