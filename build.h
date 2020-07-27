//Contributors: Sibo Wang, Renchi Yang
#ifndef BUILD_H
#define BUILD_H

#include "algo.h"
#include "graph.h"
#include "heap.h"
#include "config.h"

#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <sys/sysinfo.h>

inline size_t get_ram_size(){
    struct sysinfo si;
    sysinfo (&si);
    return si.totalram;
}

inline string get_hub_fwd_idx_file_name(){
    string prefix = config.prefix + FILESEP + config.graph_alias+FILESEP;
    prefix += config.graph_alias + ".eps-" + to_str(config.epsilon);
    // prefix += ".space-1";
    prefix += ".space-" + to_str(config.hub_space_consum);

    string suffix;

    suffix += ".compress.fwd.idx";
    string file_name = prefix + suffix;
    return file_name;
}

inline string get_hub_fwd_idx_info_file_name(){
    string idx_file = get_hub_fwd_idx_file_name();
    return replace(idx_file, "fwd.idx", "fwd.info");
}

inline string get_hub_bwd_idx_file_name(){
    string idx_file = get_hub_fwd_idx_file_name();
    return replace(idx_file, "fwd.idx", "bwd.idx");
}

inline void deserialize_hub_fwd_idx(){
    string file_name = get_hub_fwd_idx_file_name();
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    ia >> hub_fwd_idx;
    
    string rwn_file_name = get_hub_fwd_idx_file_name()+".rwn";
    assert_file_exist("rwn file", rwn_file_name);
    std::ifstream ofs_rwn(rwn_file_name);
    boost::archive::binary_iarchive ia_rwn(ofs_rwn);

    ia_rwn >> hub_sample_number;


    string info_file = get_hub_fwd_idx_info_file_name();
    assert_file_exist("info file", info_file);
    std::ifstream info_ofs(info_file);
    boost::archive::binary_iarchive info_ia(info_ofs);
    info_ia >> hub_fwd_idx_cp_pointers;
}

inline void deserialize_hub_bwd_idx(){
    string file_name = get_hub_bwd_idx_file_name();
    // assert_file_exist("index file", file_name);
    if (!exists_test(file_name)) {
        cerr << "index file " << file_name << " not find " << endl;
        return;
    }
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    ia >> hub_bwd_idx;
}

void load_hubppr_oracle(const Graph& graph){
    deserialize_hub_fwd_idx();

    // fwd_idx_size.initialize(graph.n);
    hub_fwd_idx_ptrs.resize(graph.n);
    hub_fwd_idx_size.resize(graph.n);
    std::fill(hub_fwd_idx_size.begin(), hub_fwd_idx_size.end(), 0);
    hub_fwd_idx_size_k.initialize(graph.n);


    for(auto &ptrs: hub_fwd_idx_cp_pointers){
        int node = ptrs.first;
        int size=0;

        unsigned long long ptr = ptrs.second[0];
        unsigned long long end_ptr = ptrs.second[ptrs.second.size()-1];
        for(; ptr<end_ptr; ptr+=2){
            size += hub_fwd_idx[ptr+1];
        }

        hub_fwd_idx_ptrs[node] = ptrs.second;

        // fwd_idx_size.insert(node, size);
        hub_fwd_idx_size[node] = size;

        int u = 1 + floor(log( hub_fwd_idx_size[node]*1.0 )/log(2)); //we can pre-compute to avoid reduplicate computation
        int k = pow(2, u-1)-1;
        hub_fwd_idx_size_k.insert(node, k);
    }

    hub_fwd_idx_cp_pointers.clear();

    INFO(hub_fwd_idx_size.size());

    deserialize_hub_bwd_idx();
    INFO(hub_bwd_idx.size());
}

inline string get_exact_topk_ppr_file(){
    if(!boost::algorithm::ends_with(config.exact_pprs_folder, FILESEP))
        config.exact_pprs_folder += FILESEP;
    return config.exact_pprs_folder+config.graph_alias+".topk.pprs";
}

inline void save_exact_topk_ppr(){
    string filename = get_exact_topk_ppr_file();
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << exact_topk_pprs;
}

inline void load_exact_topk_ppr(){
    string filename = get_exact_topk_ppr_file();
    if(!exists_test(filename)){
        INFO("No exact topk ppr file", filename);
        return;
    }
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> exact_topk_pprs;

    INFO(exact_topk_pprs.size());
}

inline string get_idx_file_name(){
    string file_name;
    if(config.rmax_scale==1){
        if(config.opt)
            file_name = config.graph_location+"randwalks.idx.onehopopt";
        else
            file_name = config.graph_location+"randwalks.idx";
    }
    else{
        if(config.opt)
            file_name = config.graph_location+"randwalks."+to_string(config.rmax_scale)+".idx.onehopopt";
        else
            file_name = config.graph_location+"randwalks."+to_string(config.rmax_scale)+".idx";
    }
    
    return file_name;
}

inline string get_idx_info_name(){
    string file_name;
    if(config.rmax_scale==1){
        if(config.opt)
                file_name = config.graph_location+"randwalks.info.onehopopt";
        else
                file_name = config.graph_location+"randwalks.info";
    }
    else{
        if(config.opt)
                file_name = config.graph_location+"randwalks."+to_string(config.rmax_scale)+".info.onehopopt";
        else
                file_name = config.graph_location+"randwalks."+to_string(config.rmax_scale)+".info";
        
    }
    return file_name;   
}



bool sort_by_sec_descending(const std::pair<int, int> &left, const std::pair<int, int> &right) 
{
         return ( left.second > right.second || 
               ( !( right.second > left.second ) && left.first < right.first ) );
}




inline void deserialize_idx(){
    string file_name = get_idx_file_name();
    INFO("Index file name:",file_name);
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    ia >> rw_idx;

    file_name = get_idx_info_name();
    assert_file_exist("index file", file_name);
    std::ifstream info_ifs(file_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> rw_idx_info;
}

inline void serialize_idx(){
    std::ofstream ofs(get_idx_file_name());
    boost::archive::binary_oarchive oa(ofs);
    oa << rw_idx;

    std::ofstream info_ofs(get_idx_info_name());
    boost::archive::binary_oarchive info_oa(info_ofs);
    info_oa << rw_idx_info;
}

void single_build(const Graph& graph, int start, int end, vector<int>& rw_data, unordered_map<int, pair<unsigned long long, unsigned long> >& rw_info_map, int core_id){
    unsigned long num_rw;
    for(int v=start; v<end; v++){
        num_rw = ceil(graph.g[v].size()*config.rmax*config.omega);
        rw_info_map[v] = MP(rw_data.size(), num_rw);
        for(unsigned long i=0; i<num_rw; i++){
            int des = random_walk_thd(v, graph, core_id);
            rw_data.push_back(des);
        }
    }
}

void multi_build(const Graph& graph){
    INFO("multithread building...");
    fora_setting(graph.n, graph.m);

    // rw_idx = RwIdx( graph.n, vector<int>() );
    rw_idx_info.resize(graph.n);

    unsigned NUM_CORES = std::thread::hardware_concurrency();
    assert(NUM_CORES >= 2);

    INFO(NUM_CORES);

    unsigned long long rw_max_size = graph.m*config.rmax*config.omega;
    INFO(rw_max_size, rw_idx.max_size());

    rw_idx.reserve(rw_max_size);

    vector< vector<int> > vec_rw(NUM_CORES+1);
    vector< unordered_map<int, pair<unsigned long long, unsigned long> > > vec_rw_info(NUM_CORES+1);
    std::vector< std::future<void> > futures(NUM_CORES+1);

    int num_node_per_core = graph.n/(NUM_CORES+1);
    int start=0;
    int end=0;

    {
        INFO("rand-walking...");
        Timer tm(1);
        for(int core_id=0; core_id<NUM_CORES+1; core_id++){
            end = start + num_node_per_core;
            if(core_id==NUM_CORES)
                end = graph.n;
            
            vec_rw[core_id].reserve(rw_max_size/NUM_CORES);
            futures[core_id] = std::async( std::launch::async, single_build, std::ref(graph), start, end, std::ref(vec_rw[core_id]), std::ref(vec_rw_info[core_id]), core_id );
            start = end;
        }
        std::for_each( futures.begin(), futures.end(), std::mem_fn(&std::future<void>::wait));
    }

    {
        INFO("merging...");
        Timer tm(2);
        start=0;
        end=0;
        for(int core_id=0; core_id<NUM_CORES+1; core_id++){
            end = start + num_node_per_core;
            if(core_id==NUM_CORES)
                end = graph.n;

            rw_idx.insert( rw_idx.end(), vec_rw[core_id].begin(), vec_rw[core_id].end() );

            for(int v=start; v<end; v++){
                unsigned long long p = vec_rw_info[core_id][v].first;
                unsigned long num_rw = vec_rw_info[core_id][v].second;
                rw_idx_info[v] = MP( p + rw_idx.size()-vec_rw[core_id].size(), num_rw);
            }
            start = end;
        }
    }

    {
        INFO("materializing...");
        INFO(rw_idx.size(), rw_idx_info.size());
        Timer tm(3);
        serialize_idx(); //serialize the idx to disk
    }

    cout << "Memory usage (MB):" << get_proc_memory()/1000.0 << endl << endl;
}

void build(const Graph& graph){
    // size_t space = get_ram_size();
    // size_t estimated_space = sizeof(RwIdx) + graph.n *( sizeof(vector<int>) + config.num_rw*sizeof(int) );

    // if(estimated_space > space) //if estimated raw space overflows system maximum raw space, reset number of rand-walks
    //     config.num_rw = space * config.num_rw / estimated_space;
    
    fora_setting(graph.n, graph.m);

    fwd_idx.first.nil = -1;
    fwd_idx.second.nil =-1;
    fwd_idx.first.initialize(graph.n);
    fwd_idx.second.initialize(graph.n);

    // rw_idx = RwIdx( graph.n, vector<int>() );
    rw_idx_info.resize(graph.n);

    long long tuned_index_size =0;
    long long original_index_size =0;

    {
        Timer tm(1);
        unsigned long num_rw;
        for(int source=0; source<graph.n; source++){ //from each node, do rand-walks
        
            original_index_size += ceil(graph.g[source].size()*config.rmax*config.omega);
            if(config.opt)
                num_rw = ceil(graph.g[source].size()*config.rmax*(1-config.alpha)*config.omega);
            else
                num_rw = ceil(graph.g[source].size()*config.rmax*config.omega);
            rw_idx_info[source] = MP(tuned_index_size, num_rw);
            tuned_index_size += num_rw;
        }
    }
    INFO(tuned_index_size);

    rw_idx.reserve(tuned_index_size);

    {
        INFO("rand-walking...");
        INFO(config.rmax, config.omega, config.rmax*config.omega);
        Timer tm(1);
        for(int source=0; source<graph.n; source++){ //from each node, do rand-walks
            for(unsigned long i=0; i<rw_idx_info[source].second; i++){ //for each node, do some rand-walks
                int destination = 0;
                if(config.opt)
                    destination = random_walk_no_zero_hop(source, graph);
                else
                    destination = random_walk(source, graph);
                // rw_idx[source].push_back(destination);
                rw_idx.push_back(destination);
            }
        }
        INFO(original_index_size, tuned_index_size, original_index_size*1.0/tuned_index_size);
    }

    {
        INFO("materializing...");
        INFO(rw_idx.size(), rw_idx_info.size());
        Timer tm(2);
        serialize_idx(); //serialize the idx to disk
    }

    cout << "Memory usage (MB):" << get_proc_memory()/1000.0 << endl << endl;
}

#endif
