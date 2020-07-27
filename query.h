//Contributors: Sibo Wang, Renchi Yang
#ifndef FORA_QUERY_H
#define FORA_QUERY_H

#include "algo.h"
#include "graph.h"
#include "heap.h"
#include "config.h"
#include "build.h"

//#define CHECK_PPR_VALUES 1
// #define CHECK_TOP_K_PPR 1
#define PRINT_PRECISION_FOR_DIF_K 1
// std::mutex mtx;

void montecarlo_query(int v, const Graph& graph){
    Timer timer(MC_QUERY);

    rw_counter.clean();
    ppr.reset_zero_values();

    {
        Timer tm(RONDOM_WALK);
        num_total_rw += config.omega;
        for(unsigned long i=0; i<config.omega; i++){
            int destination = random_walk(v, graph);
            if(!rw_counter.exist(destination))
                rw_counter.insert(destination, 1);
            else
                rw_counter[destination] += 1;
        }
    }

    int node_id;
    for(long i=0; i<rw_counter.occur.m_num; i++){
        node_id = rw_counter.occur[i];
        ppr[node_id] = rw_counter[node_id]*1.0/config.omega;
    }

#ifdef CHECK_PPR_VALUES
    display_ppr();
#endif
}

void montecarlo_query_topk(int v, const Graph& graph){
    Timer timer(0);

    rw_counter.clean();
    ppr.clean();

    {
        Timer tm(RONDOM_WALK);
        num_total_rw += config.omega;
        for(unsigned long i=0; i<config.omega; i++){
            int destination = random_walk(v, graph);
            if(!rw_counter.exist(destination))
                rw_counter.insert(destination, 1);
            else
                rw_counter[destination] += 1;
        }
    }

    int node_id;
    for(long i=0; i<rw_counter.occur.m_num; i++){
        node_id = rw_counter.occur[i];
        if(rw_counter.occur[i]>0)
            ppr.insert( node_id, rw_counter[node_id]*1.0/config.omega );
    }
}

void bippr_query(int v, const Graph& graph){
    Timer timer(BIPPR_QUERY);

    rw_counter.clean();
    ppr.reset_zero_values();

    {
        Timer tm(RONDOM_WALK);
        num_total_rw += config.omega;
        INFO(config.omega);
        for(unsigned long i=0; i<config.omega; i++){
            int destination = random_walk(v, graph); 
            if(!rw_counter.exist(destination))
                rw_counter.insert(destination, 1);
            else
                rw_counter[destination] += 1;
        }
    }

    INFO(config.rmax);
    if(config.rmax < 1.0){
        Timer tm(BWD_LU);
        for(long i=0; i<graph.n; i++){
            reverse_local_update_linear(i, graph);
            // if(backresult.first[v] ==0 && backresult.second.size()==0){
            if( (!bwd_idx.first.exist(v)||0==bwd_idx.first[v]) &&  0==bwd_idx.second.occur.m_num){
                continue;
            }
            ppr[i] += bwd_idx.first[v];
            // for(auto residue: backresult.second){
            for(long j=0; j<bwd_idx.second.occur.m_num; j++){
                // ppr[i]+=counts[residue.first]*1.0/config.omega*residue.second;
                int nodeid = bwd_idx.second.occur[j];
                double residual = bwd_idx.second[nodeid];
                int occur;
                if(!rw_counter.exist(nodeid))
                    occur = 0;
                else
                    occur = rw_counter[nodeid]; 

                ppr[i] += occur*1.0/config.omega*residual;
            }
        }
    }else{
        int node_id;
        for(long i=0; i<rw_counter.occur.m_num; i++){
            node_id = rw_counter.occur[i];
            ppr[node_id] = rw_counter[node_id]*1.0/config.omega;
        }
    }
#ifdef CHECK_PPR_VALUES
    display_ppr();
#endif
}

void bippr_query_topk(int v, const Graph& graph){
    Timer timer(0);

    ppr.clean();
    rw_counter.clean();

    {
        Timer tm(RONDOM_WALK);
        num_total_rw += config.omega;
        for(unsigned long i=0; i<config.omega; i++){
            int destination = random_walk(v, graph);
            if(rw_counter.notexist(destination)){
                rw_counter.insert(destination, 1);
            }
            else{
                rw_counter[destination] += 1;
            }
        }
    }

    if(config.rmax < 1.0){
        Timer tm(BWD_LU);
        for(int i=0; i<graph.n; i++){
            reverse_local_update_linear(i, graph);
            if( (!bwd_idx.first.exist(v)||0==bwd_idx.first[v]) &&  0==bwd_idx.second.occur.m_num){
                continue;
            }

            if( bwd_idx.first.exist(v) && bwd_idx.first[v]>0 )
                ppr.insert(i, bwd_idx.first[v]);

            for(long j=0; j<bwd_idx.second.occur.m_num; j++){
                int nodeid = bwd_idx.second.occur[j];
                double residual = bwd_idx.second[nodeid];
                int occur;
                if(!rw_counter.exist(nodeid)){
                    occur = 0;
                }
                else{
                    occur = rw_counter[nodeid]; 
                }

                if(occur>0){
                    if(!ppr.exist(i)){
                        ppr.insert( i, occur*residual/config.omega );
                    }
                    else{
                        ppr[i] += occur*residual/config.omega;
                    }
                }
            }
        }
    }
    else{
        int node_id;
        for(long i=0; i<rw_counter.occur.m_num; i++){
            node_id = rw_counter.occur[i];
            if(rw_counter[node_id]>0){
                if(!ppr.exist(node_id)){
                    ppr.insert( node_id, rw_counter[node_id]*1.0/config.omega );
                }
                else{
                    ppr[node_id] = rw_counter[node_id]*1.0/config.omega;
                }
            }
        }
    }
}

void hubppr_query(int s, const Graph& graph){
    Timer timer(HUBPPR_QUERY);

    ppr.reset_zero_values();

    {
        Timer tm(RONDOM_WALK);
        fwd_with_hub_oracle(graph, s);
        count_hub_dest();
        INFO("finish fwd work", hub_counter.occur.m_num, rw_counter.occur.m_num);
    }

    {
        Timer tm(BWD_LU);
        for(int t=0; t<graph.n; t++){
            bwd_with_hub_oracle(graph, t);
            // reverse_local_update_linear(t, graph);
            if( (bwd_idx.first.notexist(s) || 0==bwd_idx.first[s]) && 0==bwd_idx.second.occur.m_num ){
                continue;
            }

            if(rw_counter.occur.m_num < bwd_idx.second.occur.m_num){ //iterate on smaller-size list
                for (int i=0; i<rw_counter.occur.m_num; i++) {
                    int node = rw_counter.occur[i];
                    if (bwd_idx.second.exist(node)) {
                        ppr[t] += bwd_idx.second[node]*rw_counter[node];
                    }
                }
            }
            else{
                for (int i=0; i<bwd_idx.second.occur.m_num; i++) {
                    int node = bwd_idx.second.occur[i];
                    if (rw_counter.exist(node)) {
                        ppr[t] += rw_counter[node]*bwd_idx.second[node];
                    }
                }
            }
            ppr[t]=ppr[t]/config.omega;
            if(bwd_idx.first.exist(s))
                ppr[t] += bwd_idx.first[s];
        }
    }

#ifdef CHECK_PPR_VALUES
    display_ppr();
#endif
}

void compute_ppr_with_reserve(){
    ppr.clean();
    int node_id;
    double reserve;
    for(long i=0; i< fwd_idx.first.occur.m_num; i++){
        node_id = fwd_idx.first.occur[i];
        reserve = fwd_idx.first[ node_id ];
        if(reserve)
            ppr.insert(node_id, reserve);
    }
}

void compute_ppr_with_fwdidx(const Graph& graph, double check_rsum){
    ppr.reset_zero_values();

    int node_id;
    double reserve;
    for(long i=0; i< fwd_idx.first.occur.m_num; i++){
        node_id = fwd_idx.first.occur[i];
        reserve = fwd_idx.first[ node_id ];
        ppr[node_id] = reserve;
    }

    // INFO("rsum is:", check_rsum);
    if(check_rsum == 0.0)
        return;

    unsigned long long num_random_walk = config.omega*check_rsum;
    // INFO(num_random_walk);
    //num_total_rw += num_random_walk;

    {
        Timer timer(RONDOM_WALK); //both rand-walk and source distribution are related with num_random_walk
        //Timer tm(SOURCE_DIST);
        if(config.with_rw_idx){
            fwd_idx.second.occur.Sort();
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                double residual = fwd_idx.second[source];
                unsigned long num_s_rw = ceil(residual/check_rsum*num_random_walk);
                double a_s = residual/check_rsum*num_random_walk/num_s_rw;

                double ppr_incre = a_s*check_rsum/num_random_walk;
                
                num_total_rw += num_s_rw;
                
                //for each source node, get rand walk destinations from previously generated idx or online rand walks
                if(num_s_rw > rw_idx_info[source].second){ //if we need more destinations than that in idx, rand walk online
                    for(unsigned long k=0; k<rw_idx_info[source].second; k++){
                        int des = rw_idx[rw_idx_info[source].first + k];
                        ppr[des] += ppr_incre;
                    }
                    num_hit_idx += rw_idx_info[source].second;

                    for(unsigned long j=0; j < num_s_rw-rw_idx_info[source].second; j++){ //rand walk online
                        int des = random_walk(source, graph);
                        ppr[des] += ppr_incre;
                    }
                }else{ // using previously generated idx is enough
                    for(unsigned long k=0; k<num_s_rw; k++){
                        int des = rw_idx[rw_idx_info[source].first + k];
                        ppr[des] += ppr_incre;
                    }
                    num_hit_idx += num_s_rw;
                }
            }
        }
        else{ //rand walk online
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                double residual = fwd_idx.second[source];
                unsigned long num_s_rw = ceil(residual/check_rsum*num_random_walk);
                double a_s = residual/check_rsum*num_random_walk/num_s_rw;

                double ppr_incre = a_s*check_rsum/num_random_walk;
                
                num_total_rw += num_s_rw;
                for(unsigned long j=0; j<num_s_rw; j++){
                    int des = random_walk(source, graph);
                    ppr[des] += ppr_incre;
                }
            }
        }
    }
}






void compute_ppr_with_fwdidx_opt(const Graph& graph, double check_rsum){
    ppr.reset_zero_values();

    int node_id;
    double reserve;
    for(long i=0; i< fwd_idx.first.occur.m_num; i++){
        node_id = fwd_idx.first.occur[i];
        reserve = fwd_idx.first[ node_id ];
        ppr[node_id] = reserve;
    }

    // INFO("rsum is:", check_rsum);
    if(check_rsum == 0.0)
        return;

    check_rsum*=(1-config.alpha);
    unsigned long long num_random_walk = config.omega*check_rsum;
    // INFO(num_random_walk);
    //num_total_rw += num_random_walk;

    {
        Timer timer(RONDOM_WALK); //both rand-walk and source distribution are related with num_random_walk
        //Timer tm(SOURCE_DIST);
        if(config.with_rw_idx){
            fwd_idx.second.occur.Sort();
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                //double residual = fwd_idx.second[source];
                if(!fwd_idx.second.exist(source)) continue;
                ppr[source]+=fwd_idx.second[source]*config.alpha;
                double residual = fwd_idx.second[source]*(1-config.alpha);


                unsigned long num_s_rw = ceil(residual/check_rsum*num_random_walk);
                double a_s = residual/check_rsum*num_random_walk/num_s_rw;

                double ppr_incre = a_s*check_rsum/num_random_walk;
                
                num_total_rw += num_s_rw;
                
                //for each source node, get rand walk destinations from previously generated idx or online rand walks
                if(num_s_rw > rw_idx_info[source].second){ //if we need more destinations than that in idx, rand walk online
                    for(unsigned long k=0; k<rw_idx_info[source].second; k++){
                        int des = rw_idx[rw_idx_info[source].first + k];
                        ppr[des] += ppr_incre;
                    }
                    num_hit_idx += rw_idx_info[source].second;

                    for(unsigned long j=0; j < num_s_rw-rw_idx_info[source].second; j++){ //rand walk online
                        int des = random_walk_no_zero_hop(source, graph);
                        ppr[des] += ppr_incre;
                    }
                }else{ // using previously generated idx is enough
                    for(unsigned long k=0; k<num_s_rw; k++){
                        int des = rw_idx[rw_idx_info[source].first + k];
                        ppr[des] += ppr_incre;
                    }
                    num_hit_idx += num_s_rw;
                }
            }
        }
        else{ //rand walk online
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                if(!fwd_idx.second.exist(source)) continue;
                ppr[source]+=fwd_idx.second[source]*config.alpha;
                double residual = fwd_idx.second[source]*(1-config.alpha);
                unsigned long num_s_rw = ceil(residual/check_rsum*num_random_walk);
                double a_s = residual/check_rsum*num_random_walk/num_s_rw;

                double ppr_incre = a_s*check_rsum/num_random_walk;
                num_total_rw += num_s_rw;
                for(unsigned long j=0; j<num_s_rw; j++){
                    int des = random_walk_no_zero_hop(source, graph);
                    ppr[des] += ppr_incre;
                }
            }
        }
    }
}




/*
void compute_ppr_with_fwdidx_opt_old(const Graph& graph, double& check_rsum){
    ppr.reset_zero_values();
    memset(destination_count, 0, graph.n*sizeof(int));
    int node_id;
    double reserve;
    double residue;
    if(check_rsum == 0.0)
        return;


    // INFO("rsum is:", check_rsum);

    check_rsum*=(1-config.alpha);
    unsigned long long num_random_walk = ceil(config.omega*check_rsum);
    INFO(num_random_walk);
    unsigned long long real_num_rand_walk = num_total_rw;
    //num_total_rw += num_random_walk;

    {
        Timer timer(RONDOM_WALK); //both rand-walk and source distribution are related with num_random_walk
        //Timer tm(SOURCE_DIST);
        if(config.with_rw_idx){
            fwd_idx.second.occur.Sort();
            INFO(fwd_idx.second.occur.m_num);
            for(int i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                if(fwd_idx.first.exist(source))
                    ppr[source]+= fwd_idx.first[source];
                if(!fwd_idx.second.exist(source)) continue;
                //int source = i;
                //int source = fwd_idx.second.occur[i];
                ppr[source]+=fwd_idx.second[source]*config.alpha;
                double residual = fwd_idx.second[source]*(1-config.alpha);
      
                //INFO(residual, pg_values[source]*config.rmax, residual/(pg_values[source]*config.rmax) );

                unsigned long num_s_rw = ceil(residual*config.omega);
                unsigned long num_s_int_rw = floor(residual*config.omega);
                num_total_rw += num_s_rw;
                
                //for each source node, get rand walk destinations from previously generated idx or online rand walks
                if(num_s_rw > rw_idx_info[source].second){ 
                    //if we need more destinations than that in idx, rand walk online
                    unsigned long k=0;
                    for(; k<rw_idx_info[source].second; k++){
                        int des = rw_idx[rw_idx_info[source].first + k];
                        destination_count[des]++;
                    }
                    num_hit_idx += rw_idx_info[source].second;

                    for(; k < num_s_int_rw; k++){ //rand walk online
                        int des = random_walk_no_zero_hop(source, graph);
                        destination_count[des]++;
                    }
                    if(num_s_int_rw< num_s_rw){
                        int des = random_walk_no_zero_hop(source, graph);
                        double ppr_incre=(1.0*residual*config.omega - num_s_int_rw*1.0)/config.omega;
                        ppr[des] += ppr_incre;
                    }
                }else{ 
                    // using previously generated idx is enough
                    unsigned long k=0;
                    for(; k<num_s_int_rw; k++){
                        int des = rw_idx[rw_idx_info[source].first + k];
                        destination_count[des] ++;
                    }
                    if(num_s_int_rw < num_s_rw){
                        int des = rw_idx[rw_idx_info[source].first + k];
                        double ppr_incre=(1.0*residual*config.omega - num_s_int_rw*1.0)/config.omega;
                        ppr[des] +=ppr_incre;
                    }
                    num_hit_idx += num_s_rw;
                }
            }
            for(long i=0; i< graph.n; i++){
                ppr[i]+= 1.0*destination_count[i]/config.omega;;
            }
        }
        else{ //rand walk online
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                double residual = fwd_idx.second[source];
                unsigned long num_s_rw = ceil(residual/check_rsum*num_random_walk);
                double a_s = residual/check_rsum*num_random_walk/num_s_rw;

                double ppr_incre = a_s*check_rsum/num_random_walk;
                
                num_total_rw += num_s_rw;
                for(unsigned long j=0; j<num_s_rw; j++){
                    int des = random_walk(source, graph);
                    ppr[des] += ppr_incre;
                }
            }
        }
        real_num_rand_walk = num_total_rw - real_num_rand_walk;
        INFO(real_num_rand_walk, num_random_walk);
    }
}*/




void compute_ppr_with_fwdidx_topk(const Graph& graph, double check_rsum){
    // ppr.clean();
    // // ppr.reset_zero_values();

    // int node_id;
    // double reserve;
    // for(long i=0; i< fwd_idx.first.occur.m_num; i++){
    //     node_id = fwd_idx.first.occur[i];
    //     reserve = fwd_idx.first[ node_id ];
    //     ppr.insert(node_id, reserve);
    //     // ppr[node_id] = reserve;
    // }
    compute_ppr_with_reserve();

    // INFO("rsum is:", check_rsum);
    if(check_rsum == 0.0)
        return;

    check_rsum*= (1-config.alpha);
    unsigned long long num_random_walk = config.omega*check_rsum;
    // INFO(num_random_walk);
    //num_total_rw += num_random_walk;
    {
        Timer timer(RONDOM_WALK); //both rand-walk and source distribution are related with num_random_walk
        //Timer tm(SOURCE_DIST);
        int source;
        double residual;
        unsigned long num_s_rw;
        double a_s;
        double ppr_incre;
        unsigned long num_used_idx;
        unsigned long num_remaining_idx;
        int des;
        //INFO(num_random_walk, fwd_idx.second.occur.m_num);
        if(config.with_rw_idx){ //rand walk with previously generated idx
            fwd_idx.second.occur.Sort();
            //for each source node, get rand walk destinations from previously generated idx or online rand walks
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                source = fwd_idx.second.occur[i];
                residual = fwd_idx.second[source];
                if(ppr.exist(source)){
                    ppr[source] += residual*config.alpha;
                }else{
                    ppr.insert(source, residual*config.alpha);
                }

                residual*=(1-config.alpha);
                num_s_rw = ceil(residual*config.omega);
                a_s = residual*config.omega/num_s_rw;

                ppr_incre = a_s/config.omega;

                num_total_rw += num_s_rw;

                num_used_idx = rw_counter[source];
                num_remaining_idx = rw_idx_info[source].second-num_used_idx;
                
                if(num_s_rw <= num_remaining_idx){
                    // using previously generated idx is enough
                    for(unsigned long k=0; k<num_s_rw; k++){
                        des = rw_idx[ rw_idx_info[source].first+ num_used_idx + k ];
                        if(!ppr.exist(des))
                            ppr.insert(des, ppr_incre);
                        else
                            ppr[des] += ppr_incre;
                    }

                    rw_counter[source] = num_used_idx + num_s_rw;

                    num_hit_idx += num_s_rw;
                }
                else{
                    //we need more destinations than that in idx, rand walk online
                    for(unsigned long k=0; k<num_remaining_idx; k++){
                        des = rw_idx[ rw_idx_info[source].first + num_used_idx + k ];
                        if(!ppr.exist(des))
                            ppr.insert(des, ppr_incre);
                        else
                            ppr[des] += ppr_incre;
                    }

                    num_hit_idx += num_remaining_idx;
                    rw_counter[source] = num_used_idx + num_remaining_idx;

                    for(unsigned long j=0; j < num_s_rw-num_remaining_idx; j++){ //rand walk online
                        des = random_walk_no_zero_hop(source, graph);
                        if(!ppr.exist(des))
                            ppr.insert(des, ppr_incre);
                        else
                            ppr[des] += ppr_incre;
                    }
                }
            }
        }
        else{ //rand walk online
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                source = fwd_idx.second.occur[i];
                residual = fwd_idx.second[source];
                num_s_rw = ceil(residual*config.omega);
                a_s = residual*config.omega/num_s_rw;

                ppr_incre = a_s/config.omega;
                num_total_rw += num_s_rw;

                for(unsigned long j=0; j<num_s_rw; j++){
                    des = random_walk(source, graph);
                    if(!ppr.exist(des))
                        ppr.insert(des, ppr_incre);
                    else
                        ppr[des] += ppr_incre;
                }
            }
        }
    }

}


void compute_ppr_with_fwdidx_topk_with_bound(const Graph& graph, double check_rsum){
    compute_ppr_with_reserve();

    if(check_rsum == 0.0)
        return;

    long num_random_walk = config.omega*check_rsum;
    long real_num_rand_walk=0;

    {
        Timer timer(RONDOM_WALK); //both rand-walk and source distribution are related with num_random_walk

        //Timer tm(SOURCE_DIST);
        if(config.with_rw_idx){ //rand walk with previously generated idx
            fwd_idx.second.occur.Sort();
            INFO(fwd_idx.second.occur.m_num);
            //for each source node, get rand walk destinations from previously generated idx or online rand walks
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                double residual = fwd_idx.second[source];
                long num_s_rw = ceil(residual*config.omega);
                double a_s = residual/check_rsum*num_random_walk/num_s_rw;

                double ppr_incre = a_s*check_rsum/num_random_walk;

                num_total_rw += num_s_rw;
                real_num_rand_walk += num_s_rw;

                long num_used_idx = 0;
                bool source_cnt_exist = rw_counter.exist(source);
                if( source_cnt_exist )
                    num_used_idx = rw_counter[source];
                long num_remaining_idx = rw_idx_info[source].second-num_used_idx;
                
                if(num_s_rw <= num_remaining_idx){
                    // using previously generated idx is enough
                    long k=0;
                    for(; k<num_remaining_idx; k++){
                        if( k < num_s_rw){
                            int des = rw_idx[rw_idx_info[source].first + k];
                            if(ppr.exist(des))
                                ppr[des] += ppr_incre;
                            else
                                ppr.insert(des, ppr_incre);
                        }else
                            break;
                    }

                    if(source_cnt_exist){
                        rw_counter[source] += k;
                    }
                    else{
                        rw_counter.insert(source, k);
                    }

                    num_hit_idx += k;
                }else{
                    //we need more destinations than that in idx, rand walk online
                    for(long k=0; k<num_remaining_idx; k++){
                        int des = rw_idx[ rw_idx_info[source].first + k ];
                        if(!ppr.exist(des))
                            ppr.insert(des, ppr_incre);
                        else
                            ppr[des] += ppr_incre;
                    }
                    num_hit_idx += num_remaining_idx;

                    if(!source_cnt_exist){
                        rw_counter.insert( source, num_remaining_idx );
                    }
                    else{
                        rw_counter[source] += num_remaining_idx;
                    }

                    for(long j=0; j < num_s_rw-num_remaining_idx; j++){ //rand walk online
                        int des = random_walk(source, graph);
                        if(!ppr.exist(des))
                            ppr.insert(des, ppr_incre);
                        else 
                            ppr[des] += ppr_incre;
                    }
                }
            }
        }
        else{ //rand walk online
            for(long i=0; i < fwd_idx.second.occur.m_num; i++){
                int source = fwd_idx.second.occur[i];
                double residual = fwd_idx.second[source];
                long num_s_rw = ceil(residual/check_rsum*num_random_walk);
                double a_s = residual/check_rsum*num_random_walk/num_s_rw;

                real_num_rand_walk += num_s_rw;

                double ppr_incre = a_s*check_rsum/num_random_walk;
                for(long j=0; j<num_s_rw; j++){
                    int des = random_walk(source, graph);
                    if(!ppr.exist(des))
                        ppr.insert(des, ppr_incre);
                    else
                        ppr[des] += ppr_incre;

                }
            }
        }
    }

    if(config.delta < threshold)
        set_ppr_bounds(graph, check_rsum, real_num_rand_walk);
    else{
        zero_ppr_upper_bound = calculate_lambda(check_rsum,  config.pfail, zero_ppr_upper_bound, real_num_rand_walk);
    }
}


void icompute_ppr_with_fwdidx(const Graph& graph, double check_rsum){
    //ppr.clean();
    int node_id;

    // INFO("rsum is:", check_rsum);
    if(check_rsum == 0.0)
        return;

    INFO("sampling random walk");

    {
        Timer timer(RONDOM_WALK); //both rand-walk and source distribution are related with num_random_walk
        //Timer tm(SOURCE_DIST);
        ifwd_idx.occur.Sort();
        if(config.with_rw_idx){

            INFO("Using index");
            INFO(ifwd_idx.occur.m_num);
            for(long i=0; i < ifwd_idx.occur.m_num; i++){
                int source = ifwd_idx.occur[i];
                unsigned long num_rw = ifwd_idx[source].second;
                if(num_rw == 0) continue;

                num_total_rw += num_rw;
                
                num_rw = min(num_rw, rw_idx_info[source].second);

                unsigned long start = rw_idx_info[source].first;

                for(unsigned long k=start; k<start+num_rw; k++){
                    int des = rw_idx[k];
                    if(!ifwd_idx.exist(des)){
                        ifwd_idx.insert(des, make_pair(1,0));
                    }else{
                        ifwd_idx[des].first++;
                    }
                }
                num_hit_idx += num_rw;

                    /*
                    for(unsigned long j=0; j < num_rw-rw_idx_info[source].second; j++){ //rand walk online
                        int des = random_walk(source, graph);

                        if(ifwd_idx.first.exist(des)){
                            ifwd_idx.first[des]++;
                        }else{
                            ifwd_idx.first.insert(des, 1);
                        }
                    }*/
            }
        }else{
            INFO("random walk online");
            for(long i=0; i < ifwd_idx.occur.m_num; i++){
                int source = ifwd_idx.occur[i];
                long num_s_rw = ifwd_idx[source].second;
                for(unsigned long j=0; j<num_s_rw; ++j){
                    int des = random_walk(source, graph);
                    if(ifwd_idx.exist(des)){
                        ifwd_idx[des].first++;
                    }else{
                        ifwd_idx.insert(des, make_pair(1, 0));
                    }
                }
            }
        }
    }
}

double total_rsum = 0.0;
double random_walk_time = 0.0000004;
double random_walk_index_time = random_walk_time/140;
double previous_rmax = 0;

double estimated_random_walk_cost(double rsum, double rmax){
    double estimated_random_walk_cost = 0.0;
    if(!config.with_rw_idx){
        estimated_random_walk_cost = config.omega*rsum*(1-config.alpha)*random_walk_time;
    }else{
        if(rmax >= config.rmax){
            estimated_random_walk_cost = config.omega*rsum*(1-config.alpha)*random_walk_time;
        }else{
            estimated_random_walk_cost = config.omega*rsum*(1-config.alpha)*random_walk_index_time;
        }
    }
    INFO(rmax, config.rmax, estimated_random_walk_cost);
    return estimated_random_walk_cost;
}
//map<double, int> count_ratio;
void fora_query_basic(int v, const Graph& graph){
    Timer timer(FORA_QUERY);
    double rsum = 1.0;

    {
        Timer timer(FWD_LU);

        if(config.balanced){
            static vector<int> forward_from;
            forward_from.clear();
            forward_from.reserve(graph.n);
            forward_from.push_back(v);

            fwd_idx.first.clean();  //reserve
            fwd_idx.second.clean();  //residual
            fwd_idx.second.insert( v, rsum );

            const static double min_delta = 1.0 / graph.n;

            const static double lowest_delta_rmax = config.opt?config.epsilon*sqrt(min_delta/3/graph.m/log(2/config.pfail))/(1-config.alpha):config.epsilon*sqrt(min_delta/3/graph.m/log(2/config.pfail));
            double used_time = 0;
            double rmax = 0;
            rmax = config.rmax*8;
            double random_walk_cost = 0;
            INFO(graph.g[v].size()>0);
            if(graph.g[v].size()>0){
                while(estimated_random_walk_cost(rsum, rmax)> used_time){
                    INFO(config.omega*rsum*random_walk_time, used_time);
                    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
                    forward_local_update_linear_topk( v, graph, rsum, rmax, lowest_delta_rmax, forward_from ); 
                    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - startTime).count();
                    INFO(rsum);
                    used_time +=duration/TIMES_PER_SEC;
                    double used_time_this_iteration = duration/TIMES_PER_SEC;
                    INFO(used_time_this_iteration);
                    rmax /=2;
                }
                rmax*=2;
                INFO("Adpaitve total forward push time: ", used_time);
                INFO(config.rmax, rmax, config.rmax/rmax);
                //count_ratio[config.rmax/rmax]++;
            }else{
                forward_local_update_linear(v, graph, rsum, config.rmax);
            }
        }else{
            forward_local_update_linear(v, graph, rsum, config.rmax); //forward propagation, obtain reserve and residual
        }
        //forward_local_update_linear_with_prune(v, graph, rsum, config.rmax); 
    }

    INFO(config.omega, config.omega*rsum, rsum);
    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
    //return;
    if(config.opt){
        compute_ppr_with_fwdidx_opt(graph, rsum);
        total_rsum+= rsum*(1-config.alpha);
    }else{
        compute_ppr_with_fwdidx(graph, rsum);
        total_rsum+= rsum*(1-config.alpha);
    }
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - startTime).count();
    INFO("Total random walk time: ", duration/TIMES_PER_SEC);

#ifdef CHECK_PPR_VALUES
    display_ppr();
#endif
}

void fora_query_topk_with_bound(int v, const Graph& graph){
    Timer timer(0);
    const static double min_delta = 1.0 / graph.n;
    const static double init_delta = 1.0 / 4;
    threshold = (1.0-config.ppr_decay_alpha)/pow(500, config.ppr_decay_alpha) / pow(graph.n, 1-config.ppr_decay_alpha);

    const static double new_pfail = 1.0 / graph.n / graph.n/log(graph.n);

    config.pfail = new_pfail;  // log(1/pfail) -> log(1*n/pfail)
    config.delta = init_delta;

    const static double lowest_delta_rmax = config.epsilon*sqrt(min_delta/3/graph.m/log(2/new_pfail));

    double rsum = 1.0;

    static vector<int> forward_from;
    forward_from.clear();
    forward_from.reserve(graph.n);
    forward_from.push_back(v);

    fwd_idx.first.clean();  //reserve
    fwd_idx.second.clean();  //residual
    fwd_idx.second.insert( v, rsum );

    zero_ppr_upper_bound = 1.0;

    if(config.with_rw_idx)
        rw_counter.reset_zero_values(); //store the pointers of each node's used rand-walk idxs 

    // for delta: try value from 1/4 to 1/n
    int iteration = 0;
    upper_bounds.reset_one_values();
    lower_bounds.reset_zero_values();

    while( config.delta >= min_delta ){
        fora_setting(graph.n, graph.m);
        num_iter_topk++;

        {
            Timer timer(FWD_LU);

            if(graph.g[v].size()==0){
                rsum = 0.0;
                fwd_idx.first.insert(v, 1);
                compute_ppr_with_reserve();
                return;
            }else{
                forward_local_update_linear_topk( v, graph, rsum, config.rmax, lowest_delta_rmax, forward_from ); //forward propagation, obtain reserve and residual
            }
        }

        compute_ppr_with_fwdidx_topk_with_bound(graph, rsum);
        {
            Timer timer(STOP_CHECK);
            if(if_stop() || config.delta <= min_delta){
                break;
            }else
                config.delta = max( min_delta, config.delta/2.0 );  // otherwise, reduce delta to delta/2
        }
    }
}


void fora_query_topk_new(int v, const Graph& graph ){
    Timer timer(0);
    const static double min_delta = 1.0 / graph.n;
    if(config.k ==0) config.k = 500;
    const static double init_delta = 1.0/config.k/10;//(1.0-config.ppr_decay_alpha)/pow(500, config.ppr_decay_alpha) / pow(graph.n, 1-config.ppr_decay_alpha);
    const static double new_pfail = 1.0 / graph.n / graph.n;

    config.pfail = new_pfail;  // log(1/pfail) -> log(1*n/pfail)
    config.delta = init_delta;

    const static double lowest_delta_rmax = config.epsilon*sqrt(min_delta/3/graph.m/log(2/new_pfail));

    double rsum = 1.0;

    static vector<int> forward_from;
    forward_from.clear();
    forward_from.reserve(graph.n);
    forward_from.push_back(v);

    fwd_idx.first.clean();  //reserve
    fwd_idx.second.clean();  //residual
    fwd_idx.second.insert( v, rsum );

    

    if(config.with_rw_idx)
        rw_counter.reset_zero_values(); //store the pointers of each node's used rand-walk idxs 

    // for delta: try value from 1/4 to 1/n
    while( config.delta >= min_delta ){
        fora_topk_setting(graph.n, graph.m);
        num_iter_topk++;
        {
            Timer timer(FWD_LU);
            //INFO(config.rmax, graph.m*config.rmax, config.omega);
            if(graph.g[v].size()==0){
                rsum = 0.0;
                fwd_idx.first.insert(v, 1);
                compute_ppr_with_reserve();
                return;
            }else{
                forward_local_update_linear_topk( v, graph, rsum, config.rmax, lowest_delta_rmax, forward_from ); //forward propagation, obtain reserve and residual
            }
        }

        //i_destination_count.clean();
        //compute_ppr_with_fwdidx_new(graph, rsum);
        //compute_ppr_with_fwdidx_topk(graph, rsum);
        compute_ppr_with_fwdidx_topk(graph, rsum);

        
        {
            double kth_ppr_score = kth_ppr();

            //topk_ppr();

            //double kth_ppr_score = topk_pprs[config.k-1].second;
            if( kth_ppr_score >= (1+config.epsilon)*config.delta || config.delta <= min_delta ){  // once k-th ppr value in top-k list >= (1+epsilon)*delta, terminate
                INFO(kth_ppr_score, config.delta, rsum);
                break;
            }
            else{
                /*int j=0;
                for(; j<config.k; j++){
                    //INFO(topk_pprs[j].second, (1+config.epsilon)*config.delta);
                    if(topk_pprs[j].second<(1+config.epsilon)*config.delta)
                        break;
                }
                INFO("Our current accurate top-j", j);*/
                config.delta = max( min_delta, config.delta/4.0 );  // otherwise, reduce delta to delta/4
            }
        }
    }
}



iMap<int> updated_pprs;
void hubppr_query_topk_martingale(int s, const Graph& graph) {
    unsigned long long the_omega = 2*config.rmax*log(2*config.k/config.pfail)/config.epsilon/config.epsilon/config.delta;
    static double bwd_cost_div = 1.0*graph.m/graph.n/config.alpha;

    static double min_ppr = 1.0/graph.n;
    static double new_pfail = config.pfail/2.0/graph.n/log2(1.0*graph.n*config.alpha*graph.n*graph.n);
    static double pfail_star = log(new_pfail/2);

    static std::vector<bool> target_flag(graph.n);
    static std::vector<double> m_omega(graph.n);
    static vector<vector<int>> node_targets(graph.n);
    static double cur_rmax=1;

    // rw_counter.clean();
    for(int t=0; t<graph.n; t++){
        map_lower_bounds[t].second = 0;//min_ppr;
        upper_bounds[t] = 1.0;
        target_flag[t] = true;
        m_omega[t]=0;
    }

    int num_iter = 1;
    int target_size=graph.n;
    if(cur_rmax>config.rmax){
        cur_rmax=config.rmax;
        for(int t=0; t<graph.n; t++){
            if(target_flag[t]==false)
                continue;
            reverse_local_update_topk(s, t, reserve_maps[t], cur_rmax, residual_maps[t], graph);
            for(const auto &p: residual_maps[t]){
                node_targets[p.first].push_back(t);
            }
        }
    }
    while( target_size > config.k && num_iter<=64 ){ //2^num_iter <= 2^64 since 2^64 is the largest unsigned integer here
        unsigned long long num_rw = pow(2, num_iter);
        rw_counter.clean();
        generate_accumulated_fwd_randwalk(s, graph, num_rw);
        updated_pprs.clean();
        // update m_omega
        {
            for(int x=0; x<rw_counter.occur.m_num; x++){
                int node = rw_counter.occur[x];
                for(const int t: node_targets[node]){
                    if(target_flag[t]==false)
                        continue;
                    m_omega[t] += rw_counter[node]*residual_maps[t][node];
                    if(!updated_pprs.exist(t))
                        updated_pprs.insert(t, 1);
                }
            }
        }

        double b = (2*num_rw-1)*pow(cur_rmax/2.0, 2);
        double lambda = sqrt(pow(cur_rmax*pfail_star/3, 2) - 2*b*pfail_star) - cur_rmax*pfail_star/3;
        {
            for(int i=0; i<updated_pprs.occur.m_num; i++){
                int t = updated_pprs.occur[i];
                if( target_flag[t]==false )
                    continue;

                double reserve = 0;
                if(reserve_maps[t].find(s)!=reserve_maps[t].end()){
                    reserve = reserve_maps[t][s];
                }
                set_martingale_bound(lambda, 2*num_rw-1, t, reserve, cur_rmax, pfail_star, min_ppr, m_omega[t]);
            }
        }

        topk_pprs.clear();
        topk_pprs.resize(config.k);
        partial_sort_copy(map_lower_bounds.begin(), map_lower_bounds.end(), topk_pprs.begin(), topk_pprs.end(), 
            [](pair<int, double> const& l, pair<int, double> const& r){return l.second > r.second;});

        double k_bound = topk_pprs[config.k-1].second;
        if( k_bound*(1+config.epsilon) >= upper_bounds[topk_pprs[config.k-1].first] || (num_rw >= the_omega && cur_rmax <= config.rmax) ){
            break;
        }

        for(int t=0; t<graph.n; t++){
            if(target_flag[t]==true && upper_bounds[t] <= k_bound){
                target_flag[t] = false;
                target_size--;
            }
        }
        num_iter++;
    }
}

void get_topk(int v, Graph &graph){
    display_setting();
    if(config.algo == MC){
        montecarlo_query_topk(v, graph);
        topk_ppr();
    }
    else if(config.algo == BIPPR){
        bippr_query_topk(v, graph);
        topk_ppr();
    }
    else if(config.algo == FORA){
        //fora_query_topk_new(v, graph);
        if(config.opt)
            fora_query_topk_new(v, graph);
        else
            fora_query_topk_with_bound(v, graph);
        topk_ppr();
    }
    else if(config.algo == FWDPUSH){
        Timer timer(0);
        double rsum = 1;
        
        {
            Timer timer(FWD_LU);
            forward_local_update_linear(v, graph, rsum, config.rmax);
        }
        compute_ppr_with_reserve();
        topk_ppr();
    }
    else if(config.algo == HUBPPR){
        Timer timer(0);
        hubppr_query_topk_martingale(v, graph);
    }

     // not FORA, so it's single source
     // no need to change k to run again
     // check top-k results for different k
    if(config.algo != FORA  && config.algo != HUBPPR){
        compute_precision_for_dif_k(v);
    }

    compute_precision(v);

#ifdef CHECK_TOP_K_PPR
    vector<pair<int, double>>& exact_result = exact_topk_pprs[v];
    INFO("query node:", v);
    for(int i=0; i<topk_pprs.size(); i++){
        cout << "Estimated k-th node: " << topk_pprs[i].first << " PPR score: " << topk_pprs[i].second << " " << map_lower_bounds[topk_pprs[i].first].first<< " " << map_lower_bounds[topk_pprs[i].first].second
             <<" Exact k-th node: " << exact_result[i].first << " PPR score: " << exact_result[i].second << endl;
    }
#endif
}

void fwd_power_iteration(const Graph& graph, int start, unordered_map<int, double>& map_ppr){
    static thread_local unordered_map<int, double> map_residual;
    map_residual[start] = 1.0;

    int num_iter=0;
    double rsum = 1.0;
    while( num_iter < config.max_iter_num ){
        num_iter++;
        // INFO(num_iter, rsum);
        vector< pair<int,double> > pairs(map_residual.begin(), map_residual.end());
        map_residual.clear();
        for(const auto &p: pairs){
            if(p.second > 0){
                map_ppr[p.first] += config.alpha*p.second;
                int out_deg = graph.g[p.first].size();

                double remain_residual = (1-config.alpha)*p.second;
                rsum -= config.alpha*p.second;
                if(out_deg==0){
                    map_residual[start] += remain_residual;
                }
                else{
                    double avg_push_residual = remain_residual / out_deg;
                    for (int next : graph.g[p.first]) {
                        map_residual[next] += avg_push_residual;
                    }
                }
            }
        }
        pairs.clear();
    }
    map_residual.clear();
}

void multi_power_iter(const Graph& graph, const vector<int>& source, unordered_map<int, vector<pair<int ,double>>>& map_topk_ppr ){
    static thread_local unordered_map<int, double> map_ppr;
    for(int start: source){
        fwd_power_iteration(graph, start, map_ppr);

        vector<pair<int ,double>> temp_top_ppr(config.k);
        partial_sort_copy(map_ppr.begin(), map_ppr.end(), temp_top_ppr.begin(), temp_top_ppr.end(), 
            [](pair<int, double> const& l, pair<int, double> const& r){return l.second > r.second;});
        
        map_ppr.clear();
        map_topk_ppr[start] = temp_top_ppr;
    }
}

void gen_exact_topk(const Graph& graph){
    // config.epsilon = 0.5;
    // montecarlo_setting();

    vector<int> queries;
    load_ss_query(queries);
    INFO(queries.size());
    unsigned int query_size = queries.size();
    query_size = min( query_size, config.query_size );
    INFO(query_size);
    string exact_top_file_str = get_exact_topk_ppr_file();
    if(exists_test(exact_top_file_str)){
        INFO("exact top k exists");
        return;
    }
    assert(config.k < graph.n-1);
    assert(config.k > 1);
    INFO(config.k);

    split_line();
    // montecarlo_setting();

    unsigned NUM_CORES = std::thread::hardware_concurrency()-1;
    assert(NUM_CORES >= 2);

    int num_thread = min(query_size, NUM_CORES);
    int avg_queries_per_thread = query_size/num_thread;

    vector<vector<int>> source_for_all_core(num_thread);
    vector<unordered_map<int, vector<pair<int ,double>>>> ppv_for_all_core(num_thread);

    for(int tid=0; tid<num_thread; tid++){
        int s = tid*avg_queries_per_thread;
        int t = s+avg_queries_per_thread;

        if(tid==num_thread-1)
            t+=query_size%num_thread;

        for(;s<t;s++){
            // cout << s+1 <<". source node:" << queries[s] << endl;
            source_for_all_core[tid].push_back(queries[s]);
        }
    }


    {
        Timer timer(PI_QUERY);
        INFO("power itrating...");
        std::vector< std::future<void> > futures(num_thread);
        for(int tid=0; tid<num_thread; tid++){
            futures[tid] = std::async( std::launch::async, multi_power_iter, std::ref(graph), std::ref(source_for_all_core[tid]), std::ref(ppv_for_all_core[tid]) );
        }
        std::for_each( futures.begin(), futures.end(), std::mem_fn(&std::future<void>::wait));
    }

    // cout << "average iter times:" << num_iter_topk/query_size << endl;
    cout << "average generation time (s): " << Timer::used(PI_QUERY)*1.0/query_size << endl;

    INFO("combine results...");
    for(int tid=0; tid<num_thread; tid++){
        for(auto &ppv: ppv_for_all_core[tid]){
            exact_topk_pprs.insert( ppv );
        }
        ppv_for_all_core[tid].clear();
    }

    save_exact_topk_ppr();
}

void topk(Graph& graph){
    vector<int> queries;
    load_ss_query(queries);
    INFO(queries.size());
    unsigned int query_size = queries.size();
    query_size = min( query_size, config.query_size );
    int used_counter=0;

    assert(config.k < graph.n-1);
    assert(config.k > 1);
    INFO(config.k);

    split_line();
    
    load_exact_topk_ppr();

     // not FORA, so it's single source
     // no need to change k to run again
     // check top-k results for different k
    if(config.algo != FORA && config.algo != HUBPPR){
        unsigned int step = config.k/5;
        if(step > 0){
            for(unsigned int i=1; i<5; i++){
                ks.push_back(i*step);
            }
        }
        ks.push_back(config.k);
        for(auto k: ks){	
            PredResult rst(0,0,0,0,0);
            pred_results.insert(MP(k, rst));
        }
    }

    used_counter = 0; 
    if(config.algo == FORA){
        fwd_idx.first.nil = -9;
        fwd_idx.first.initialize(graph.n);
        fwd_idx.second.nil = -9;
        fwd_idx.second.initialize(graph.n);
        rw_counter.nil =-9;
        rw_counter.init_keys(graph.n);
        upper_bounds.nil = -9;
        upper_bounds.init_keys(graph.n);
        lower_bounds.nil = -9;
        lower_bounds.init_keys(graph.n);
        ppr.nil = -9;
        ppr.initialize(graph.n);
        topk_filter.nil = -9;
        topk_filter.initialize(graph.n);
        //i_destination_count.nil = -9; 
        //i_destination_count.initialize(graph.n);
    }
    else if(config.algo == MC){
        rw_counter.initialize(graph.n);
        ppr.initialize(graph.n);
        montecarlo_setting();
    }
    else if(config.algo == BIPPR){
        bippr_setting(graph.n, graph.m);
        rw_counter.initialize(graph.n);
        bwd_idx.first.initialize(graph.n);
        bwd_idx.second.initialize(graph.n);
        ppr.initialize(graph.n); 
    }
    else if(config.algo == FWDPUSH){
        fwdpush_setting(graph.n, graph.m);
        fwd_idx.first.initialize(graph.n);
        fwd_idx.second.initialize(graph.n);
        ppr.initialize(graph.n);
    }
    else if(config.algo == HUBPPR){
        hubppr_topk_setting(graph.n, graph.m);
        rw_counter.initialize(graph.n);
        upper_bounds.init_keys(graph.n);
        if(config.with_rw_idx){
            hub_counter.initialize(graph.n);
            load_hubppr_oracle(graph);
        }
        residual_maps.resize(graph.n);
        reserve_maps.resize(graph.n);
        map_lower_bounds.resize(graph.n);
        for(int v=0; v<graph.n; v++){
            residual_maps[v][v]=1.0;
            map_lower_bounds[v] = MP(v, 0);
        }
        updated_pprs.initialize(graph.n);
    }

    for(int i=0; i<query_size; i++){
        cout << i+1 <<". source node:" << queries[i] << endl;
        get_topk(queries[i], graph);
        split_line();
    }

    cout << "average iter times:" << num_iter_topk/query_size << endl;
    display_time_usage(used_counter, query_size);
    set_result(graph, used_counter, query_size);

     //not FORA, so it's single source
     //no need to change k to run again
     // check top-k results for different k
    if(config.algo != FORA && config.algo != HUBPPR){
        display_precision_for_dif_k();
    }
}

void query(Graph& graph){
    INFO(config.algo);
    vector<int> queries;
    load_ss_query(queries);
    unsigned int query_size = queries.size();
    query_size = min( query_size, config.query_size );
    INFO(query_size);
    int used_counter=0;

    assert(config.rmax_scale >= 0);
    INFO(config.rmax_scale);

    ppr.init_keys(graph.n);

    if(config.algo == BIPPR){ //bippr
        bippr_setting(graph.n, graph.m);
        display_setting();
        used_counter = BIPPR_QUERY;

        bwd_idx.first.initialize(graph.n);
        bwd_idx.second.initialize(graph.n);

        rw_counter.initialize(graph.n);
        for(int i=0; i<query_size; i++){
            cout << i+1 <<". source node:" << queries[i] << endl;
            bippr_query(queries[i], graph);
            split_line();
        }
    }else if(config.algo == HUBPPR){
        bippr_setting(graph.n, graph.m);
        display_setting();
        used_counter = HUBPPR_QUERY;

        bwd_idx.first.initialize(graph.n);
        bwd_idx.second.initialize(graph.n);
        hub_counter.initialize(graph.n);
        rw_counter.initialize(graph.n);
        
        load_hubppr_oracle(graph);
        for(int i=0; i<query_size; i++){
            cout << i+1 <<". source node:" << queries[i] << endl;
            hubppr_query(queries[i], graph);
            split_line();
        }
    }
    else if(config.algo == FORA){ //fora
        fora_setting(graph.n, graph.m);
        display_setting();
        used_counter = FORA_QUERY;
        fwd_idx.first.nil = -1;
        fwd_idx.second.nil =-1;
        fwd_idx.first.initialize(graph.n);
        fwd_idx.second.initialize(graph.n);


        // rw_counter.initialize(graph.n);
        for(int i=0; i<query_size; i++){
            int source =  queries[i];
            cout << i+1 <<". source node:" << source << endl;
            fora_query_basic(source, graph);
            split_line();
        }
        double avg_rsum = total_rsum/query_size;
        INFO(avg_rsum*config.omega);
        /*for(auto x: count_ratio){
            INFO(x.first, x.second);
        }*/
    }else if(config.algo == MC){ //mc
        montecarlo_setting();
        display_setting();
        used_counter = MC_QUERY;

        rw_counter.initialize(graph.n);

        for(int i=0; i<query_size; i++){
            cout << i+1 <<". source node:" << queries[i] << endl;
            montecarlo_query(queries[i], graph);
            split_line();
        }
    }
    else if(config.algo == FWDPUSH){
        fwdpush_setting(graph.n, graph.m);
        display_setting();
        used_counter = FWD_LU;

        fwd_idx.first.initialize(graph.n);
        fwd_idx.second.initialize(graph.n);

        for(int i=0; i<query_size; i++){
            cout << i+1 <<". source node:" << queries[i] << endl;
            Timer timer(used_counter);
            double rsum = 1;
            forward_local_update_linear(queries[i], graph, rsum, config.rmax);
            compute_ppr_with_reserve();
            split_line();
        }
    }

    display_time_usage(used_counter, query_size);
    set_result(graph, used_counter, query_size);
}

void batch_topk(Graph& graph){
    vector<int> queries;
    load_ss_query(queries);
    INFO(queries.size());
    unsigned int query_size = queries.size();
    query_size = min( query_size, config.query_size );
    int used_counter=0;
    assert(config.k < graph.n-1);
    assert(config.k > 1);
    INFO(config.k);

    split_line();
    load_exact_topk_ppr();

    used_counter = 0; 
    if(config.algo == FORA){
        fwd_idx.first.nil = -9;
        fwd_idx.first.initialize(graph.n);
        fwd_idx.second.nil = -9;
        fwd_idx.second.initialize(graph.n);
        rw_counter.nil =-9;
        rw_counter.init_keys(graph.n);
        upper_bounds.nil = -9;
        upper_bounds.init_keys(graph.n);
        lower_bounds.nil = -9;
        lower_bounds.init_keys(graph.n);
        ppr.nil = -9;
        ppr.initialize(graph.n);
        topk_filter.nil = -9;
        topk_filter.initialize(graph.n);
        //i_destination_count.nil = -9; 
        //i_destination_count.initialize(graph.n);
    }
    else if(config.algo == MC){
        rw_counter.initialize(graph.n);
        ppr.initialize(graph.n);
        montecarlo_setting();
    }
    else if(config.algo == BIPPR){
        bippr_setting(graph.n, graph.m);
        rw_counter.initialize(graph.n);
        bwd_idx.first.initialize(graph.n);
        bwd_idx.second.initialize(graph.n);  
        ppr.initialize(graph.n); 
    }
    else if(config.algo == FWDPUSH){
        fwdpush_setting(graph.n, graph.m);
        fwd_idx.first.initialize(graph.n);
        fwd_idx.second.initialize(graph.n);
        ppr.initialize(graph.n);
    }
    else if(config.algo == HUBPPR){
        hubppr_topk_setting(graph.n, graph.m);
        rw_counter.initialize(graph.n); 
        upper_bounds.init_keys(graph.n);
        if(config.with_rw_idx){
            hub_counter.initialize(graph.n);
            load_hubppr_oracle(graph);
        }
        residual_maps.resize(graph.n);
        reserve_maps.resize(graph.n);
        map_lower_bounds.resize(graph.n);
        for(int v=0; v<graph.n; v++){
            residual_maps[v][v]=1.0;
            map_lower_bounds[v] = MP(v, 0);
        }
        updated_pprs.initialize(graph.n);
    }

    unsigned int step = config.k/5;
    if(step > 0){
        for(unsigned int i=1; i<5; i++){
            ks.push_back(i*step);
        }
    }
    ks.push_back(config.k);
    for(auto k: ks){
        PredResult rst(0,0,0,0,0);
        pred_results.insert(MP(k, rst));
    }

    // not FORA, so it's of single source
    // no need to change k to run again
    // check top-k results for different k
    if(config.algo != FORA && config.algo != HUBPPR ){
        for(int i=0; i<query_size; i++){
            cout << i+1 <<". source node:" << queries[i] << endl;
            get_topk(queries[i], graph);
            split_line();
        }

        display_time_usage(used_counter, query_size);
        set_result(graph, used_counter, query_size);

        display_precision_for_dif_k();
    }
    else{ // for FORA, when k is changed, run algo again
        for(unsigned int k: ks){
            config.k = k;
            INFO("========================================");
            INFO("k is set to be ", config.k);
            result.topk_recall=0;
            result.topk_precision=0;
            result.real_topk_source_count=0;
            Timer::clearAll();
            for(int i=0; i<query_size; i++){
                cout << i+1 <<". source node:" << queries[i] << endl;
                get_topk(queries[i], graph);
                split_line();
            }
            pred_results[k].topk_precision=result.topk_precision;
            pred_results[k].topk_recall=result.topk_recall;
            pred_results[k].real_topk_source_count=result.real_topk_source_count;

            cout << "k=" << k << " precision=" << result.topk_precision/result.real_topk_source_count 
                              << " recall=" << result.topk_recall/result.real_topk_source_count << endl;
            cout << "Average query time (s):"<<Timer::used(used_counter)/query_size<<endl;
            Timer::reset(used_counter);
        }

        // display_time_usage(used_counter, query_size);
        display_precision_for_dif_k();
    }
}

#endif //FORA_QUERY_H
