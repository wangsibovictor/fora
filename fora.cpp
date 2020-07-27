//Contributors: Sibo Wang, Renchi Yang
#define _CRT_SECURE_NO_DEPRECATE
#define HEAD_INFO

#include "mylib.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <set>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "graph.h"
#include "config.h"
#include "algo.h"
#include "query.h"
#include "build.h"

#include <boost/progress.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <chrono>


using namespace std::chrono;

using namespace boost;
using namespace boost::property_tree;

using namespace std;


string get_time_path() {
    using namespace boost::posix_time;
    auto tm = second_clock::local_time();
#ifdef WIN32
    return  "../../execution/" + to_iso_string(tm);
#else
    return parent_folder+FILESEP+"execution/" + to_iso_string(tm);
#endif
}

#include <boost/program_options.hpp>

namespace po = boost::program_options;


using namespace std;

int main(int argc, char *argv[]) {
    ios::sync_with_stdio(false);
    program_start(argc, argv);

    // this is main entry
    Saver::init();
    srand(time(NULL));
    config.graph_alias = "nethept";
    for (int i = 0; i < argc; i++) {
        string help_str = ""
                "fora query --algo <algo> [options]\n"
                "fora topk  --algo <algo> [options]\n"
                "fora batch-topk --algo <algo> [options]\n"
                "fora build [options]\n"
                "fora generate-ss-query [options]\n"
                "fora gen-exact-topk [options]\n"
                "fora\n"
                "\n"
                "algo: \n"
                "  bippr\n"
                "  montecarlo\n"
                "  fora\n"
                "  hubppr\n"
                "options: \n"
                "  --prefix <prefix>\n"
                "  --epsilon <epsilon>\n"
                "  --dataset <dataset>\n"
                "  --query_size <queries count>\n"
                "  --k <top k>\n"
                "  --with_idx\n"
                "  --hub_space <hubppr oracle space-consumption>\n"
                "  --exact_ppr_path <eaact-topk-pprs-path>\n"
                "  --rw_ratio <rand-walk cost ratio>\n"
                "  --result_dir <directory to place results>"
                "  --rmax_scale <scale of rmax>\n";

        if (string(argv[i]) == "--help") {
            cout << help_str << endl;
            exit(0);
        }
    }

    config.action = argv[1];
    cout << "action: " << config.action << endl;

    // init graph first
    for (int i = 0; i < argc; i++) {
        if (string(argv[i]) == "--prefix") {
            config.prefix = argv[i + 1];
        }
        if (string(argv[i]) == "--dataset") {
            config.graph_alias = argv[i + 1];
        }
    }

    for (int i = 0; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--algo") {
    		config.algo = string(argv[i + 1]);
        }
       
        else if (arg == "--epsilon") {
            config.epsilon = atof(argv[i + 1]);
            INFO(config.epsilon);
        }else if(arg == "--multithread"){
             config.multithread = true;
        }
        else if(arg == "--result_dir"){
            config.exe_result_dir = string(argv[i + 1]);
        }
        else if( arg == "--exact_ppr_path"){
            config.exact_pprs_folder = string(argv[i + 1]);
        }
        else if(arg == "--with_idx"){
            config.with_rw_idx = true;
        }else if(arg == "--rmax_scale"){
            config.rmax_scale = atof(argv[i+1]);
        }else if(arg == "--force-rebuild"){
            config.force_rebuild = true;
        }else if (arg == "--query_size"){
            config.query_size = atoi(argv[i+1]);
        }else if(arg == "--hub_space"){
            config.hub_space_consum = atoi(argv[i+1]);
        }else if (arg == "--version"){
            config.version = argv[i+1];
        }
        else if(arg == "--k"){
            config.k = atoi(argv[i+1]);
        }
        else if(arg == "--rw_ratio"){
            config.rw_cost_ratio = atof(argv[i + 1]);
        }
        else if (arg == "--prefix" || arg == "--dataset") {
            // pass
        }else if(arg == "--opt"){
            config.opt = 1;
        }else if(arg == "--balanced"){
            config.balanced = true;
        }

        else if (arg.substr(0, 2) == "--") {
            cerr << "command not recognize " << arg << endl;
            exit(1);
        }
    }

    INFO(config.version);
    vector<string> possibleAlgo = {BIPPR, FORA, FWDPUSH, MC, HUBPPR};
    
    INFO(config.action);

    srand (time(NULL));
    if (config.action == QUERY) {
        auto f = find(possibleAlgo.begin(), possibleAlgo.end(), config.algo);
        assert (f != possibleAlgo.end());
        if(f == possibleAlgo.end()){
            INFO("Wrong algo param: ", config.algo);
            exit(1);
        }

        config.graph_location = config.get_graph_folder();
        Graph graph(config.graph_location);
        INFO("load graph finish");
        init_parameter(config, graph);
   
        INFO("finished initing parameters");
        INFO(graph.n, graph.m);

        // if(config.multithread){
        //     init_multi_setting(graph.n);
        // }

        if(config.with_rw_idx)
            deserialize_idx();
        query(graph);
        INFO("finished query");
    }else if (config.action == GEN_SS_QUERY){
        config.graph_location = config.get_graph_folder();
        Graph graph(config.graph_location);
        INFO(graph.n, graph.m);

        generate_ss_query(graph.n);
    }else if (config.action == TOPK){
        auto f = find(possibleAlgo.begin(), possibleAlgo.end(), config.algo);
        assert (f != possibleAlgo.end());
        if(f == possibleAlgo.end()){
            INFO("Wrong algo param: ", config.algo);
            exit(1);
        }
        
        config.graph_location = config.get_graph_folder();
        Graph graph(config.graph_location);
        INFO("load graph finish");
        init_parameter(config, graph);
   
        if(config.exact_pprs_folder=="" || !exists_test(config.exact_pprs_folder))
            config.exact_pprs_folder = config.graph_location;

        INFO("finihsed initing parameters");
        INFO(graph.n, graph.m);

        // if(config.multithread){
        //     init_multi_setting(graph.n);
        // }

        if(config.with_rw_idx)
            deserialize_idx();
            
        topk(graph);
    }else if(config.action == BATCH_TOPK){
        auto f = find(possibleAlgo.begin(), possibleAlgo.end(), config.algo);
        assert(f != possibleAlgo.end());
        if(f == possibleAlgo.end()){
            INFO("Wrong algo param: ", config.algo);
            exit(1);
        }
        
        config.graph_location = config.get_graph_folder();
        Graph graph(config.graph_location);
        INFO("load graph finish");
        init_parameter(config, graph);
   
        if(config.exact_pprs_folder=="" || !exists_test(config.exact_pprs_folder))
            config.exact_pprs_folder = config.graph_location;

        INFO("finihsed initing parameters");
        INFO(graph.n, graph.m);

        if(config.with_rw_idx)
            deserialize_idx();
            
        batch_topk(graph);
    }else if (config.action == GEN_EXACT_TOPK){
        config.graph_location = config.get_graph_folder();
        Graph graph(config.graph_location);
        INFO("load graph finish");
        init_parameter(config, graph);

        if(config.exact_pprs_folder=="" || !exists_test(config.exact_pprs_folder))
            config.exact_pprs_folder = config.graph_location;
   
        INFO("finihsed initing parameters");
        INFO(graph.n, graph.m);

        gen_exact_topk(graph);
    }else if(config.action == BUILD){
        config.graph_location = config.get_graph_folder();
        Graph graph(config.graph_location);
        INFO("load graph finish");
        init_parameter(config, graph);


        INFO("finihsed initing parameters");
        INFO(graph.n, graph.m);

        {
            Timer tm(0);
            if(config.multithread)
                multi_build(graph);
            else
                build(graph);
        }
    }else{
        cerr << "sub command not regoznized" << endl;
        exit(1);
    }
    Timer::show();
    if(config.action == QUERY || config.action == TOPK){
        Counter::show();
        auto args = combine_args(argc, argv);
        Saver::save_json(config, result, args);
    }
    else

    program_stop();
    return 0;
}
