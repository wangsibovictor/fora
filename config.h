#ifndef __CONFIG_H__
#define __CONFIG_H__


#ifdef WIN32
#define FILESEP "\\"
#else
#define FILESEP "/"
#endif

#include <boost/progress.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <unordered_map>
#include <list>
using namespace boost;
using namespace boost::property_tree;


const double ALPHA_DEFAULT = 0.2;
const int ALPHA_DENOMINATOR = 5;
const int ALPHA_NUMERATOR = 1;

const int NUM_OF_QUERY = 20;

const string QUERY = "query";
const string GEN_SS_QUERY = "generate-ss-query";
const string TOPK = "topk";
const string BUILD = "build";
const string GEN_EXACT_TOPK = "gen-exact-topk";
//const string ASSVERSION = "version";
const string BATCH_TOPK = "batch-topk";

const string BIPPR = "bippr";
const string FORA = "fora";
const string FWDPUSH = "fwdpush";
const string MC = "montecarlo";
const string HUBPPR = "hubppr";

const int MC_QUERY = 1;
const int BIPPR_QUERY = 2;
const int FORA_QUERY = 3;
const int HUBPPR_QUERY = 4;
const int FWD_LU = 5;
const int RONDOM_WALK = 6;
const int SOURCE_DIST = 7;
const int SORT_MAP = 8;
const int BWD_LU = 9;
const int PI_QUERY = 10;
const int STOP_CHECK=11;

const int DEG_NORM = 7;
const double RG_COST = 1;

// 0.2 for webstanford, 0.15 for dblp, 0.126 for pokec, 0.128 for LP, 0.097 for orkut 
// const double SG_PUSH_COST = 0.125;
// const double SG_RW_COST = 1.0;

const double SG_RW_COST = 8.0;

#ifdef WIN32
const string parent_folder = "../../";
#else
const string parent_folder = string("./") + FILESEP;
#endif

typedef pair<map<int, double>, map<int, double>> HubBwdidx;
//              pi               residual
typedef pair<double, HubBwdidx> HubBwdidxWithResidual;

// typedef pair<unordered_map<int, double>, unordered_map<int, double>> Bwdidx;
typedef pair<iMap<double>, iMap<double>> Bwdidx;
typedef pair<iMap<double>, iMap<double>> Fwdidx;
typedef iMap<pair<long, long> > IFwdidx;


typedef std::vector< std::vector<int> > RwIdx; //random walk idx

class Config {
public:
    string graph_alias;
    string graph_location;

    string action = ""; // query/generate index, etc..

    string prefix = "d:\\dropbox\\research\\data\\";

    string version = "vector";

    string exe_result_dir = parent_folder;

    string get_graph_folder() {
        return prefix + graph_alias + FILESEP;
    }
    bool multithread = false;
    bool with_rw_idx = false;
    bool opt = false;
    bool remap = false;
    bool force_rebuild = false;
    bool balanced = false;
    // int num_rw = 10;

    double omega; // 1/omega  omega = # of random walk
    double rmax; // identical to r_max

    unsigned int query_size = 1000;

    unsigned int max_iter_num = 100;

    double pfail = 0;
    double dbar = 0;
    double epsilon = 0;
    double delta = 0;

    unsigned int k = 500;
    double ppr_decay_alpha = 0.77;

    double rw_cost_ratio = 8.0;//8.0;

    double rmax_scale = 1;
    double multithread_param = 1.0;

    string algo;

    double alpha = ALPHA_DEFAULT;
    int alpha_numerator = 1;
    int alpha_denominator = 5;

    string exact_pprs_folder;

    unsigned int hub_space_consum = 1;

    ptree get_data() {
        ptree data;
        data.put("graph_alias", graph_alias);
        data.put("action", action);
        data.put("alpha", alpha);
        data.put("pfail", pfail);
        data.put("epsilon", epsilon);
        data.put("delta", delta);
        data.put("idx", with_rw_idx);
        // data.put("avg-idx-count", num_rw);
        data.put("k", k);
        data.put("rand-walk & push cost ratio", rw_cost_ratio);
        data.put("query-size", query_size);
        data.put("algo", algo);
        data.put("rmax", rmax);
        data.put("rmax-scale", rmax_scale);
        data.put("omega", omega);
        data.put("result-dir", exe_result_dir);
        return data;
    }
};

class Result {
public:
    int n;
    long long m;
    double avg_query_time;

    double total_mem_usage;
    double total_time_usage;

    double num_randwalk;
    double num_rw_idx_use;
    double hit_idx_ratio;

    double randwalk_time;
    double randwalk_time_ratio;

    double propagation_time;
    double propagation_time_ratio;

    double source_dist_time;
    double source_dist_time_ratio;

    double topk_sort_time;
    double topk_sort_time_ratio;

    // double topk_precision;
    double topk_max_abs_err;
    double topk_avg_abs_err;
    double topk_max_relative_err;
    double topk_avg_relative_err;
    double topk_precision;
    double topk_recall;
    // double topk_max_add_err;
    // double topk_avg_add_err;
    int real_topk_source_count;

    ptree get_data() {
        ptree data;
        data.put("n", n);
        data.put("m", m);
        data.put("avg query time(s/q)", avg_query_time);

        data.put("total memory usage(MB)", total_mem_usage);

        data.put("total time usage(s)", total_time_usage);
        data.put("total time on rand-walks(s)", randwalk_time);
        data.put("total time on propagation(s)", propagation_time);
        // data.put("total time on source distribution(s)", source_dist_time);
        data.put("total time on sorting top-k ppr(s)", topk_sort_time);


        data.put("total time ratio on rand-walks(%)", randwalk_time_ratio);
        data.put("total time ratio on propagation(%)", propagation_time_ratio);
        // data.put("total time ratio on source distribution(%)", source_dist_time_ratio);
        // data.put("total time ratio on sorting top-k ppr(%)", topk_sort_time_ratio);

        data.put("total number of rand-walks", num_randwalk);
        data.put("total number of rand-walk idx used", num_rw_idx_use);
        data.put("total usage ratio of rand-walk idx", hit_idx_ratio);

        // data.put("avg topk precision", topk_precision);
        data.put("topk max absolute error", topk_max_abs_err/real_topk_source_count);
        data.put("topk avg absolute error", topk_avg_abs_err/real_topk_source_count);
        data.put("topk max relative error", topk_max_relative_err/real_topk_source_count);
        data.put("topk avg relative error", topk_avg_relative_err/real_topk_source_count);
        data.put("topk precision", topk_precision/real_topk_source_count);
        data.put("topk recall", topk_recall/real_topk_source_count);
        return data;
    }

};

extern Config config;
extern Result result;

bool exists_test(const std::string &name);

void assert_file_exist(string desc, string name);

namespace Saver {
    static string get_current_time_str() {
        time_t rawtime;
        struct tm *timeinfo;
        char buffer[80];

        time(&rawtime);
        timeinfo = localtime(&rawtime);

        strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
        std::string str(buffer);

        return str;

    }

    static string get_time_path() {
        // using namespace boost::posix_time;
        // auto tm = second_clock::local_time();
        if(!boost::algorithm::ends_with(config.exe_result_dir, FILESEP))
            config.exe_result_dir += FILESEP;
        config.exe_result_dir += "execution/";
        if(!boost::filesystem::exists(config.exe_result_dir)){
            boost::filesystem::path dir(config.exe_result_dir);
            boost::filesystem::create_directories(dir);
        }

        string filename = config.graph_alias+"."+config.action+"."+config.algo;
        if(config.algo == "assppr")
            filename = filename + "." + to_string((int)config.rw_cost_ratio);

        string idx_flag = config.with_rw_idx?"with_idx":"without_idx";
        filename = filename+"."+idx_flag+".";
        filename += "k-"+to_string(config.k)+".";
        filename += "rmax-"+to_string(config.rmax_scale);

        return config.exe_result_dir + filename;
        // return config.exe_result_dir + to_iso_string(tm);
    }

    static ptree combine;

    static void init() {
        combine.put("start_time", get_current_time_str());
    }


    static void save_json(Config &config, Result &result, vector<string> args) {
        ofstream fout(get_time_path() + ".json");
        string command_line = "";
        for (int i = 1; i < args.size(); i++) {
            command_line += " " + args[i];
        }
        combine.put("end_time", get_current_time_str());
        combine.put("command_line", command_line);
        combine.put_child("config", config.get_data());
        combine.put_child("result", result.get_data());
        ptree timer;
        for (int i = 0; i < (int) Timer::timeUsed.size(); i++) {
            if (Timer::timeUsed[i] > 0) {
                timer.put(to_str(i), Timer::timeUsed[i] / TIMES_PER_SEC);
            }
        }
        combine.put_child("timer", timer);

        write_json(fout, combine, true);
    }
};


#endif
