#include "mylib.h"
#include "config.h"


Config config;
Result result;

bool exists_test(const std::string &name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    }
    else {
        f.close();
        return false;
    }
}

void assert_file_exist(string desc, string name) {

    if (!exists_test(name)) {
        cerr << desc << " " << name << " not find " << endl;
        exit(1);
    }
}

