#include "mylib.h"


template<class T>
string toStr(T t) {
    stringstream ss;
    ss << t;
    return ss.str();
}

string __n_variable(string t, int n) {
    t = t + ',';
    int i = 0;
    if (n) for (; i < SIZE(t) && n; i++) if (t[i] == ',') n--;
    n = i;
    for (; t[i] != ','; i++);
    t = t.substr((unsigned long) n, (unsigned long) (i - n));
    trim(t);
    if (t[0] == '"') return "";
    return t + "=";
}

int Counter::cnt[1000] = {0};

vector<double> Timer::timeUsed;
vector<string> Timer::timeUsedDesc;


#ifndef WIN32
#ifdef __CYGWIN__

//CYGWIN
uint64 rdtsc() {
    uint64 t0;
    asm volatile("rdtsc" : "=A"(t0));
    return t0;
}

#else
//LINUX
uint64 rdtsc(void)
{
    unsigned a, d;
    //asm("cpuid");
    asm volatile("rdtsc" : "=a" (a), "=d" (d));
    return (((uint64)a) | (((uint64)d) << 32));
}
#endif
#endif


