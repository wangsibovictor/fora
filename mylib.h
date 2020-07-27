#ifndef __MYLIB_H__
#define __MYLIB_H__


#include <iostream>
#include <set>
#include <list>
#include <sstream>
#include <cmath>
//#include <queue>
//#include <boost/lockfree/queue.hpp>
#include <fstream>
#include <string>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <climits>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
//#include <deque>
#include <queue>
//#include <chrono>
//#include <boost/atomic.hpp>
#include <atomic>
#include <thread>
#include <future>
#include <queue>
#include <mutex>
#include <sys/resource.h>
#include <condition_variable>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;
typedef unsigned int uint;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int, int> ipair;
typedef pair<double, double> dpair;
#define MP make_pair
#define F first
#define S second


#ifndef TIMES_PER_SEC
#define TIMES_PER_SEC (1.0e9)
#endif

typedef char int8;
typedef unsigned char uint8;
typedef long long int64;
typedef unsigned long long uint64;

#define SIZE(t) (int)(t.size())
#define ALL(t) (t).begin(), (t).end()
#define FOR(i, n) for(int (i)=0; (i)<((int)(n)); (i)++)


typedef pair<int, int> ipair;

const int MAXN = 1000000;

const int VectorDefaultSize=20;
const int TOPNUM = 1;


template <typename _T>
class iVector
{
public:
    unsigned int m_size;
    _T* m_data;
    unsigned int m_num;

    void free_mem()
    {
        delete[] m_data;
    }

    iVector()
    {
        //printf("%d\n",VectorDefaultSize);
        m_size = VectorDefaultSize;
        m_data = new _T[VectorDefaultSize];
        m_num = 0;
    }
    iVector( unsigned int n )
    {
        if ( n == 0 )
        {
            n = VectorDefaultSize;
        }
//      printf("iVector allocate: %d\n",n);
        m_size = n;
        m_data = new _T[m_size];
        m_num = 0;
    }
    void push_back( _T d )
    {
        if ( m_num == m_size )
        {
            re_allocate( m_size*2 );
        }
        m_data[m_num] = d ;
        m_num++;        
    }
    void push_back( const _T* p, unsigned int len )
    {
        while ( m_num + len > m_size )
        {
            re_allocate( m_size*2 );
        }
        memcpy( m_data+m_num, p, sizeof(_T)*len );
        m_num += len;
    }

    void re_allocate( unsigned int size )
    {
        if ( size < m_num )
        {
            return;
        }
        _T* tmp = new _T[size];
        memcpy( tmp, m_data, sizeof(_T)*m_num );
        m_size = size;
        delete[] m_data;
        m_data = tmp;
    }
    void Sort()
    {
        if ( m_num < 20 )
        {
            int k ;
            _T tmp;
            for ( int i = 0 ; i < m_num-1 ; ++i )
            {
                k = i ;
                for ( int j = i+1 ; j < m_num ; ++j )
                    if ( m_data[j] < m_data[k] ) k = j ;
                if ( k != i )
                {
                    tmp = m_data[i];
                    m_data[i] = m_data[k];
                    m_data[k] = tmp;
                }
            }
        }
        else sort( m_data, m_data+m_num );
    }
    void unique()
    {
        if ( m_num == 0 ) return;
        Sort();
        unsigned int j = 0;
        for ( unsigned int i = 0 ; i < m_num ; ++i )
            if ( !(m_data[i] == m_data[j]) )
            {
                ++j;
                if ( j != i ) m_data[j] = m_data[i];
            }
        m_num = j+1;
    }
    int BinarySearch( _T& data )
    {
        for ( int x = 0 , y = m_num-1 ; x <= y ; )
        {
            int p = (x+y)/2;
            if ( m_data[p] == data ) return p;
            if ( m_data[p] < data ) x = p+1;
            else y = p-1;
        }
        return -1;
    }
    void clean()
    {
        m_num = 0;
    }
    void assign( iVector& t )
    {
        m_num = t.m_num;
        m_size = t.m_size;
        delete[] m_data;
        m_data = t.m_data;
    }

    bool remove( _T& x )
    {
        for ( int l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;

            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memmove( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
            else if ( m_data[m] < x ) l = m+1;
            else r = m;
        }
        return false;
    }

    void sorted_insert( _T& x )
    {
        if ( m_num == 0 )
        {
            push_back( x );
            return;
        }

        if ( m_num == m_size ) re_allocate( m_size*2 );

        int l,r;

        for ( l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;
            if ( m_data[m] < x ) l = m+1;
            else r = m;
        }

        if ( l < m_num && m_data[l] == x )
        {
            //printf("Insert Duplicate....\n");
            //cout<<x<<endl;
    //      break;
        }
        else
        {
            if ( m_num > l )
            {
                memmove( m_data+l+1, m_data+l, sizeof(_T)*(m_num-l) );
            }
            m_num++;
            m_data[l] = x;
        }
    }

    bool remove_unsorted( _T& x )
    {
        for ( int m = 0 ; m < m_num ; ++m )
        {
            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memcpy( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
        }
        return false;
    }

    _T& operator[]( unsigned int i )
    {
        //if ( i < 0 || i >= m_num ) 
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[i];
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //close range check for [] in iVector if release

};




template <typename _T>
struct iMap
{   
    _T* m_data;
    int m_num;
    int cur;
    iVector<int> occur;
    _T nil; 
    iMap()
    {
        m_data = NULL;
        m_num = 0;
        //nil = std::make_pair((long)-9,(long)-9);
        //nil = 1073741834;
    }
    iMap(int size){
        initialize(size);
    }
    void free_mem()
    {
        delete[] m_data;
        occur.free_mem();
    }

    void initialize( int n )
    {
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        //nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = nil;
        cur = 0;
    }
    void clean()
    {
        for ( int i = 0 ; i < occur.m_num ; ++i )
        {
            m_data[occur[i]] = nil;
        }
        occur.clean();
        cur = 0;
    }
    
    //init keys 0-n, value as 0
    void init_keys(int n){
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        //nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i ){
            m_data[i] = 0;
            occur.push_back( i );
            cur++;
        }
    }
    //reset all values to be zero
    void reset_zero_values(){
        // for ( int i = 0 ; i < m_num ; ++i )
            // m_data[i] = 0.0;
        memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    void reset_one_values(){
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = 1.0;
        // memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    _T get( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap get out of range!!!\n");
        //  return -8;
        //}
        return m_data[p];
    }
    _T& operator[](  int p )
    {
        //if ( i < 0 || i >= m_num ) 
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[p];
    }
    void erase( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap get out of range!!!\n");
        //}
        m_data[p] = nil;
        cur--;
    }
    bool notexist( int p )
    {
        return m_data[p] == nil ;
    }
    bool exist( int p )
    {
        return !(m_data[p] == nil);
    }
    void insert( int p , _T d )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap insert out of range!!!\n");
        //}
        if ( m_data[p] == nil )
        {
            occur.push_back( p );
            cur++;
        }
        m_data[p] = d;
    }
    void inc( int p )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("inc some unexisted point\n");
        //}
        m_data[p]++;
    }
    void inc( int p , int x )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("inc some unexisted point\n");
        //}
        m_data[p] += x;
    }
    void dec( int p )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("dec some unexisted point\n" );
        //}
        m_data[p]--;
    }
    //close range check when release!!!!!!!!!!!!!!!!!!!!    
};



struct PendingQueue
{
    iVector<int> queue;
    iMap<int> pos;
    int point;

    PendingQueue(){
    }
    PendingQueue( int n )
    {
        point = 0;
        pos.initialize( n );
    }

    void clean()
    {
        queue.clean();
        point = 0;
        pos.clean();
    }

    void push_back( int x )
    {
        pos.insert( x , queue.m_num );
        queue.push_back( x );
    }

    int pop()
    {
        for ( ; point < queue.m_num ; point++ )
        {
            if ( pos.get(queue[point]) == point )
            {
                point++;
                return queue[point-1];
            }
        }
        return -1;
    }
};



static inline string &ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
    return s;
}

static inline string &rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
    return s;
}

static inline string &trim(string &s) { return ltrim(rtrim(s)); }

string __n_variable(string t, int n);

#define __expand_nv(x) __n_variable(t, x)<< t##x << " "


template<class T0>
void ___debug(string t, T0 t0, ostream &os) {
    os << __expand_nv(0);
}

template<class T0, class T1>
void ___debug(string t, T0 t0, T1 t1, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1);
}

template<class T0, class T1, class T2>
void ___debug(string t, T0 t0, T1 t1, T2 t2, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1) << __expand_nv(2);
}

template<class T0, class T1, class T2, class T3>
void ___debug(string t, T0 t0, T1 t1, T2 t2, T3 t3, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1) << __expand_nv(2) << __expand_nv(3);
}

template<class T0, class T1, class T2, class T3, class T4>
void ___debug(string t, T0 t0, T1 t1, T2 t2, T3 t3, T4 t4, ostream &os) {
    os << __expand_nv(0) << __expand_nv(1) << __expand_nv(2) << __expand_nv(3) << __expand_nv(4);
}

template<class T0>
void ___debug(string t, deque<T0> t0, ostream &os) {
    os << __n_variable(t, 0);
    FOR(i, SIZE(t0))os << t0[i] << " ";
}

template<class T0>
void ___debug(string t, vector<T0> t0, ostream &os) {
    os << __n_variable(t, 0);
    FOR(i, SIZE(t0))os << t0[i] << " ";
}

template<class T0, class T1>
void ___debug(string t, vector<pair<T0, T1> > t0, ostream &os) {
    os << __n_variable(t, 0);
    FOR(i, SIZE(t0))os << t0[i].F << "," << t0[i].S << " ";
}

#define RUN_TIME(...) { int64 t=rdtsc();  __VA_ARGS__; t=rdtsc()-t; cout<<  #__VA_ARGS__ << " : " << t/TIMES_PER_SEC <<"s"<<endl;  }

#ifdef HEAD_TRACE
#define TRACE(...) {{ ___debug( #__VA_ARGS__,  __VA_ARGS__,cerr); cerr<<endl;  } }
#define IF_TRACE(args) args
#define TRACE_LINE(...) { ___debug( #__VA_ARGS__,  __VA_ARGS__,cerr); cerr<<"                    \033[100D";  }
#define TRACE_SKIP(a, ...) { static int c=-1; c++; if(c%a==0)TRACE( __VA_ARGS__); }
#define TRACE_LINE_SKIP(a, ...) { static int c=-1; c++; if(c%a==0) TRACE_LINE(__VA_ARGS__);  }
#define TRACE_LINE_END(...) {cerr<<endl; }
ofstream __HEAD_H__LOG("log.txt");
#define TRACE_LOG(...) { __HEAD_H__LOG.close(); ofstream cerr("log.txt", ofstream::out|ofstream::app); ___debug( #__VA_ARGS__,  __VA_ARGS__, cerr); cerr<<endl;  }
#else
#define TRACE(...) ;
#define IF_TRACE(args) ;
#define TRACE_LINE(...) ;
#define TRACE_SKIP(a, ...) ;
#define TRACE_LINE_SKIP(a, ...) ;
#define TRACE_LINE_END(...) ;
#define TRACE_LOG(...) ;
#endif //HEAD_TRACE


//void setInfoFile(string s) { __HEAD_H_FOUT.open(s.c_str()); }

#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}
#define INFO(...) do {\
     ___debug( #__VA_ARGS__,  __VA_ARGS__,cout); cout<<endl; \
    } while(0)
#define ASSERTT(v, ...) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; INFO(__VA_ARGS__); exit(1);}}

template<class T>
string toStr(T t);


class Counter {
public:
    static int cnt[1000];
    int myid = 0;

    Counter(int id = 0) {
        myid = id;
        cnt[id]++;
    }

    void add(int x) {
        cnt[myid] += x;
    }

    ~Counter() {
    }

    static void show() {
        for (int i = 0; i < 1000; i++)
            if (cnt[i] > 0)
                INFO("Counter", i, cnt[i]);
    }
};


uint64 rdtsc();


class Timer {
public:
    static vector<double> timeUsed;
    static vector<string> timeUsedDesc;
    int id;
    std::chrono::steady_clock::time_point startTime;
    bool showOnDestroy;

    Timer(int id, string desc = "", bool showOnDestroy = false) {
        this->id = id;
        while ((int) timeUsed.size() <= id) {
            timeUsed.push_back(0);
            timeUsedDesc.push_back("");
        }
        timeUsedDesc[id] = desc;
        startTime = std::chrono::steady_clock::now();
        this->showOnDestroy = showOnDestroy;
    }

    static double used(int id) {
        return timeUsed[id] / TIMES_PER_SEC;
    }

    ~Timer() {
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - startTime).count();
        if (showOnDestroy) {
            cout << "time spend on " << timeUsedDesc[id] << ":" << duration / TIMES_PER_SEC << "s" << endl;
        }
        timeUsed[id] += duration;
    }

    static void show(bool debug = false) {
        if (debug) { TRACE("### Timer");
        }
        else {
            INFO("### Timer");
        }
        for (int i = 0; i < (int) timeUsed.size(); i++) {
            if (timeUsed[i] > 0) {
                char str[1000];
                sprintf(str, "%.6lf", timeUsed[i] / TIMES_PER_SEC);
                string s = str;
                if ((int) s.size() < 15) s = " " + s;
                char t[1000];
                memset(t, 0, sizeof t);
                sprintf(t, "%4d %s %s", i, s.c_str(), timeUsedDesc[i].c_str());
                if (debug) { TRACE(t);
                }
                else {
                    INFO(t);
                }
            }
        }
    }

    static void reset(int id){
        timeUsed[id] = 0;
    }

    static void clearAll() {
        timeUsed.clear();
        timeUsedDesc.clear();
    }
};


static vector<string> combine_args(int argc, char **argv) {
    vector<string> args;
    for (int i = 0; i < argc; ++i) {
        args.push_back(argv[i]);
    }
    return args;
}

static string to_str(double t) {
    stringstream ss;
    ss << t;
    return ss.str();
}

inline bool file_exists_test(const std::string &name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

static string replace(std::string str, const std::string &from, const std::string &to) {
    //original string will not be modified
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos)
        return str;
    str.replace(start_pos, from.length(), to);
    return str;
}

// unit is kilobyte (KB)
inline int get_proc_memory(){
    struct rusage r_usage;
    getrusage(RUSAGE_SELF,&r_usage);
    return r_usage.ru_maxrss;
}

const string Green = "\033[0;32m";
const string Reset = "\033[0m";
const string Red = "\033[0;31m";


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

static void program_start(int argc, char **argv) {

    cout << Green << "--------------start------------" << get_current_time_str() << Reset << endl;
    string combine = "";
    for (int i = 1; i < argc; i++) {
        combine += argv[i];
        combine += " ";
    }
    cout << Green << "args:" << combine << Reset << endl;
}

static void program_stop() {
    cout << Red << "--------------stop------------" << get_current_time_str() << Reset << endl;
    cout << endl;
    cout << endl;
    cout << endl;

}

#ifndef NDEBUG
#   define ASSERTMSG(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERTMSG(condition, message) do { } while (false)
#endif

#endif //__MYLIB_H__
