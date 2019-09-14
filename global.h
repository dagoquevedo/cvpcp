/*
 File       : Global.h
 Project    : Heuristic for CVPCP
 Developer  : Dagoberto Quevedo
 Date       : Mar 13 2013
 Institution: PISIS - FIME - UANL
 */

#ifndef CVPCP_Global_h
#define CVPCP_Global_h

//Library
#include <algorithm>
#include <fcntl.h>
#include <float.h>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <map>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <unistd.h>

#include "getCPUTime.c"
#include "getRSS.c"

using namespace std;

//Solution
typedef struct Tsolution {
    double                  f;
    bool                    F;
    vector<int>             A;
    vector<int>             B;
    vector<int>             J;
    vector<int>             P;
    vector< vector<int> >   X;
} solution;

//Global Variables
int             type;
int             id;
double          bestObj;
int             p;
int             n;
int             **d;
int             *w;
int             *s;
int             **cxy;
double          D_MAX;
vector<int>     V;
solution        result;
int             i_best;
unsigned int    seed;
float           alpha;
int             e;
int             r_max;
bool            __DEBUG__  = false;

double          time_start,time_stop;
double          euclidean(double, double, double, double);
void            print(solution);
void            save(solution);
void            finalize();
int             indexOf(vector<int>, int);
unsigned int    get_seed();
int             random_int(int, int);
float           random_float(float, float);
int             random_choice(vector<int>);

inline void     free_array2int(int);
inline int      **matrix2int(int, int);
inline void     free_matrix2int(int **,int);
inline double   **matrix2(int, int);
inline void     free_matrix2(double **,int);

vector<int>     setDifference(vector<int>,vector<int>);
vector<int>     setUnion(vector<int>,vector<int>);
vector<int>     setIntersection(vector<int>,vector<int>);

vector<pair<int, double> > sortMapDouble(map<int,double>);

int indexOf(vector<int> A, int value)
{
    double index = find(A.begin(), A.end(), value) - A.begin();
    return (int) (index < A.size() ? index : -1);
}

double euclidean(double x1, double x2, double y1, double y2)
{
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

unsigned int get_seed()
{
    
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__) || defined(__APPLE__)
    
    int randomData = open("/dev/random", O_RDONLY);
    int seed;
    ssize_t result = read(randomData, &seed, sizeof seed);
    if (result == sizeof seed)
    {}
    else
    {
        seed = (unsigned int)time(NULL);
    }
    close(randomData);
#else
    seed=(unsigned int)time(NULL);
#endif
    return (unsigned int)seed;
}

int random_int(int min, int max)
{
    return min + (rand() % ( max - min + 1 ) );
}

float random_float(float min, float max)
{
    return ((max - min) * ((float)rand()/RAND_MAX)) + min;
}

int random_choice(vector<int> A)
{
    int i = random_int(0, (int)A.size() - 1);
    return A[i];
}

void print(solution result)
{
    size_t      peakSize;
    
#if defined(__APPLE__)
    peakSize = getPeakRSS();
#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    peakSize = getCurrentRSS();
#endif
    float gap = -1;
    
    if(!result.F)
        result.f = -1;
    else
        gap = bestObj > -1 ? (((result.f - bestObj) / bestObj) * 100) : -1;
    
    if (__DEBUG__) {
        printf("%d\t%d\t%d\t%d\t%.2lf\t%.2f\t%.2f\t%.4lf\t%zu\t%d\t%d\t%d\t%.2f\t%d\t%u\n",
               type,id,n,p,bestObj,result.f,gap,(getCPUTime() - time_start),
               peakSize,result.F,i_best,r_max,alpha,e,seed);
    }
    else {
        printf("%d\t%d\t%d\t%d\t%.2lf\t%.2f\t%.2f\t%.4lf\t%zu\t%d\n",
               type,id,n,p,bestObj,result.f,gap,(time_stop - time_start),
               peakSize,result.F);
    }
}

inline void free_array2int(int *a)
{
    delete [] a;
}

//Integer
inline int **matrix2int(int up1, int up2)
{
    int    i,j;
    int **tmp;
    
    tmp = new int*[up1];
    for (i = 0; i < up1; i++)
    {
        tmp[i] = new int[up2];
        for(j = 0;j < up2; j++) tmp[i][j] = 0;
    }
    return(tmp);
}

inline void free_matrix2int(int **a,int up1)
{
    int   i;
    
    if (!a) return;
    for (i = 0; i< up1; i++) delete [] a[i];
    delete [] a;
}

//Double
inline double **matrix2(int up1, int up2)
{
    int    i;
    double **tmp;
    
    tmp = new double*[up1];
    for (i = 0; i < up1; ++i)
        tmp[i] = new double[up2];
    
    return(tmp);
}

inline void free_matrix2(double **a,int up1)
{
    int   i;
    
    if (!a) return;
    for (i = 0; i < up1; i++) delete [] a[i];
    delete [] a;
}

vector<int> setDifference(vector<int> A,vector<int> B)
{
    vector<int> result;
    
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    
    set_difference( A.begin(), A.end(), B.begin(), B.end(), back_inserter(result));
    return result;
}

vector<int> setUnion(vector<int> A,vector<int> B)
{
    vector<int> result;
    
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    
    set_union(A.begin(), A.end(), B.begin(), B.end(), back_inserter(result));
    return result;
}

vector<int> setIntersection(vector<int> A,vector<int> B)
{
    vector<int> result;
    
    sort(A.begin(), A.end());
    sort(B.begin(), B.end());
    
    set_intersection(A.begin(), A.end(), B.begin(), B.end(), back_inserter(result));
    return result;
}

void finalize()
{
    free_matrix2int(d,n);
    free_matrix2int(cxy,n);
    free_array2int(w);
    free_array2int(s);
    
    result.X.clear();
    result.P.clear();
    result.A.clear();
    result.B.clear();
    result.J.clear();
}

//Sort Map, unique values
template<typename A, typename B>
pair<B,A> flip_pair(const pair<A,B> &p)
{
    return pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
map<B,A> flip_map(const map<A,B> &src)
{
    map<B,A> dst;
    transform(src.begin(), src.end(), inserter(dst, dst.begin()), flip_pair<A,B>);
    return dst;
}

//Sort map, non-unique values
struct DoubleCompar {
    bool operator()(const pair<int, double> &lhs, const pair<int, double> &rhs) {
        return lhs.second < rhs.second;
    }
};

vector<pair<int, double> > sortMapDouble(map<int,double> A)
{
    vector<pair<int, double> > B(A.begin(),A.end());
    partial_sort(B.begin(), B.end(), B.end(), DoubleCompar());
    return B;
}

void save(const char * path_out, solution result)
{
    ofstream result_file;
    vector<int>::iterator j;
    int k;
    
    result_file.open (path_out);
    
    for (k = 0; k < p; ++k) {
        result_file << result.P[k];
        result_file << "|";
        for (j = result.X[k].begin(); j != result.X[k].end(); ++j) {
            result_file << *j;
            result_file << " ";
        }
        result_file << "\n";
    }
    
    result_file.close();
}

#endif
