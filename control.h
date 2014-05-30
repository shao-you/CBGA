#define limit_distance 2001//map size
#define max_ext_set 5000

#define population_sz 300
#define generation 500
#define crossover_rate 1.0
#define mutation_rate 0.2
#define migration_rate 0.0
#define migration_interval 100
#define max_common_tours 15
#define RADIUS 50
#define R_NUM 100
#define RUN_TIMES_AVG 15
#define convergence 15
#define AFDX_LOW_BOUND 0.4 
// 首先我們要定義π 這個常數:
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

// 把度度量的和弧度量的轉換式寫為巨集，並利用 #ifdef 避免重複定義:
#ifndef DEGREEOF
#define DEGREEOF(a) ((a*180.0)/M_PI)
#endif
 
#ifndef RADIANOF
#define RADIANOF(a) ((a*M_PI)/180.0)
#endif
#include <iostream>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <climits>//INT_MAX, INT_MIN
#include <cmath>//sqrt, abs
#include <iomanip>//setw
#include <cstring>//memcpy
using namespace std; 