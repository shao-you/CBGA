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
// �����ڭ̭n�w�q�k �o�ӱ`��:
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

// ��׫׶q���M���׶q���ഫ���g�������A�çQ�� #ifdef �קK���Ʃw�q:
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