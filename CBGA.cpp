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
ofstream footprint;
ofstream print_waypoint;

//joint-cluster(依直線二分法判斷是否交集) + spanning
struct node_list
{
    float* set_node_x;
    float* set_node_y;
    int set_node_sz;
};
struct chromosome
{
    int* chromo;
    float fitness_val;
	float rotate_degree;
	float max_step_len;
};
struct node
{
    node(int x, node* y):sensor(x),next(y){}
    int sensor;
    node* next;
};
struct info
{
    info(float m, float n, float w, float x, int y, node* z):
        set_x(m),set_y(n),l_line(w),r_line(x),set_size(y),ptr(z){}
    float set_x;
    float set_y;
    float l_line;
    float r_line;
    int set_size;
    node* ptr;
};
struct coordinate
{
    int x_axis;
    int y_axis;
    info* ptr;
	int r;//radius
};
struct X_Y
{
    float x_axis;
    float y_axis;
    int size;
    int* ary;//成員sensors
};

void random(coordinate* ary, int RandomNode)
{
print_waypoint<<"Print out sensor points:"<<endl;
	//srand(time(NULL));//set seed
    for(int i=0;i<RandomNode;i++)
    {
        int x = rand()%limit_distance;//0~2000
        int y = rand()%limit_distance;
print_waypoint<<" ("<<x<<", "<<y<<")";
        ary[i].x_axis = x;
        ary[i].y_axis = y;
        ary[i].ptr = NULL;
		ary[i].r = RADIUS;//--------1
		//int rr = rand()%RADIUS+1;
		//ary[i].r = rr;//------------2
          //cout<<"=="<<rr<<"==";
          //cout<<"("<<x<<", "<<y<<") ";
    }
print_waypoint<<endl;
}

inline double Euclideandistance(double a_x, double a_y, double b_x, double b_y)
{
    return sqrt((a_x-b_x)*(a_x-b_x) + (a_y-b_y)*(a_y-b_y));
}

void calculate_node_dist(coordinate* ary, float** vector, int RandomNode)
{
    int col = 1, row;
    for(int i=0;i<RandomNode-1;i++)
    {
        row = 0;
        for(int j=0;j<col;j++)
        {
            vector[col][row] = vector[row][col] = Euclideandistance(ary[col].x_axis,ary[col].y_axis,ary[row].x_axis,ary[row].y_axis);
              //cout<<vector[col][row]<<"   ";
            row++;
        }
        col++;
    }
}
void free_memory2(float** dist_vector_joint, X_Y* joint_vector, int joint_count)
{
    for(int i=0;i<joint_count;i++)
    {
        delete [] joint_vector[i].ary;
        delete [] dist_vector_joint[i];
    }
    delete [] joint_vector;
    delete [] dist_vector_joint;
}
void free_memory(coordinate* ary, float** vector, int RandomNode, float** dist_vector,
                int set_count, node_list* cache)
{
    for(int i=0;i<RandomNode;i++)
    {
        if(ary[i].ptr != NULL)
        {
            node* tmp = ary[i].ptr->ptr;
            delete ary[i].ptr;
            while(tmp != NULL)
            {
                node* tmp2 = tmp;
                tmp = tmp->next;
                ary[tmp2->sensor].ptr = NULL;
                delete tmp2;
            }
        }
        delete [] vector[i];
    }
    delete [] ary;
    delete [] vector;

    for(int i=0;i<set_count;i++)
    {
        delete [] dist_vector[i];
        delete [] cache[i].set_node_x;
        delete [] cache[i].set_node_y;
    }
    delete [] dist_vector;
    delete [] cache;
}

float find_int(int a, int b, float line, int Radius)//給圓心和line, 算出交點
{
    float z = Radius*Radius - (line-a)*(line-a);
    if(z < 0) return -1;//無解
	else if(z == 0) return 0;//相切
    else return sqrt(z);
}

bool check_set(coordinate* ary, info* current_set, int which_node)
{
    int ary_sz = (current_set->set_size+1);
    float upper[ary_sz];
    float lower[ary_sz];
    int count;

    for(float line=current_set->l_line ; line<=current_set->r_line ; line+=0.1)
    {
        //==========clear upper and lower
        count = 0;
        //==========start to find a line
        float result = find_int(ary[which_node].x_axis,ary[which_node].y_axis,line,ary[which_node].r);
        if(result == -1) continue;
        else
        {
            upper[count] = ary[which_node].y_axis + result;
            lower[count] = ary[which_node].y_axis - result;
            count++;
        }
        node* tmp = current_set->ptr;
        while(tmp != NULL)
        {
            result = find_int(ary[tmp->sensor].x_axis,ary[tmp->sensor].y_axis,line,ary[tmp->sensor].r);
            if(result == -1) break;//impossible to be a cluster
            else
            {
                upper[count] = ary[tmp->sensor].y_axis + result;
                lower[count] = ary[tmp->sensor].y_axis - result;
                count++;
                tmp = tmp->next;
            }
        }
        //==========compare upper and lower
    fail:
        if(count != ary_sz) continue;
        else
        {
            float min_upper = INT_MAX, max_lower = INT_MIN;
            for(int i=0;i<ary_sz;i++)
            {
                if(upper[i] < min_upper) min_upper = upper[i];
                if(lower[i] > max_lower) max_lower = lower[i];
            }
			if(max_lower > min_upper)
            {
                count = 0;
                goto fail;
            }
            current_set->set_x = line;
            current_set->set_y = (min_upper+max_lower)/2;
            return true;
        }
    }
    return false;
}

bool check_ext_set(int* s1, int p, int* s2, int q, int set1, int set2,
                   X_Y* sets_vector, X_Y* ext_set, int ext_count, int Radius, float** vector, coordinate* ary)
{
      //cout<<p<<","<<q<<" = ";
    //==========檢查s1和s2中sensors是否兩兩相交
    for(int g=0;g<p;g++) for(int h=0;h<q;h++)
    {
        if(vector[sets_vector[set1].ary[s1[g]]][sets_vector[set2].ary[s2[h]]] > 2*Radius) return false;
    }
    //==========有相交,檢查是否為subset
        for(int s=0;s<ext_count;s++)
        {
                int find_1 = 0, find_2 = 0;
                for(int g=0;g<p;g++)
                {
                    for(int e=0;e<ext_set[s].size;e++)
                    if(sets_vector[set1].ary[s1[g]] == ext_set[s].ary[e])
                    {
                        find_1++;
                        break;
                    }
                    if(find_1 != g+1) break;
                }
                if(find_1 != p) continue;
                for(int h=0;h<q;h++)
                {
                    for(int e=0;e<ext_set[s].size;e++)
                    if(sets_vector[set2].ary[s2[h]] == ext_set[s].ary[e])
                    {
                        find_2++;
                        break;
                    }
                    if(find_2 != h+1) break;
                }
                if(find_2 != q) continue;
                else return false;//already has a superset
        }
    //==========spacial case to speed up
    if(p == 1 && q == 1)
    {
        ext_set[ext_count].x_axis =
            (float)(ary[sets_vector[set1].ary[s1[0]]].x_axis + ary[sets_vector[set2].ary[s2[0]]].x_axis)/2;
        ext_set[ext_count].y_axis =
            (float)(ary[sets_vector[set1].ary[s1[0]]].y_axis + ary[sets_vector[set2].ary[s2[0]]].y_axis)/2;
        return true;
    }
    //==========決定line區間, 至少含兩個sensor以上的set, 取前兩個sensor決定interval
    float l_line, r_line;
    if(p >= 2)
    {
        if(ary[sets_vector[set1].ary[s1[0]]].x_axis < ary[sets_vector[set1].ary[s1[1]]].x_axis)
        {
            r_line = ary[sets_vector[set1].ary[s1[0]]].x_axis + Radius;
            l_line = ary[sets_vector[set1].ary[s1[1]]].x_axis - Radius;
        }
        else
        {
            r_line = ary[sets_vector[set1].ary[s1[1]]].x_axis + Radius;
            l_line = ary[sets_vector[set1].ary[s1[0]]].x_axis - Radius;
        }
    }
    else //if(q >= 2)
    {
        if(ary[sets_vector[set2].ary[s2[0]]].x_axis < ary[sets_vector[set2].ary[s2[1]]].x_axis)
        {
            r_line = ary[sets_vector[set2].ary[s2[0]]].x_axis + Radius;
            l_line = ary[sets_vector[set2].ary[s2[1]]].x_axis - Radius;
        }
        else
        {
            r_line = ary[sets_vector[set2].ary[s2[1]]].x_axis + Radius;
            l_line = ary[sets_vector[set2].ary[s2[0]]].x_axis - Radius;
        }
    }
    if(l_line < 0) l_line = 0;
    if(r_line > limit_distance-1) r_line = limit_distance-1;
    //==========進行直線二分法
    int ary_sz = p+q;
    float upper[ary_sz];
    float lower[ary_sz];
    int count;
    for(float line=l_line ; line<=r_line ; line+=0.1)
    {
        //==========clear upper and lower
        count = 0;
        //==========start to find a line
        for(int g=0;g<p;g++)
        {
            float result = find_int(ary[sets_vector[set1].ary[s1[g]]].x_axis,ary[sets_vector[set1].ary[s1[g]]].y_axis,line,Radius);
            if(result == -1) break;
            else
            {
                upper[count] = ary[sets_vector[set1].ary[s1[g]]].y_axis + result;
                lower[count] = ary[sets_vector[set1].ary[s1[g]]].y_axis - result;
                count++;
            }
        }

        if(count != p) continue;
        for(int h=0;h<q;h++)
        {
            float result = find_int(ary[sets_vector[set2].ary[s2[h]]].x_axis,ary[sets_vector[set2].ary[s2[h]]].y_axis,line,Radius);
            if(result == -1) break;
            else
            {
                upper[count] = ary[sets_vector[set2].ary[s2[h]]].y_axis + result;
                lower[count] = ary[sets_vector[set2].ary[s2[h]]].y_axis - result;
                count++;
            }
        }
        //==========compare upper and lower
    fail:
        if(count != ary_sz) continue;
        else
        {
            float min_upper = INT_MAX, max_lower = INT_MIN;
            for(int i=0;i<ary_sz;i++)
            {
                if(upper[i] < min_upper) min_upper = upper[i];
                if(lower[i] > max_lower) max_lower = lower[i];
                for(int j=0;j<ary_sz;j++)
                {
                    if(upper[i] < lower[j])
                    {
                        count = 0;
                        goto fail;
                    }
                }
            }
            ext_set[ext_count].x_axis = line;
            ext_set[ext_count].y_axis = (min_upper+max_lower)/2;
            return true;
        }
    }
    return false;
}
int cluster(coordinate* ary, float** vector, int RandomNode)
{
    int col = 0, row, set_count = RandomNode;
    for(int i=0;i<RandomNode-1;i++)
    {
        row = col+1;
        for(int j=0;j<RandomNode-col-1;j++)
        {
            if(vector[col][row] <= ary[col].r + ary[row].r)
            {
                if(ary[col].ptr == NULL && ary[row].ptr == NULL)
                {
                    //float xx = (float)(ary[col].x_axis + ary[row].x_axis)/2;
                    //float yy = (float)(ary[col].y_axis + ary[row].y_axis)/2;
                    ary[col].ptr = new info(0,0,0,0,1,new node(col,NULL));

                    float right, left;
                    if(ary[col].x_axis < ary[row].x_axis)
                    {
                        right = ary[col].x_axis + ary[col].r;
                        left = ary[row].x_axis - ary[row].r;
                    }
                    else
                    {
                        right = ary[row].x_axis + ary[row].r;
                        left = ary[col].x_axis - ary[col].r;
                    }
                    //set interval
                    if(left < 0) ary[col].ptr->l_line = 0;
                    else ary[col].ptr->l_line = left;
                    if(right > limit_distance-1) ary[col].ptr->r_line = limit_distance-1;
                    else ary[col].ptr->r_line = right;

					check_set(ary,ary[col].ptr,row);//set point is done
					ary[col].ptr->ptr = new node(row,ary[col].ptr->ptr);
					ary[col].ptr->set_size++;
                    ary[row].ptr = ary[col].ptr;
                    set_count--;
                }
                else if(ary[col].ptr != NULL && ary[row].ptr == NULL)
                {
                    bool Y_N = check_set(ary,ary[col].ptr,row);
                    if(Y_N)
                    {
                        node* tmp = ary[col].ptr->ptr;
                        ary[col].ptr->ptr = new node(row,tmp);
                        ary[col].ptr->set_size++;
                        ary[row].ptr = ary[col].ptr;
                        set_count--;
                    }
                }
                else if(ary[col].ptr == NULL && ary[row].ptr != NULL)
                {
                    bool Y_N = check_set(ary,ary[row].ptr,col);
                    if(Y_N)
                    {
						cout<<"=========================================="<<endl;getchar();//impossible to occur?
                        node* tmp = ary[row].ptr->ptr;
                        ary[row].ptr->ptr = new node(col,tmp);
                        ary[row].ptr->set_size++;
                        ary[col].ptr = ary[row].ptr;
                        set_count--;
                    }
                }
				//else if(ary[col].ptr != NULL && ary[row].ptr != NULL) {}//each sensor already belongs to a cluster, no matter the clusters are different or not		
            }
            row++;
        }
        col++;
    }
    return set_count;
}
int c(int m,int n,vector<vector<int> >& cache)//'>' 與 '>' 要分開
{
    if (cache[m-n][n]==0)
        if (m==n||n==0) cache[m-n][n]=1;
        else cache[m-n][n]=c(m-1,n,cache)+c(m-1,n-1,cache);
    return cache[m-n][n];
}
int combination_number(int m, int n)
{
    vector<vector<int> > cache(m-n+1,vector<int>(n+1,0));//'>' 與 '>' 要分開
    return c(m,n,cache);
}
void combination(int m, int n, int* result, const int len, int** set, int& offset)
{
    if(n==0) {
		for(int i=0;i<len;i++)
		{
		    set[offset][i] = result[i];
		    //cout<<result[i]<<" ";
        }
		//cout<<endl;
		offset++;
		return;
	}
	for(int i=m; i>=n; --i) {
		result[n-1] = i-1;
		combination(i-1,n-1,result,len,set,offset);
	}
}
int extend_set(X_Y* ext_set, const int set_count, X_Y* sets_vector, coordinate* ary, int Radius, float** vector)
{
    int ext_count = 0;
    for(int i=0;i<set_count;i++)
    {
        for(int j=i+1;j<set_count;j++)
        {
            if(sets_vector[i].size != 1 || sets_vector[j].size != 1)//兩set皆單獨，在之前已檢查不會構成set
            {
                for(int p=sets_vector[i].size;p>=1;p--)
                {
                    int a = combination_number(sets_vector[i].size,p);
                    //int set1[a][p];//set1的index
                    if(a <= 0) {cout<<"Combination number is error!"<<endl;return ext_count;}
                    int** set1 = new int* [a];
                    for(int r=0;r<a;r++) set1[r] = new int [p];
                    int cc=0;
                    int result1[p];
                    combination(sets_vector[i].size,p,result1,p,set1,cc);
                    for(int q=sets_vector[j].size;q>=1;q--)
                    {
                        int b = combination_number(sets_vector[j].size,q);
                        //int set2[b][q];//set2的index
                        if(b <= 0) {cout<<"Combination number is error!"<<endl;return ext_count;}
                        int** set2 = new int* [b];
                        for(int r=0;r<b;r++) set2[r] = new int [q];
                        int dd=0;
                        int result2[q];
                        combination(sets_vector[j].size,q,result2,q,set2,dd);
                        for(int m=0;m<a;m++) for(int n=0;n<b;n++)
                        {
                            if(check_ext_set(set1[m],p,set2[n],q,i,j,sets_vector,ext_set,ext_count,Radius,vector,ary))
                            {
                                ext_set[ext_count].size = p+q;
                                ext_set[ext_count].ary = new int [ext_set[ext_count].size];
                                int cnt = 0;
                                for(int t=0;t<p;t++) ext_set[ext_count].ary[cnt++]
                                    = sets_vector[i].ary[set1[m][t]];
                                for(int t=0;t<q;t++) ext_set[ext_count].ary[cnt++]
                                    = sets_vector[j].ary[set2[n][t]];
                                ext_count++;
                                if(ext_count >= max_ext_set)
                                {
                                    cout<<"Achieve the max allowable extended sets!"<<endl;
                                    return ext_count;
                                }
                            }
                        }
                        for(int r=0;r<b;r++) delete [] set2[r];
                        delete [] set2;
                    }
                    for(int r=0;r<a;r++) delete [] set1[r];
                    delete [] set1;
                }
            }
        }
    }
    return ext_count;
}

void calculate_set_dist(X_Y* sets_vector, int set_count, float** dist_vector, int RandomNode, coordinate* ary)
{
    int set_number = 0;//目前累積set數
    for(int i=0;i<RandomNode;i++)
    { 
        if(ary[i].ptr == NULL)
        {
            sets_vector[set_number].x_axis = ary[i].x_axis;
            sets_vector[set_number].y_axis = ary[i].y_axis;
            sets_vector[set_number].size = 1;
            sets_vector[set_number].ary = new int [1];
            sets_vector[set_number].ary[0] = i;
            set_number++;
        }
        else
        {
            bool clear = true;
            for(int k=0;k<set_number;k++)
            {
                if((ary[i].ptr->set_x == sets_vector[k].x_axis) &&
                      (ary[i].ptr->set_y == sets_vector[k].y_axis)) {clear = false;break;}
            }
            if(clear)
            {
                sets_vector[set_number].x_axis = ary[i].ptr->set_x;
                sets_vector[set_number].y_axis = ary[i].ptr->set_y;
				int sz = ary[i].ptr->set_size;
                sets_vector[set_number].size = sz;
                sets_vector[set_number].ary = new int [sz];
                node* ite = ary[i].ptr->ptr;
                for(int h=0;h<sz;h++)
                {
                    sets_vector[set_number].ary[h] = ite->sensor;
                    ite = ite->next;
                }
                set_number++;
            }
        }
    }
      //cout<<set_number<<"=="<<set_count<<endl;
    if(set_number != set_count) {cerr<<"Some errors occurred!"<<endl;exit(-1);}
    int col = 1, row;
    for(int i=0;i<set_count-1;i++)
    {
        row = 0;
        for(int j=0;j<col;j++)
        {
            dist_vector[col][row] = dist_vector[row][col] =
            Euclideandistance(sets_vector[col].x_axis,sets_vector[col].y_axis,
                              sets_vector[row].x_axis,sets_vector[row].y_axis);
            row++;
        }
        col++;
    }
}
void merge(X_Y* sets_vector, int set_count, X_Y* ext_set, int ext_count, float** dist_vector_joint, X_Y* joint_vector)
{
    int offset = 0;
    for(int w=0;w<set_count;w++)
        joint_vector[offset++] = sets_vector[w];
    for(int w=0;w<ext_count;w++)
        joint_vector[offset++] = ext_set[w];
    delete [] sets_vector;
    delete [] ext_set;

    int col = 1, row;
    for(int i=0;i<(set_count+ext_count-1);i++)
    {
        row = 0;
        for(int j=0;j<col;j++)
        {
            dist_vector_joint[col][row] = dist_vector_joint[row][col] =
            Euclideandistance(joint_vector[col].x_axis,joint_vector[col].y_axis,
                              joint_vector[row].x_axis,joint_vector[row].y_axis);
            row++;
        }
        col++;
    }
}
float greedy(int set_count, float** dist_vector, int* chromo)//return total cost
{
    int gene = 0;
    int mark[set_count];
	for(int i=0;i<set_count;i++) mark[i] = 1;

    int WhichSet, Start_Set = 0, End_Set = Start_Set;//default
    float min, Accumulated_Cost = 0;
    chromo[gene++] = Start_Set;
    mark[Start_Set] = 0;
    //cout<<Start_Set<<" ";
    for(int j=0;j<set_count-1;j++)
    {
        min = INT_MAX;
        for(int i=0;i<set_count;i++)
        {
            if(i != Start_Set && mark[i] == 1)
              if(dist_vector[i][Start_Set] < min)
              {
                  min = dist_vector[i][Start_Set];
                  WhichSet = i;
              }
        }
         //cout<<WhichSet<<" ";
        chromo[gene++] = WhichSet;
        mark[WhichSet] = 0;//picked
        Start_Set = WhichSet;
        Accumulated_Cost += min;
    }
    float total_cost = Accumulated_Cost+dist_vector[End_Set][Start_Set];
    return total_cost;
}
void permutation(int n, int* solution, bool* used, int& count, int num, int** all_permu)
{
    if (n == num)
    {
        if(count < 0) return;//reach population size
        for (int i=0; i<num; i++)
        {
            all_permu[count][i] = solution[i];
            //cout << solution[i] << " ";
        }
        //cout << endl;
        count--;
        return;
    }
    for (int i=0; i<num; i++)
        if (!used[i])
        {
            used[i] = true;
            solution[n] = i;
            permutation(n+1,solution,used,count,num,all_permu);
            used[i] = false;
        }
}
bool check_valid(chromosome& pp, int sz)
{
    int flag[sz];
	for(int i=0;i<sz;i++) flag[i] = 1;
    for(int i=0;i<sz;i++)
    {
        if(flag[pp.chromo[i]] == 0) {cout<<"invalid chromosome!"<<endl;getchar();return false;}
        else flag[pp.chromo[i]] = 0;
    }
    return true;
}
chromosome* CGA(int* first_chromo, int set_count, int& population)//Chromosome Generation Algorithm
{
    int factorial = set_count-1, num = set_count-1;
    while(factorial < population && num != 2) factorial *= --num;
    if(num == 2) population = factorial;//population is too large
//cout<<"population: "<<population<<endl;
    
	factorial = 1, num = 1;
    while(factorial < population) factorial *= ++num;//取出用作permutation的set數可產生大於population的量

    int** all_permu = new int* [population];//0~(population-1)
    for(int i=0;i<population;i++) all_permu[i] = new int [num];
//cout<<"num: "<<num<<endl;
    int pos;
    int permu_index[num];//record which positions are permutated from the first chromosome
    int permu_value[num];
    for(int i=0;i<num;i++)
    {//不可重複
    random_again:
        pos = rand()%(set_count-1) + 1;//1~(set_count-1), start point cannot be chosen
        for(int n=0;n<i;n++) if(permu_index[n] == pos) goto random_again;
        permu_index[i] = pos;
        permu_value[i] = first_chromo[permu_index[i]]; 
    }
    int solution[num];    // 用來存放一組可能的答案
    bool used[num];   // 紀錄數字是否使用過，用過為 true
    for (int i=0; i<num; i++) used[i] = false;
    int count = population-1;
    permutation(0,solution,used,count,num,all_permu);
//cout<<count<<"======";getchar();
    chromosome* all_chromosome = new chromosome [population];//0~(population-1)
	all_chromosome[population-1].chromo = new int [set_count];
	for(int j=0;j<set_count;j++) all_chromosome[population-1].chromo[j] = first_chromo[j];
    for(int i=0;i<population-1;i++)
    {
        all_chromosome[i].chromo = new int [set_count];
        for(int j=0;j<set_count;j++) all_chromosome[i].chromo[j] = first_chromo[j];//default start from set 0
        for(int h=0;h<num;h++) all_chromosome[i].chromo[permu_index[h]] = permu_value[all_permu[i][h]];
    }

    for(int i=0;i<population;i++) delete [] all_permu[i];
    delete [] all_permu;
//for(int i=0;i<population;i++) check_valid(all_chromosome[i],set_count);
    return all_chromosome;
}
chromosome* random_initial(int* first_chromo, int set_count, int population)
{
    chromosome* all_chromosome = new chromosome [population];//0~(population-1)
    for(int i=0;i<population;i++)
    {
        all_chromosome[i].chromo = new int [set_count];
        for(int j=0;j<set_count;j++) all_chromosome[i].chromo[j] = first_chromo[j];//default start from set 0
        for(int j=1;j<set_count;j++)
        {
            int index = rand()%(set_count-1) + 1;
            swap(all_chromosome[i].chromo[j],all_chromosome[i].chromo[index]);
        } 
    }
//for(int i=0;i<population;i++) check_valid(all_chromosome[i],set_count);
    return all_chromosome;
}
bool check_prime(int num)
{
	if(num < 2) {cout<<"Input number is too small!!"<<endl;return false;}
	for(int i=2;i<=sqrt(num);i++) if(num%i == 0) return false;
	return true;
}
chromosome* balance_SD(int set_count, int& population)
{
	int prime_num = set_count;
	while(!check_prime(prime_num)) prime_num--;
	//cout<<" prime: "<<prime_num<<"";
	//population = (population/(prime_num-1)+1)*(prime_num-1);
	//cout<<population;getchar();
	 
	int add = 0;//1~(prime_num-1)
    chromosome* all_chromosome = new chromosome [population];//0~(population-1)
    for(int i=0;i<population;i++)
    {
        //add = rand()%(prime_num-1)+1;//-------1
		add = add%(prime_num-1)+1;//-------2
		all_chromosome[i].chromo = new int [set_count];
		all_chromosome[i].chromo[0] = 0;
        for(int j=1;j<prime_num;j++) all_chromosome[i].chromo[j] = (all_chromosome[i].chromo[j-1] + add)%prime_num;
		for(int m=prime_num;m<set_count;m++) all_chromosome[i].chromo[m] = m; 
		for(int n=prime_num;n<set_count;n++)
        {
			int index = rand()%(set_count-1) + 1;
            swap(all_chromosome[i].chromo[n],all_chromosome[i].chromo[index]);
        }
    }
	//============================
	/*for(int i=0;i<population;i++) 
	{
		check_valid(all_chromosome[i],set_count);
		//for(int h=0;h<set_count;h++) cout<<all_chromosome[i].chromo[h]<<" ";cout<<endl;getchar();
	}*/
    return all_chromosome;
}
chromosome* balance_SD_2(int set_count, int& population)
{
	int prime_num = set_count;
	while(!check_prime(prime_num)) prime_num--;
	//cout<<" prime: "<<prime_num<<"";
	//population = (population/(prime_num-1)+1)*(prime_num-1);
	//cout<<population;getchar();
	 
	int add = 0;//1~(prime_num-1)
    chromosome* all_chromosome = new chromosome [population];//0~(population-1)
    for(int i=0;i<population;i++)
    {
        //add = rand()%(prime_num-1)+1;//-------1
		add = add%(prime_num-1)+1;//-------2
		all_chromosome[i].chromo = new int [set_count];
		for(int q=0;q<set_count;q++) all_chromosome[i].chromo[q] = 0;
        //for(int j=1;j<prime_num;j++) all_chromosome[i].chromo[j] = (all_chromosome[i].chromo[j-1] + add)%prime_num;
		//for(int m=prime_num;m<set_count;m++) all_chromosome[i].chromo[m] = m; 
		for(int n=prime_num;n<set_count;n++)
        {
			int index;
			do{index = rand()%(set_count-1) + 1;}while(all_chromosome[i].chromo[index]!=0);
			all_chromosome[i].chromo[index] = n;
            //swap(all_chromosome[i].chromo[n],all_chromosome[i].chromo[index]);
        }
		int record = 0;
		for(int j=1;j<set_count;j++) 
		{	
			if(all_chromosome[i].chromo[j]!=0) continue;
			else record = all_chromosome[i].chromo[j] = (record + add)%prime_num;
		}
    }
	//============================
	/*for(int i=0;i<population;i++) 
	{
		check_valid(all_chromosome[i],set_count);
		//for(int h=0;h<set_count;h++) cout<<all_chromosome[i].chromo[h]<<" ";cout<<endl;getchar();
	}*/
    return all_chromosome;
}
void fitness(int population, chromosome* all_chromosome, int set_count, float** dist_vector)
{//calculate fitness
    for(int i=0;i<population;i++)
    {
        int start = 0, end = start;
        float cost = 0;
        for(int j=1;j<set_count;j++)
        {
            cost += dist_vector[start][all_chromosome[i].chromo[j]];
            start = all_chromosome[i].chromo[j];
        }
        cost += dist_vector[start][end];
        all_chromosome[i].fitness_val = cost;
          //cout<<"=="<<cost<<" "<<endl;
    }
}
int circle_int(float x1, float y1, float x2, float y2, int r1, int r2, X_Y* list)//回傳幾個交點數
{
    int cnt = 0;
    if(y1 != y2)
    {
        float m = (x1-x2)/(y2-y1);
        float k = (r1*r1 - r2*r2 + x2*x2 - x1*x1 + y2*y2 - y1*y1)/(2*(y2-y1));
        float a = 1 + m*m;
        float b = 2*(m*k - m*y2 - x2);
        float c = x2*x2 + y2*y2 + k*k - 2*k*y2 - r2*r2;

        float res = b*b - 4*a*c;
        if(res > 0)
        {
            list[cnt].x_axis = (-1*b + sqrt(res))/(2*a);
            list[cnt].y_axis = m*(list[cnt].x_axis) + k;
            cnt++;

            list[cnt].x_axis = (-1*b - sqrt(res))/(2*a);
            list[cnt].y_axis = m*(list[cnt].x_axis) + k;
            cnt++;
        }
        else if(res == 0)//相切
        {
            list[cnt].x_axis = -1*b/(2*a);
            list[cnt].y_axis = m*(list[cnt].x_axis) + k;
            cnt++;
        }
    }
    else//y1 == y2
    {
        if(x1 != x2)
        {
            float x = (r1*r1 - r2*r2 + x2*x2 - x1*x1)/(2*(x2-x1));
		    float a = 1;
            float b = -2*y1;
            float c = x*x + x1*x1 - 2*x1*x + y1*y1 - r1*r1;

            float res = b*b - 4*a*c;
            if(res > 0)
            {
                list[cnt].x_axis = x;
                list[cnt].y_axis = (-1*b + sqrt(res))/(2*a);
                cnt++;

                list[cnt].x_axis = x;
                list[cnt].y_axis = (-1*b - sqrt(res))/(2*a);
                cnt++;
            }
            else if(res == 0)//相切
            {
                list[cnt].x_axis = x;
                list[cnt].y_axis = y1;//== y2 == (-1*b)/(2*a)
                cnt++;
            }
        }
        else ;//x1 == x2, y1 == y2
    }
    return cnt;
}
bool check_node_vaild(float x_, float y_, X_Y& now_set, coordinate* ary)
{
    for(int i=0;i<now_set.size;i++)
    {
        float dist = Euclideandistance(ary[now_set.ary[i]].x_axis,ary[now_set.ary[i]].y_axis,x_,y_);
		//if(dist > ary[now_set.ary[i]].r ) {cout<<"=="<<(dist - ary[now_set.ary[i]].r)<<"==";getchar();return false;}
		if(dist > ary[now_set.ary[i]].r && (dist - ary[now_set.ary[i]].r) > 0.05) return false;//浮點數有誤差值
    }
    return true;
}
void calculate_set_node(node_list* cache, coordinate* ary, X_Y* sets_vector, int set_count)
{
    for(int m=0;m<set_count;m++)
    {
        int SZ = sets_vector[m].size;
        int index = 0;
        cache[m].set_node_x = new float [SZ*(SZ-1)+1];//include original waypoint
        cache[m].set_node_y = new float [SZ*(SZ-1)+1];
        cache[m].set_node_x[index] = sets_vector[m].x_axis;
        cache[m].set_node_y[index] = sets_vector[m].y_axis;
        index++;
        cache[m].set_node_sz = 1;

        int col = 1, row;
        for(int i=0;i<SZ-1;i++)
        {
            row = 0;
            for(int j=0;j<col;j++)
            {
                //find intersection
                int sensor1 = sets_vector[m].ary[col];
                int sensor2 = sets_vector[m].ary[row];
                X_Y list[2];
                int int_number = circle_int(ary[sensor1].x_axis,ary[sensor1].y_axis,
                           ary[sensor2].x_axis,ary[sensor2].y_axis,ary[sensor1].r,ary[sensor2].r,list);//至少一交點 or 完全相同(0)
                //chech valid
                for(int k=0;k<int_number;k++)//一次或兩次
                {
                    if(check_node_vaild(list[k].x_axis,list[k].y_axis,sets_vector[m],ary))//檢查某點是否在共同交集區
                    {
                        //store into list
                        cache[m].set_node_x[index] = list[k].x_axis;
                        cache[m].set_node_y[index] = list[k].y_axis;
                        index++;
                        cache[m].set_node_sz++;
                    }
                }
                row++;
            }
            col++;
        }
    }
}
int find_cross_node(float x_, float y_, X_Y& now_set, X_Y* list, coordinate* ary, int strat_from=0)
{
    int cnt = 0;
    for(int i=strat_from;i<now_set.size;i++)
    {
        float pp = ary[now_set.ary[i]].x_axis;
        float qq = ary[now_set.ary[i]].y_axis;

        float dist = Euclideandistance(x_,y_,pp,qq);
		float rr = ary[now_set.ary[i]].r;
        if(dist >= rr)//有交點
        {
            if(x_ != pp)
            {
                float m = (y_ - qq)/(x_ - pp);
                float n = y_ - m*x_;
                float a = 1 + m*m;
                float b = -2*pp + 2*m*n - 2*qq*m;
                float c = pp*pp + qq*qq - rr*rr - 2*qq*n + n*n;

                float res = b*b - 4*a*c;
                if(res > 0)
                {
                    list[cnt].x_axis = (-1*b + sqrt(res))/(2*a);
                    list[cnt].y_axis = m*(list[cnt].x_axis) + n;
                    if(Euclideandistance(list[cnt].x_axis,list[cnt].y_axis,pp,qq) <= dist
                       &&
                       Euclideandistance(list[cnt].x_axis,list[cnt].y_axis,x_,y_) <= dist) cnt++;
                    else
                    {
                        list[cnt].x_axis = (-1*b - sqrt(res))/(2*a);
                        list[cnt].y_axis = m*(list[cnt].x_axis) + n;
                        cnt++;
                    }
                }
            }
            else
            {
                list[cnt].x_axis = pp;
                if(y_ > qq) list[cnt].y_axis = qq + rr;
                else list[cnt].y_axis = qq - rr;
                cnt++;
            }
        }
    }
    return cnt;
}
float degree_case(float xx, float yy)
{
	//"atan" Range: (-pi/2, pi/2), Domain: all real numbers
	if((xx>0 && yy>0) || (xx>0 && yy<0)) return DEGREEOF(atan(yy/xx));//quadrant 1, 4
	else if((xx<0 && yy>0) || (xx<0 && yy<0)) return DEGREEOF(atan(yy/xx))+180.0;//quadrant 2, 3
	else if(xx==0 && yy>0) return 90.0;//upper line
	else if(xx==0 && yy<0) return 270.0;//lower line
	else if(xx>0 && yy==0) return 0.0;//right line
	else if(xx<0 && yy==0) return 180.0;//left line
	else if(xx==0 && yy==0) return 0.0;//center, impossible here
}
float calculate_rotate(float a_x,float a_y,float b_x,float b_y,float c_x,float c_y)
{
	float len_a, len_b, len_c;
	len_a = Euclideandistance(b_x,b_y,c_x,c_y);//BC
	len_b = Euclideandistance(a_x,a_y,c_x,c_y);//AC
	len_c = Euclideandistance(a_x,a_y,b_x,b_y);//AB
	if(len_a == 0 || len_c == 0) return 0;//some points are repeat
	
	/*//by arccos
	//"acos" Range: [0, pi], Domain: [-1, 1]
	float cos_B = ((len_a*len_a)+(len_c*len_c)-(len_b*len_b))/(2*len_a*len_c);//cos(B)
	if(cos_B>1 || cos_B<-1) {cout<<" => "<<(double)cos_B<<" = "<<len_a<<" "<<len_b<<" "<<len_c<<endl;getchar();}
	float radian = acos(cos_B);
	return 180.0 - DEGREEOF(radian);//degree*/
	
	//by arctan
	float a = (a_y-b_y);
	float b = (a_x-b_x);
	float c = (c_y-b_y);
	float d = (c_x-b_x);
	float degree1 = degree_case(a,b);//degree
	float degree2 = degree_case(c,d);//degree
	float diff = max(degree1,degree2)-min(degree1,degree2);
	if(diff > 180) diff = 360.0 - diff;
    return 180.0-diff;//degree
}
bool check_segment_circle(float x_, float y_, float x__, float y__, float pp, float qq, float rr)
{//檢查線段是否和圓有交點
    //(x_,y_), (x__,y__)表示線段兩端點
	//(pp,qq)表示圓心
	float m = (y_ - y__)/(x_ - x__);//斜率
    float n = y_ - m*x_;
    float a = 1 + m*m;
    float b = -2*pp + 2*m*n - 2*qq*m;
    float c = pp*pp + qq*qq - rr*rr - 2*qq*n + n*n;
	
    float res = b*b - 4*a*c;
	
	float XX, YY;
	bool cnt = false;
	
	if(res == 0)//相切
	{
		XX = (-1*b)/(2*a);
        YY = m*XX + n;
		if(!((XX > max(x_,x__) || XX < min(x_,x__)) || (YY > max(y_,y__) || YY < min(y_,y__))))//檢查點是否在線段上
		{
			cnt = true;
		}
	}
	else if(res > 0)//直線割圓
	{
		XX = (-1*b + sqrt(res))/(2*a);
        YY = m*XX + n;
		if(!((XX > max(x_,x__) || XX < min(x_,x__)) || (YY > max(y_,y__) || YY < min(y_,y__))))//檢查點是否在線段上
		{
			cnt = true;
		}
		//==================================
		XX = (-1*b - sqrt(res))/(2*a);
        YY = m*XX + n;
		if(!((XX > max(x_,x__) || XX < min(x_,x__)) || (YY > max(y_,y__) || YY < min(y_,y__))))//檢查點是否在線段上
		{
			cnt = true;
		}
		//==================================
		if( ((pp-x_)*(pp-x_) + (qq-y_)*(qq-y_) - rr*rr) < 0 && ((pp-x__)*(pp-x__) + (qq-y__)*(qq-y__) - rr*rr) < 0)//線段在圓內
		{//線段包含在圓內
			cnt = true;
		}
	}
    return cnt;
}
void fitness_CP(int population, chromosome* all_chromosome, int set_count, X_Y* sets_vector,
                coordinate* ary, node_list* cache, X_Y* path=NULL)//fitness_val表示chromosome走訪長度, 越小越好
{//calculate fitness
	for(int i=0;i<population;i++)
    {
        float start_x = sets_vector[0].x_axis;
        float start_y = sets_vector[0].y_axis;
        float end_x = start_x;
        float end_y = start_y;
        float cost = 0;
		float max_step = INT_MIN;
		float rotate = 0, a_x, a_y, b_x, b_y;
		a_x = b_x = start_x;
		a_y = b_y = start_y;
		int index = 0;
		
if(path != NULL) print_waypoint<<"Print out visited points: (No LLA)"<<endl;
        for(int jj=1;jj<set_count;jj++)
        {
if(path != NULL) print_waypoint<<" ("<<start_x<<", "<<start_y<<")";
if(path != NULL)
{
	path[index].x_axis = start_x;
	path[index].y_axis = start_y;
	index++;
}
            //find線段和圓的交點
            int j = all_chromosome[i].chromo[jj];//next city
            X_Y cross_list[sets_vector[j].size];//所有可能性
            int int_number = find_cross_node(start_x,start_y,sets_vector[j],cross_list,ary);//可能無cross
            //check valid, choose best in cross_list
            float Min = INT_MAX;
            float XX, YY;
            for(int k=0;k<int_number;k++)
            {
                if(check_node_vaild(cross_list[k].x_axis,cross_list[k].y_axis,sets_vector[j],ary))//檢查某點是否在共同交集區
                {
                    float dist = Euclideandistance(cross_list[k].x_axis,cross_list[k].y_axis,start_x,start_y);
                    if(dist < Min)
                    {
                        XX = cross_list[k].x_axis;
                        YY = cross_list[k].y_axis;
                        Min = dist;
                    }
                }
            }//cout<<Min<<"--";getchar();
            //compare with set_node
            for(int s=0;s<cache[j].set_node_sz;s++)//nodes in cache are already checked by "check_node_vaild"
            {
                float dist = Euclideandistance(cache[j].set_node_x[s],cache[j].set_node_y[s],start_x,start_y);
                if(dist < Min)
                {
                    XX = cache[j].set_node_x[s];
                    YY = cache[j].set_node_y[s];
                    Min = dist;
                }
            }//cout<<Min<<"==";getchar();
			rotate += calculate_rotate(a_x,a_y,b_x,b_y,XX,YY);
			a_x = b_x;
			a_y = b_y;
			b_x = XX;
			b_y = YY;
			//update
            start_x = XX;
            start_y = YY;
            cost += Min;
			if(Min > max_step) max_step = Min; 
        } 
if(path != NULL) print_waypoint<<" ("<<start_x<<", "<<start_y<<")"<<" ("<<end_x<<", "<<end_y<<")"<<endl;
if(path != NULL)
{
	path[index].x_axis = start_x;
	path[index].y_axis = start_y;
	index++;
}
		
        float last_step = Euclideandistance(start_x,start_y,end_x,end_y);
		cost += last_step;
		rotate += calculate_rotate(a_x,a_y,b_x,b_y,end_x,end_y);
		if(last_step > max_step) max_step = last_step;
        all_chromosome[i].fitness_val = cost;
		all_chromosome[i].rotate_degree = rotate;
		all_chromosome[i].max_step_len = max_step;
		  //cout<<"=="<<max_step<<" "<<endl;
		  //cout<<"=="<<rotate<<" "<<endl;
          //cout<<"=="<<cost<<" "<<endl; 
    }
}
int final_fitness_CP(int best_chro, chromosome* all_chromosome, int set_count, X_Y* sets_vector,
                coordinate* ary, X_Y* path)//return # of shortcut success 	
{
	float start_x = path[0].x_axis;
    float start_y = path[0].y_axis;
    float end_x = start_x;
    float end_y = start_y;
    float cost = 0;
	float previous_step_length = 1; 
	bool fg=true;
	float max_step = INT_MIN;
	float rotate = 0;
	float a_x, a_y, b_x, b_y;
	a_x = b_x = start_x;
	a_y = b_y = start_y;	
	int success_count = 0;	
	int error_cnt = 2;//0:Advanced-LLA, 2:LLA 
	
//footprint<<"Print out segments: "<<endl;			
if(path != NULL) print_waypoint<<"Print out visited points: (LLA)"<<endl;
if(path != NULL) print_waypoint<<" ("<<start_x<<", "<<start_y<<")";	
	int jj;
	for(jj=1;jj<set_count;jj++)
	{
		int mid_set = all_chromosome[best_chro].chromo[jj];
		int latter_index = (jj+1)%set_count;
		int end_set = all_chromosome[best_chro].chromo[latter_index];
		
		float x2 = path[latter_index].x_axis;
		float y2 = path[latter_index].y_axis;
		float base = Euclideandistance(start_x,start_y,x2,y2);
		
		int num_sensor = sets_vector[mid_set].size;

		float x_detour, y_detour;//只"嘗試"多繞一次路
		for(int j=0;j<num_sensor;j++)
		{
			int mid_sensor_id = sets_vector[mid_set].ary[j];
			float x3 = ary[mid_sensor_id].x_axis;//(x3,y3)為中間圓心
			float y3 = ary[mid_sensor_id].y_axis;
			float radius = ary[mid_sensor_id].r;

			float area = 0.5*(abs((start_x*y2+x2*y3+x3*start_y)-(x2*start_y+x3*y2+start_x*y3)));//行列式三點求面積
			float hight = 2*area/base;
			if(!check_segment_circle(start_x,start_y,x2,y2,x3,y3,radius))
			{
				if(error_cnt == 0)//---------------
				{
					/*float seg1, seg2;
					float tmp;
					tmp = Euclideandistance(x1,y1,x3,y3);
					seg1 = sqrt((tmp*tmp)-(hight*hight));
					
					tmp = Euclideandistance(x2,y2,x3,y3);
					seg2 = sqrt((tmp*tmp)-(hight*hight));
					
					x_detour = (seg1*x2+seg2*x1)/(seg1+seg2);//hight和base的垂直點
					y_detour = (seg1*y2+seg2*y1)/(seg1+seg2);	
//cout<<"=="<<hight<<"=="<<Euclideandistance(x3,y3,x_detour,y_detour)<<"==";getchar();
					X_Y cross_list[sets_vector[mid_set].size];//所有可能性
					int int_number = find_cross_node(x_detour,y_detour,sets_vector[mid_set],cross_list,ary,j);//可能無cross					
					x_detour = cross_list[0].x_axis;
					y_detour = cross_list[0].y_axis;
//cout<<"=="<<ary[mid_sensor_id].r<<"=="<<Euclideandistance(x3,y3,x_detour,y_detour)<<"==";getchar();	
//cout<<"=="<<ary[mid_sensor_id].r<<"=="<<hight<<"==";getchar();
//cout<<"=="<<(hight-ary[mid_sensor_id].r)<<"==";getchar();
					if(hight-ary[mid_sensor_id].r >= 30) {fg = false;break;}	//需要繞太遠		
					else error_cnt++;*/			
				}
				else
				{
					fg = false;
					break;
				}
			}
		}
		//cout<<"=="<<error_cnt<<"==";getchar();

		float Min;
		if(fg)//shortcut is successful
		{
			if(error_cnt == 1)//只接受多繞一次路//-------------------
			{
				float detour_length = Euclideandistance(start_x,start_y,x_detour,y_detour) + Euclideandistance(x_detour,y_detour,x2,y2);
				Min = detour_length;
				//cout<<"=="<<detour_length<<"=="<<base<<"==";getchar();
			}
			else Min = base;
			jj++;
			start_x = x2;
			start_y = y2;
			success_count++;
		}
		else
		{
			float new_start_x = path[jj].x_axis;
			float new_start_y = path[jj].y_axis;
			Min = Euclideandistance(start_x,start_y,new_start_x,new_start_y);
			start_x = new_start_x;
			start_y = new_start_y;
		}
if(path != NULL) print_waypoint<<" ("<<start_x<<", "<<start_y<<")";
		cost += Min;
		if(Min > max_step) max_step = Min; 
		
		float this_time_rotate = calculate_rotate(a_x,a_y,b_x,b_y,start_x,start_y);
		rotate += this_time_rotate;
		//cout<<this_time_rotate/previous_step_length<<" ";
		//cout<<this_time_rotate<<" ";	
		//previous_step_length = Min;
//footprint<<previous_step_length<<" ";
		a_x = b_x;
		a_y = b_y;
		b_x = start_x;
		b_y = start_y;
	}

	if(jj==set_count)//補加最後一段距離
	{
		float last_step = Euclideandistance(start_x,start_y,end_x,end_y);
		cost += last_step;
		float this_time_rotate = calculate_rotate(a_x,a_y,b_x,b_y,end_x,end_y);
		//cout<<this_time_rotate/previous_step_length<<" ";
		//cout<<this_time_rotate<<" ";	
//footprint<<last_step<<" ";
if(path != NULL) print_waypoint<<" ("<<end_x<<", "<<end_y<<")";
		rotate += this_time_rotate;
		if(last_step > max_step) max_step = last_step;
	}
if(path != NULL) print_waypoint<<endl;
	all_chromosome[best_chro].fitness_val = cost;
	all_chromosome[best_chro].rotate_degree = rotate;
	all_chromosome[best_chro].max_step_len = max_step;
	return success_count;
}
void IOpt(chromosome& chro ,int set_count, float** dist_vector,
          X_Y* sets_vector, coordinate* ary, node_list* cache)
{
    if(set_count < 4) {cerr<<"chromosome is too small!"<<endl;return;}
    int point1, point2;

        //find a NMTS to be exchanged
        float max_diff = 0;
        int NMTS;
        for(int i=0;i<set_count-1;i++)
        {
            int now_node = chro.chromo[i];
            float min = INT_MAX;

            for(int j=1;j<set_count;j++)
              if(j!=now_node && dist_vector[now_node][j]<min) min = dist_vector[now_node][j];

            float dd = dist_vector[now_node][chro.chromo[i+1]] - min;
            if(dd > max_diff)
            {
                max_diff = dd;
                NMTS = i;
            }
        }
        if(max_diff == 0) return;//no NMTS

        //find another edge
        float best = INT_MAX;
        int match;
        const float original = chro.fitness_val;//original
        int* adjust = new int [set_count];

        point2 = NMTS+1;
        for(int i=0;i<=NMTS-2;i++)
        {
            point1 = i+1;

            /*int _p1 = chro.chromo[point1-1];
            int p1  = chro.chromo[point1];
            int _p2 = chro.chromo[point2-1];
            int p2  = chro.chromo[point2];
            float after = dist_vector[_p1][_p2] + dist_vector[p1][p2];//new
            float before = dist_vector[p1][_p1] + dist_vector[p2][_p2];//old
            if(after<before && (after-before)<best)//越小表減少越多cost
            {
                best = after - before;
                match = point1;
            }*/
                int pos = 0;
                for(int cnt=0;cnt<=point1-1;cnt++) adjust[pos++] = chro.chromo[cnt];
                for(int cnt=point2-1;cnt>=point1;cnt--) adjust[pos++] = chro.chromo[cnt];
                for(int cnt=point2;cnt<set_count;cnt++) adjust[pos++] = chro.chromo[cnt];
                swap(chro.chromo,adjust);
                fitness_CP(1,&chro,set_count,sets_vector,ary,cache);
                if(chro.fitness_val<original && (chro.fitness_val-original)<best)//find one candidate
                {
                    best = chro.fitness_val - original;
                    match = point1;
                }
                swap(chro.chromo,adjust);//restore chro
                chro.fitness_val = original;
        }
        //=====
        point1 = NMTS+1;
        for(int i=NMTS+2;i<=set_count-2;i++)
        {
            point2 = i+1;

            /*int _p1 = chro.chromo[point1-1];
            int p1  = chro.chromo[point1];
            int _p2 = chro.chromo[point2-1];
            int p2  = chro.chromo[point2];
            float after = dist_vector[_p1][_p2] + dist_vector[p1][p2];//new
            float before = dist_vector[p1][_p1] + dist_vector[p2][_p2];//old
            if(after<before && (after-before)<best)//越小表減少越多cost
            {
                best = after - before;
                match = point2;
            }*/
                int pos = 0;
                for(int cnt=0;cnt<=point1-1;cnt++) adjust[pos++] = chro.chromo[cnt];
                for(int cnt=point2-1;cnt>=point1;cnt--) adjust[pos++] = chro.chromo[cnt];
                for(int cnt=point2;cnt<set_count;cnt++) adjust[pos++] = chro.chromo[cnt];
                swap(chro.chromo,adjust);
                fitness_CP(1,&chro,set_count,sets_vector,ary,cache);
                if(chro.fitness_val<original && (chro.fitness_val-original)<best)//find one candidate
                {
                    best = chro.fitness_val - original;
                    match = point2;
                }
                swap(chro.chromo,adjust);//restore chro
                chro.fitness_val = original;
        }
        //=====
        if(best != INT_MAX)//exist an exchangeable edge
        {
            point1 = NMTS+1;
            point2 = match;
            if(point1 > point2) swap(point1,point2);

            int pos = 0;
            for(int cnt=0;cnt<=point1-1;cnt++) adjust[pos++] = chro.chromo[cnt];
            for(int cnt=point2-1;cnt>=point1;cnt--) adjust[pos++] = chro.chromo[cnt];
            for(int cnt=point2;cnt<set_count;cnt++) adjust[pos++] = chro.chromo[cnt];
            swap(chro.chromo,adjust);
        }
        delete [] adjust;
}
void _2_Opt(chromosome& chro ,int set_count, float** dist_vector,
            X_Y* sets_vector, coordinate* ary, node_list* cache)
{//do until to find a better chromosome
    if(set_count < 4) {cerr<<"chromosome is too small!"<<endl;return;}
    int point1, point2;
    int try_times = 10;//chances to find a better chromosome

    while(try_times > 0)
    {
        //find two valid edges
        point1 = rand()%(set_count-1) + 1;
        do{point2 = rand()%(set_count-1) + 1;}while(point1 == point2);
        if(point1 > point2 && point1-point2 <= 1) continue;
        else if(point1 < point2 && point2-point1 <= 1) continue;
        else if(point1 > point2) swap(point1,point2);//set point1 the first edge

        /*int _p1 = chro.chromo[point1-1];
        int p1  = chro.chromo[point1];
        int _p2 = chro.chromo[point2-1];
        int p2  = chro.chromo[point2];

        if(dist_vector[_p1][_p2] + dist_vector[p1][p2]
           <
           dist_vector[p1][_p1] + dist_vector[p2][_p2])
        {
            int* adjust = new int [set_count];
            int pos = 0;
            for(int cnt=0;cnt<=point1-1;cnt++) adjust[pos++] = chro.chromo[cnt];
            for(int cnt=point2-1;cnt>=point1;cnt--) adjust[pos++] = chro.chromo[cnt];
            for(int cnt=point2;cnt<set_count;cnt++) adjust[pos++] = chro.chromo[cnt];
            swap(chro.chromo,adjust);
            delete [] adjust;
            break;//found one
        }
		else try_times--;*/
            float tmp = chro.fitness_val;//original
            int* adjust = new int [set_count];
            int pos = 0;
            for(int cnt=0;cnt<=point1-1;cnt++) adjust[pos++] = chro.chromo[cnt];
            for(int cnt=point2-1;cnt>=point1;cnt--) adjust[pos++] = chro.chromo[cnt];
            for(int cnt=point2;cnt<set_count;cnt++) adjust[pos++] = chro.chromo[cnt];
            swap(chro.chromo,adjust);
            fitness_CP(1,&chro,set_count,sets_vector,ary,cache);
            if(chro.fitness_val >= tmp)//fail
            {
                swap(chro.chromo,adjust);
                chro.fitness_val = tmp;
                try_times--;
                delete [] adjust;
            }
            else
            {
                delete [] adjust;
                break;
            }
    }
}
void EM(chromosome& chro ,int set_count, float** dist_vector,
        X_Y* sets_vector, coordinate* ary, node_list* cache)
{
    if(set_count < 4) {cerr<<"chromosome is too small!"<<endl;return;}
    int point1, point2;
    int try_times = 10;//chances to find a better chromosome

    while(try_times > 0)
    {
        point1 = rand()%(set_count-1) + 1;
        do{point2 = rand()%(set_count-1) + 1;}while(point1 == point2);

        float tmp = chro.fitness_val;
        swap((chro.chromo[point1]),(chro.chromo[point2]));
        fitness_CP(1,&chro,set_count,sets_vector,ary,cache);
        if(chro.fitness_val >= tmp)
        {
            swap((chro.chromo[point1]),(chro.chromo[point2]));//fail
            chro.fitness_val = tmp;
            try_times--;
        }
        else break;//find one
    }
}
void _Best_2_Opt(chromosome& chro ,int set_count, float** dist_vector,
                 X_Y* sets_vector, coordinate* ary, node_list* cache)
{
    if(set_count < 4) {cerr<<"chromosome is too small!"<<endl;return;}
    const int repeat_times = 1;//can be modified
    int point1 = rand()%(set_count-3) + 1;//choose a start point to mutate

    for(int i=0;i<repeat_times;i++)
    {
        for(int point2=point1+1; point2<=set_count-2; point2++)
        {
            /*int _p1 = chro.chromo[point1-1];
            int p1  = chro.chromo[point1];
            int p1_ = chro.chromo[point1+1];
            int _p2 = chro.chromo[point2-1];
            int p2  = chro.chromo[point2];
            int p2_ = chro.chromo[point2+1];

            if((dist_vector[_p1][p2] + dist_vector[p1_][p2]
               + dist_vector[p1][_p2] + dist_vector[p1][p2_])
                    <
               (dist_vector[_p1][p1] + dist_vector[p1][p1_]
               + dist_vector[_p2][p2] + dist_vector[p2][p2_]))
            swap((chro.chromo[point1]),(chro.chromo[point2]));*/

                float tmp = chro.fitness_val;
                swap((chro.chromo[point1]),(chro.chromo[point2]));
                fitness_CP(1,&chro,set_count,sets_vector,ary,cache);
                if(chro.fitness_val >= tmp)
                {
                    swap((chro.chromo[point1]),(chro.chromo[point2]));//fail
                    chro.fitness_val = tmp;
                }
        }
        point1++;
        if(point1 > set_count-3) break;
    }
}
int min_value_chromosome(int population, chromosome* ary, const char* ptr)
{
    int which_chro;
    float metric_value = INT_MAX;
	
	if(strcmp(ptr,"fitness_val")==0)
	{
		for(int h=0;h<population;h++)
		{
			if(ary[h].fitness_val < metric_value)
			{
				which_chro = h;
				metric_value = ary[h].fitness_val;
			}
		}
	}
	else if(strcmp(ptr,"rotate_degree")==0)
	{
		for(int h=0;h<population;h++)
		{
			if(ary[h].rotate_degree < metric_value)
			{
				which_chro = h;
				metric_value = ary[h].rotate_degree;
			}
		}
	}
	else if(strcmp(ptr,"max_step_len")==0)
	{
		for(int h=0;h<population;h++)
		{
			if(ary[h].max_step_len < metric_value)
			{
				which_chro = h;
				metric_value = ary[h].max_step_len;
			}
		}
	}
    return which_chro;
}
int max_value_chromosome(int population, chromosome* ary, const char* ptr)
{
    int which_chro;
    float metric_value = INT_MIN;
	
	if(strcmp(ptr,"fitness_val")==0)
	{
		for(int h=0;h<population;h++)
		{
			if(ary[h].fitness_val > metric_value)
			{
				which_chro = h;
				metric_value = ary[h].fitness_val;
			}
		}
	}
	else if(strcmp(ptr,"rotate_degree")==0)
	{
		for(int h=0;h<population;h++)
		{
			if(ary[h].rotate_degree > metric_value)
			{
				which_chro = h;
				metric_value = ary[h].rotate_degree;
			}
		}
	}
	else if(strcmp(ptr,"max_step_len")==0)
	{
		for(int h=0;h<population;h++)
		{
			if(ary[h].max_step_len > metric_value)
			{
				which_chro = h;
				metric_value = ary[h].max_step_len;
			}
		}
	}
    return which_chro;
}
void CSEX(chromosome* next_generation, chromosome& p1, chromosome& p2,
         int set_count, float** dist_vector, int& cnt,
         X_Y* sets_vector, coordinate* ary, node_list* cache)
{
    int flag = 0, num = 0, len = 1, loop = 1;
    int a_ahead, b_ahead, left_a, right_a, left_b, right_b;
    struct interval
    {
        int p1_l, p1_r;
        int p2_l, p2_r;
    };
    interval* common_tour = new interval [max_common_tours];

    while(loop < set_count-1)//find all common subtours
    {
        a_ahead = loop;
        for(int m=1;m<set_count;m++) if(p2.chromo[m] == p1.chromo[loop]) {b_ahead = m;break;}
        if(flag == 0)
        {
            left_a = a_ahead;
            left_b = right_b = b_ahead;
            flag = 1;
            loop++;
        }
        else
        {
            if(abs(right_b-b_ahead) == 1)
            {
                right_a = a_ahead;
                right_b = b_ahead;
                loop++;
                len++;
            }
            else
            {
                if(len > 1 && right_b < left_b)//only reverse common tours
                {
                    swap(right_b,left_b);
                    common_tour[num].p1_l = left_a;
                    common_tour[num].p1_r = right_a;
                    common_tour[num].p2_l = left_b;
                    common_tour[num].p2_r = right_b;
                    num++;//總共common tours
                    if(num >= max_common_tours) {/*cout<<"too much common tours!"<<endl;*/break;}
                }
                flag = 0; 
                len = 1;
            }
        }
    }

    //generate 2*2^num-2 offsprings
    int off_num = 1;
    for(int u=0;u<num;u++) off_num *= 2;
    off_num = 2*off_num;
//cout<<"  "<<num<<"  ";//getchar();
    chromosome* off = new chromosome [off_num];
    for(int i=0;i<off_num/2;i++)
    {
        off[i].chromo = new int [set_count];
        for(int j=0;j<set_count;j++) off[i].chromo[j] = p1.chromo[j];
    }
    for(int i=off_num/2;i<off_num;i++)
    {
        off[i].chromo = new int [set_count];
        for(int j=0;j<set_count;j++) off[i].chromo[j] = p2.chromo[j];
    }

    //all combinations
    int choose[num];
    for(int ct=0;ct<off_num/2;ct++)
    {
        int tmp = ct;
        int ite = num;
        while(ite != 0)
        {
              //cout<<tmp%2<<" ";
            choose[ite-1] = tmp%2;
            tmp /= 2;
            ite--;
        }
        for(int f=0;f<num;f++)
        {
            if(choose[f] == 0)//B
            {
                int pivot = common_tour[f].p1_l;
                for(int d=common_tour[f].p2_l; d<=common_tour[f].p2_r; d++)
                  off[ct+off_num/2].chromo[d] = p1.chromo[pivot++];
            }
            else//A
            {
                int pivot = common_tour[f].p2_l;
                for(int d=common_tour[f].p1_l; d<=common_tour[f].p1_r; d++)
                  off[ct].chromo[d] = p2.chromo[pivot++];
            }
        }
    }
    delete [] common_tour;

    //choose best chromosome
    //fitness(off_num,off,set_count,dist_vector);//-----------------------1
    fitness_CP(off_num,off,set_count,sets_vector,ary,cache);//---2
//for(int v=0;v<off_num;v++) cout<<off[v].fitness_val<<" * ";
    int best1 = min_value_chromosome(off_num,off,"fitness_val");
    next_generation[cnt].chromo = new int [set_count];
    for(int y=0;y<set_count;y++) next_generation[cnt].chromo[y] = off[best1].chromo[y];
    next_generation[cnt].fitness_val = off[best1].fitness_val;
    cnt++;

    for(int i=0;i<off_num;i++) delete [] off[i].chromo;
    delete [] off;
}
void SCX(chromosome* next_generation, chromosome& p1, chromosome& p2,
        int set_count, float** dist_vector, int& cnt)//每次產生一個chromosome
{
    int node = 0, count = 0;//start set
    int flag[set_count];//index as set number
    for(int i=0;i<set_count;i++) flag[i] = 0;
    int* child = new int [set_count];
    flag[node] = 1;
    child[count++] = node;

    while(count != set_count)
    {
        int node1, node2;
        int pos1, pos2;

        node1 = node2 = set_count;//impossible set number
        for(int m=0;m<set_count;m++) if(p1.chromo[m] == node) {pos1 = m;break;}
        for(int i=pos1+1;i<set_count;i++) if(flag[p1.chromo[i]] == 0)
        {
            node1 = p1.chromo[i];
            break;
        }
        if(node1 == set_count) for(int i=1;i<set_count;i++) if(flag[i] == 0) {node1 = i;break;}

        for(int m=0;m<set_count;m++) if(p2.chromo[m] == node) {pos2 = m;break;}
        for(int i=pos2+1;i<set_count;i++) if(flag[p2.chromo[i]] == 0)
        {
            node2 = p2.chromo[i];
            break;
        }
        if(node2 == set_count) for(int i=1;i<set_count;i++) if(flag[i] == 0) {node2 = i;break;}
        if(dist_vector[node][node1] < dist_vector[node][node2])
        {
            child[count++] = node = node1;
            flag[node1] = 1;
        }
        else
        {
            child[count++] = node = node2;
            flag[node2] = 1;
        }
    }

    next_generation[cnt++].chromo = child;
    child = NULL;
    fitness(1,next_generation+cnt-1,set_count,dist_vector);
    /*cout<<next_generation[cnt-1].fitness_val<<endl;
    cout<<p1.fitness_val<<"  "<<p2.fitness_val<<endl;*/
    float MIN = min(p1.fitness_val,p2.fitness_val);
    if(MIN < next_generation[cnt-1].fitness_val)
    {
        if(p1.fitness_val < p2.fitness_val)
          for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p1.chromo[u];
        else for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p2.chromo[u];
    }
    /*fitness(1,next_generation+cnt-1,set_count,dist_vector);
    cout<<next_generation[cnt-1].fitness_val<<endl;getchar();*/
}
void SCX_CP(chromosome* next_generation, chromosome& p1, chromosome& p2,
        int set_count, int& cnt, X_Y* sets_vector, coordinate* ary,
        node_list* cache, int method)//每次產生一個chromosome
{
    int node = 0, count = 0;//start set
    int flag[set_count];//index as set number
    for(int i=0;i<set_count;i++) flag[i] = 0;
    int* child = new int [set_count];
    flag[node] = 1;
    child[count++] = node;
    float start_x = sets_vector[node].x_axis;
    float start_y = sets_vector[node].y_axis;

    while(count != set_count)
    {
        int node1, node2;
        int pos1, pos2;

        node1 = node2 = set_count;//impossible set number
        for(int m=0;m<set_count;m++) if(p1.chromo[m] == node) {pos1 = m;break;}
        for(int i=pos1+1;i<set_count;i++) if(flag[p1.chromo[i]] == 0)
        {
            node1 = p1.chromo[i];
            break;
        }
        if(node1 == set_count)
        {
            if(method == 0) for(int i=1;i<set_count;i++) {if(flag[i] == 0) {node1 = i;break;}}
            else for(int i=1;i<set_count;i++) {if(flag[p1.chromo[i]] == 0) {node1 = p1.chromo[i];break;}}
        }

        for(int m=0;m<set_count;m++) if(p2.chromo[m] == node) {pos2 = m;break;}
        for(int i=pos2+1;i<set_count;i++) if(flag[p2.chromo[i]] == 0)
        {
            node2 = p2.chromo[i];
            break;
        }
        if(node2 == set_count)
        {
            if(method == 0) for(int i=1;i<set_count;i++) {if(flag[i] == 0) {node2 = i;break;}}
            else for(int i=1;i<set_count;i++) {if(flag[p2.chromo[i]] == 0) {node2 = p2.chromo[i];break;}}
        }
        //==========
        int j = node1;//next city
        float Min = INT_MAX;
        float XX, YY;
        int n1_n2;
        for(int jj=0;jj<2;jj++)
        {
            //find線段和圓的交點
            X_Y cross_list[sets_vector[j].size];//所有可能性
            int int_number = find_cross_node(start_x,start_y,sets_vector[j],cross_list,ary);//可能無cross
            //check valid, choose best in cross_list
            for(int k=0;k<int_number;k++)
            {
                if(check_node_vaild(cross_list[k].x_axis,cross_list[k].y_axis,sets_vector[j],ary))//檢查某點是否在共同交集區
                {
                    float dist = Euclideandistance(cross_list[k].x_axis,cross_list[k].y_axis,start_x,start_y);
                    if(dist < Min)
                    {
                        XX = cross_list[k].x_axis;
                        YY = cross_list[k].y_axis;
                        Min = dist;
                        n1_n2 = jj;
                    }
                }
            }
            //compare with set_node
            for(int s=0;s<cache[j].set_node_sz;s++)
            {
                float dist = Euclideandistance(cache[j].set_node_x[s],cache[j].set_node_y[s],start_x,start_y);
                if(dist < Min)
                {
                    XX = cache[j].set_node_x[s];
                    YY = cache[j].set_node_y[s];
                    Min = dist;
                    n1_n2 = jj;
                }
            }
            j = node2;
        }
        //==========
        if(n1_n2 == 0)//node1
        {
            child[count++] = node = node1;
            flag[node1] = 1;
        }
        else
        {
            child[count++] = node = node2;
            flag[node2] = 1;
        }
        start_x = XX;
        start_y = YY;
    }

    next_generation[cnt++].chromo = child;
    child = NULL;
    fitness_CP(1,next_generation+cnt-1,set_count,sets_vector,ary,cache);

    float MIN = min(p1.fitness_val,p2.fitness_val);
    if(MIN < next_generation[cnt-1].fitness_val)
    {
        if(p1.fitness_val < p2.fitness_val)
        {
            for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p1.chromo[u];
            next_generation[cnt-1].fitness_val = p1.fitness_val;
        }
        else
        {
            for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p2.chromo[u];
            next_generation[cnt-1].fitness_val = p2.fitness_val;
        }
    }
}
void CX(chromosome* next_generation, chromosome& p1, chromosome& p2, int set_count, int& cnt,
        X_Y* sets_vector, coordinate* ary, node_list* cache)
{
    int* child = new int [set_count];
    for(int i=0;i<set_count;i++) child[i] = p2.chromo[i];

    int next;
    int start = 1;//index in parent2
    int end = p1.chromo[start];//final city
    while(1)
    {
        child[start] = p1.chromo[start];
        next = p2.chromo[start];
        if(next == end) break;
        for(int i=1;i<set_count;i++) if(p1.chromo[i] == next) {start = i;break;}
    }

    next_generation[cnt++].chromo = child;
    child = NULL;
    fitness_CP(1,next_generation+cnt-1,set_count,sets_vector,ary,cache);

    float MIN = min(p1.fitness_val,p2.fitness_val);
    if(MIN < next_generation[cnt-1].fitness_val)
    {
        if(p1.fitness_val < p2.fitness_val)
        {
            for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p1.chromo[u];
            next_generation[cnt-1].fitness_val = p1.fitness_val;
        }
        else
        {
            for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p2.chromo[u];
            next_generation[cnt-1].fitness_val = p2.fitness_val;
        }
    }
}
void simple(chromosome* next_generation, chromosome& p1, chromosome& p2, int set_count, int& cnt,
        X_Y* sets_vector, coordinate* ary, node_list* cache)
{
        int point = rand()%(set_count-3)+1;
        chromosome tmp1, tmp2;
        tmp1.chromo = new int [set_count];
        tmp2.chromo = new int [set_count];
        int pos = 0;
        for(int i=0;i<=point;i++)
        {
            tmp1.chromo[pos] = p1.chromo[i];
            tmp2.chromo[pos] = p2.chromo[i];
            pos++;
        }
        for(int i=point+1;i<set_count;i++)
        {
            tmp1.chromo[pos] = p2.chromo[i];
            tmp2.chromo[pos] = p1.chromo[i];
            pos++;
        }
        if(check_valid(tmp1,set_count))
        {
            fitness_CP(1,&tmp1,set_count,sets_vector,ary,cache);
            fitness_CP(1,&tmp2,set_count,sets_vector,ary,cache);

            if(tmp1.fitness_val < tmp2.fitness_val)
            {
                next_generation[cnt++] = tmp1;
                delete [] tmp2.chromo;
            }
            else
            {
                next_generation[cnt++] = tmp2;
                delete [] tmp1.chromo;
            }
        }
        else
        {
            delete [] tmp1.chromo;
            delete [] tmp2.chromo;
            next_generation[cnt].chromo = new int [set_count];
            next_generation[cnt].fitness_val = INT_MAX;
            cnt++;
        }

        float MIN = min(p1.fitness_val,p2.fitness_val);
        if(MIN < next_generation[cnt-1].fitness_val)
        {
            if(p1.fitness_val < p2.fitness_val)
            {
                for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p1.chromo[u];
                next_generation[cnt-1].fitness_val = p1.fitness_val;
            }
            else
            {  
                for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p2.chromo[u];
                next_generation[cnt-1].fitness_val = p2.fitness_val;
            }
        }
}
void AFDX_initial(int set_count, float** dist_vector, double** adap_degree, double& threshold, double* up_bound, double& low_bound)
{
	//Calculate every farest distant of all set points.
	//set individual upper bounds
	for(int tmp=0;tmp<set_count;tmp++) up_bound[tmp] = 0.0;
	
	for(int i=0;i<set_count;i++)
	{
		double max = 0;
		for(int j=0;j<set_count;j++) if(max < dist_vector[i][j]) max = dist_vector[i][j];
		for(int j=0;j<set_count;j++) 
		{
			adap_degree[i][j] = 1-(dist_vector[i][j]/max);
			if(adap_degree[i][j] > up_bound[i] && adap_degree[i][j] != 1) up_bound[i] = adap_degree[i][j];
		}
	}
	//set a common lower bound
	low_bound = AFDX_LOW_BOUND;
	double tmp_low_bound[set_count];
	for(int i=0;i<set_count;i++)
	{
		double max = 0;
		for(int j=0;j<set_count;j++) 
			if(adap_degree[i][j] > max && adap_degree[i][j] < up_bound[i]) max = adap_degree[i][j];
		tmp_low_bound[i] = max;
	}
	double min = INT_MAX;
	for(int i=0;i<set_count;i++) if(tmp_low_bound[i] < min) min = tmp_low_bound[i];
	low_bound = min;/**/
	//cout<<endl<<"original lower bound: "<<low_bound<<endl;
	threshold = (1.0+low_bound)/2;//initial threshold
}
void AFDX(chromosome* next_generation, chromosome& p1, chromosome& p2,
        int set_count, int& cnt, X_Y* sets_vector, coordinate* ary, node_list* cache, 
		double** adap_degree, double threshold, double* up_bound)//每次產生一個chromosome
{
	int next = rand()%(set_count-1) + 1;//next is set ID
	int* child = new int [set_count];
	int index = 0;
	child[index++] = 0;//start from 0
	child[index++] = next;
	bool flag[set_count];
	for(int tmp=0;tmp<set_count;tmp++) flag[tmp] = true;
	flag[0] = flag[next] = false;//index means set
	
	for(int tmp=2;tmp<set_count;tmp++)
	{	
		int ind1, ind2;
		for(int i=1;i<set_count;i++)
			if(p1.chromo[i] == next) {ind1 = i;break;}
		for(int i=1;i<set_count;i++)
			if(p2.chromo[i] == next) {ind2 = i;break;}
		//check neighbors
		int neighbors[4];//save total 'next' in parents
		int count = 0;
		bool find_or_not = false;
		if(ind1 == 1) neighbors[count++] = p1.chromo[ind1+1];
		else if(ind1 == set_count-1) neighbors[count++] = p1.chromo[set_count-2];
		else 
		{
			neighbors[count++] = p1.chromo[ind1+1];
			neighbors[count++] = p1.chromo[ind1-1];
		}
		if(ind2 == 1) neighbors[count++] = p2.chromo[ind2+1];
		else if(ind2 == set_count-1) neighbors[count++] = p2.chromo[set_count-2];
		else 
		{
			neighbors[count++] = p2.chromo[ind2+1];
			neighbors[count++] = p2.chromo[ind2-1];
		}
		double tmp_thr = threshold;
		if(threshold > up_bound[next]) threshold = up_bound[next];
		
		for(int n=0;n<count;n++)
		{			
			if(flag[neighbors[n]] && adap_degree[next][neighbors[n]]>=threshold)//從親代臨近中選
			{
				child[index++] = next = neighbors[n];
				flag[neighbors[n]] = false;
				find_or_not = true;
				break;
			}
		}
		
		if(find_or_not == false)
		{
			for(int h=1;h<set_count;h++)//h is set ID
			if(flag[h] && adap_degree[next][h]>=threshold)//從符合threshold中隨意選
			{
				child[index++] = next = h;
				flag[h] = false;
				find_or_not = true;
				break;
			}
		}
		
		if(find_or_not == false)
		{		
			double max = -1;
			int which;
			for(int h=1;h<set_count;h++)//h is set ID
			if(flag[h] && adap_degree[next][h]>max)//select another closest (歸屬度最大) 
			{
				max = adap_degree[next][h];
				which = h;
			}
			child[index++] = next = which;
			flag[which] = false;
		}
		threshold = tmp_thr;//restore threshold
	}
	next_generation[cnt++].chromo = child;
	child = NULL;
    fitness_CP(1,next_generation+cnt-1,set_count,sets_vector,ary,cache);
	
	float MIN = min(p1.fitness_val,p2.fitness_val);
    if(MIN < next_generation[cnt-1].fitness_val)
    {
        if(p1.fitness_val < p2.fitness_val)
        {
            for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p1.chromo[u];
            next_generation[cnt-1].fitness_val = p1.fitness_val;
        }
        else
        {
            for(int u=0;u<set_count;u++) next_generation[cnt-1].chromo[u] = p2.chromo[u];
            next_generation[cnt-1].fitness_val = p2.fitness_val;
        }
    }
}
int evaluate_parents(int set_count, chromosome* all_chromosome, int parent1, int parent2)
{
	bool** matrix = new bool* [set_count];
	for(int i=0;i<set_count;i++) 
	{
		matrix[i] = new bool [set_count];
		for(int j=0;j<set_count;j++) matrix[i][j] = false;
	}
	//============
	int tmp1, tmp2, distance_value = 0;
	for(int k=0;k<set_count-1;k++) 
	{
		tmp1 = all_chromosome[parent1].chromo[k];
		tmp2 = all_chromosome[parent1].chromo[k+1];
		matrix[tmp1][tmp2] = matrix[tmp2][tmp1] = true;
	}
	matrix[0][tmp2] = matrix[tmp2][0] = true;
	//============
	for(int k=0;k<set_count-1;k++) 
	{
		tmp1 = all_chromosome[parent2].chromo[k];
		tmp2 = all_chromosome[parent2].chromo[k+1];
		if(matrix[tmp1][tmp2] == true) distance_value++;
	}
	if(matrix[0][tmp2] == true) distance_value++;
	//============
	for(int i=0;i<set_count;i++) delete [] matrix[i];
	delete [] matrix;
	
	return distance_value;
}
float evaluate_population_diversity_SD(int population, int set_count, chromosome* all_chromosome)
{
	int** matrix = new int* [set_count];
	for(int i=0;i<set_count;i++) 
	{
		matrix[i] = new int [set_count];
		for(int j=0;j<set_count;j++) matrix[i][j] = 0;
	}
	//============
	int tmp1, tmp2;
	float total = 0, SD_average = (float)population/(set_count-1);
	//set_count*population/set_count*(set_count-1);
	for(int num=0;num<population;num++)
	{
		for(int k=0;k<set_count-1;k++) 
		{
			tmp1 = all_chromosome[num].chromo[k];
			tmp2 = all_chromosome[num].chromo[k+1];
			matrix[tmp1][tmp2] += 1;
			//matrix[tmp2][tmp1] = matrix[tmp1][tmp2];
		}
		matrix[tmp2][0] += 1;
		//matrix[tmp2][0] = matrix[0][tmp2];
	}
	//============
	for(int i=0;i<set_count;i++) 
	{
		for(int j=0;j<set_count;j++) total += matrix[i][j]*matrix[i][j];
		delete [] matrix[i];
	}
	delete [] matrix;
	return sqrt(total/(set_count*(set_count-1))-SD_average*SD_average);
}
chromosome* GA_loop(int population, chromosome* all_chromosome, int set_count, float** dist_vector,
                    X_Y* sets_vector, coordinate* ary, node_list* cache,const int met1,const int met2,
					double** adap_degree, double threshold, double* up_bound)
{ 
    chromosome* next_generation = new chromosome [population];
    int cnt = 0;

    while(cnt != population)//至少還可容納一個chromosome
    {
int inner_loop_count = 10;
int candidate;
int mini_dist = INT_MAX; 
		//select two chromosomes as parents---------1
        int p1, p2, parent1, parent2;
        p1 = rand()%population;
        do{p2 = rand()%population;}while(p1 == p2);
        float a = all_chromosome[p1].fitness_val;
        float b = all_chromosome[p2].fitness_val;
		/*	float rotate_degree;
	float max_step_len;*/
        if(a >= b) parent1 = p2;
        else parent1 = p1;
inner_loop:	
        do{p1 = rand()%population;}while(p1 == parent1);
        do{p2 = rand()%population;}while(p1 == p2 || p2 == parent1);
        a = all_chromosome[p1].fitness_val;
        b = all_chromosome[p2].fitness_val;
        if(a >= b) parent2 = p2;
        else parent2 = p1;
//cout<<"parent: "<<parent1<<" "<<parent2<<endl;getchar();

	/*if(inner_loop_count > 0) 
	{
		int tmp = evaluate_parents(set_count,all_chromosome,parent1,parent2);
		if(mini_dist > tmp) {mini_dist = tmp;candidate = parent2;}
		inner_loop_count--;
		goto inner_loop;
	}
parent2 = candidate;*/
        //crossover,每次產生一個新chromosome,可能會重複
        if((double)rand()/(double)RAND_MAX <= crossover_rate)
        {
            switch(met1)
            {
            case 1:
            CSEX(next_generation,all_chromosome[parent1],all_chromosome[parent2],set_count,dist_vector,cnt,sets_vector,ary,cache);
            break;
            case 2:
            CX(next_generation,all_chromosome[parent1],all_chromosome[parent2],set_count,cnt,sets_vector,ary,cache);//cycle crossover
            break;
            case 3://MSCX
            SCX_CP(next_generation,all_chromosome[parent1],all_chromosome[parent2],set_count,cnt,sets_vector,ary,cache,1);//SCX:0, MSCX:1
            break;
            case 4://SCX
            SCX_CP(next_generation,all_chromosome[parent1],all_chromosome[parent2],set_count,cnt,sets_vector,ary,cache,0);//SCX:0, MSCX:1
            break;
            case 5:
            simple(next_generation,all_chromosome[parent1],all_chromosome[parent2],set_count,cnt,sets_vector,ary,cache);
            break;
            case 6:
            AFDX(next_generation,all_chromosome[parent1],all_chromosome[parent2],set_count,cnt,sets_vector,ary,cache,
			adap_degree,threshold,up_bound);
            break;			
            //SCX(next_generation,all_chromosome[parent1],all_chromosome[parent2],set_count,dist_vector,cnt);//without CP
            //EAX();
            }
        }
        else//直接複製原來兩個chromosome
        {
            next_generation[cnt].chromo = new int [set_count];
            //memcpy(next_generation[cnt].chromo,all_chromosome[parent1].chromo,set_count*sizeof(int));
            for(int y=0;y<set_count;y++) next_generation[cnt].chromo[y] = all_chromosome[parent1].chromo[y];
            next_generation[cnt].fitness_val = all_chromosome[parent1].fitness_val;
            cnt++;

            if(cnt == population) break;
            next_generation[cnt].chromo = new int [set_count];
            //memcpy(next_generation[cnt].chromo,all_chromosome[parent2].chromo,set_count*sizeof(int));
            for(int y=0;y<set_count;y++) next_generation[cnt].chromo[y] = all_chromosome[parent2].chromo[y];
            next_generation[cnt].fitness_val = all_chromosome[parent2].fitness_val;
            cnt++;
        }
    }
//fitness_CP(population,next_generation,set_count,sets_vector,ary,cache);
    //mutation
    for(int i=0;i<population;i++)
    {
        if((double)rand()/(double)RAND_MAX <= mutation_rate)
        {
            switch(met2)
            {
            case 1:
            _2_Opt(next_generation[i],set_count,dist_vector,sets_vector,ary,cache);//-----------1
            break;
            case 2:
            _Best_2_Opt(next_generation[i],set_count,dist_vector,sets_vector,ary,cache);//------2
            break;
            case 3:
            IOpt(next_generation[i],set_count,dist_vector,sets_vector,ary,cache);//-------------3
            break;
            case 4://Reciprocal exchange mutation (EM)
            EM(next_generation[i],set_count,dist_vector,sets_vector,ary,cache);//---------------4
            break;
            }
        }
        delete [] all_chromosome[i].chromo;
    }
    delete [] all_chromosome;
//for(int i=0;i<population;i++) check_valid(next_generation[i],set_count);
    return next_generation;
}
int MaxIntersection(X_Y s, int* u, bool wipe)
{
    int x = 0;
    for(int i=0;i<s.size;i++) if(u[s.ary[i]])
    {
        x++;
        if(wipe) u[s.ary[i]] = 0;
    }
    return x;
}
void Spanning_Set_Covering(int u_sz, int s_sz, X_Y* joint_vector, float** dist_vector_joint)
{
    int* universe = new int[u_sz];
    for(int i=0;i<u_sz;i++) universe[i] = 1;//initialization, full

    int cost = 0, cur = 0, WhichSet, StartSet = 0, EndSet = StartSet;
    double max = INT_MIN, alpha, Accumulated_Cost = -1;
    /*if(StartSet == 0)//挑一個涵蓋最多sensor的set當起始點
    {
        for(int i=0;i<s_sz;i++) if(joint_vector[i].size > max)//set number is from 0
        {
            max = joint_vector[i].size;
            StartSet = i;
        }
    }*/
    double const_para = 0.5;
    while(cur != u_sz)
    {
          //cout<<"Start from set: "<<StartSet+1<<endl;
        WhichSet = max = 0;
        for(int i=0;i<s_sz;i++) if(joint_vector[i].size != 0)//not yet be picked
        {
            double sz = MaxIntersection(joint_vector[i], universe, false);
            alpha = sz/dist_vector_joint[StartSet][i];                //----1
            //alpha = sz;                                             //----2
            //alpha = sz + 2500*(1/dist_vector_joint[StartSet][i]);   //----3
            //alpha = (1/dist_vector_joint[StartSet][i]);             //----4
            //cout<<"==="<<1/(dist_vector_joint[StartSet][i])<<"==";getchar();
              //cout<<"Alpha of set "<<i+1<<": "<<alpha<<endl;
            if(alpha > max)
            {
                max = alpha;
                WhichSet = i;
            }
        }
          //cout<<"Pick set: "<<WhichSet<<endl;
        Accumulated_Cost += dist_vector_joint[StartSet][WhichSet];
        StartSet = WhichSet;
          //cout<<"Accumulated cost: "<<Accumulated_Cost<<endl;
        cur += MaxIntersection(joint_vector[WhichSet], universe, true);
        joint_vector[WhichSet].size = 0;//picked
          //cout<<"Covered sensors: "<<cur<<endl<<"--------------------------"<<endl;
        cost++;
    }
    cout<<" "<<setw(10)<<Accumulated_Cost+dist_vector_joint[StartSet][EndSet]<<" ";
    cout<<""<<setw(10)<<cost<<" ";
    delete [] universe;
}
void print_coordinates(X_Y* sets_vector, int set_count, X_Y* ext_set, int ext_count)
{
    for(int i=0;i<set_count;i++)
    cout<<"("<<sets_vector[i].x_axis<<", "<<sets_vector[i].y_axis<<")"<<endl;
    cout<<"=============================================="<<endl;
    for(int i=0;i<ext_count;i++)
    cout<<"("<<ext_set[i].x_axis<<", "<<ext_set[i].y_axis<<")"<<endl;
    cout<<"=============================================="<<endl;
} 

int main(int argc, char** argv)
{
footprint.open("segment.txt", ios::out);
print_waypoint.open("waypoint.txt", ios::out);
    cout<<"==> [crossover] [mutation] [initial_generation]"<<endl;
    int method1, method2, method3;
    if(argc > 1)
    {
        method1 = atoi(argv[1]);
        method2 = atoi(argv[2]);
		method3 = atoi(argv[3]);
        cout<<argc-1<<" parameters => "<<method1<<"::"<<method2<<"::"<<method3<<endl;
    }

    int RandomNode=R_NUM, Radius=RADIUS, loop = RUN_TIMES_AVG;
    float original = 0, average = 0, avg_num = 0, best_ite = INT_MAX;
	float greedy_avg = 0, greedy_rotate_degree = 0, greedy_max_step_len = 0, conv = 0;
	float greedy_avg_LLA = 0, greedy_rotate_degree_LLA = 0, greedy_max_step_len_LLA = 0;
	float SD_avg = 0, rotate_degree = 0, max_step_len = 0;
	float initial_population_best = 0, initial_rotate = 0, initial_max_step = 0;
	float avg_after_lookahead = 0, rotate_after_lookahead = 0, max_step_after_lookahead = 0; 
	float LLA_success_count = 0, greedy_LLA_success_count = 0;
	float LLA_hit_ratio_greedy = 0, LLA_hit_ratio = 0;
	
    srand(time(NULL));//set seed
    for(int p=0;p<loop;p++)
    {
        coordinate* ary = new coordinate [RandomNode];//所有點
        float** vector = new float* [RandomNode];//點之間距離
        for(int i=0;i<RandomNode;i++)
        {
            vector[i] = new float [RandomNode];
            vector[i][i] = 0;
        }
        random(ary,RandomNode);//產生所有點座標
        //cout<<endl;
        calculate_node_dist(ary,vector,RandomNode);//計算點之間距離
        int set_count = cluster(ary,vector,RandomNode);
		/*************************************************/
		/*if(set_count!=73)
		{
			for(int i=0;i<RandomNode;i++)
			{
				if(ary[i].ptr != NULL)
				{
					node* tmp = ary[i].ptr->ptr;
					delete ary[i].ptr;//delete info
					while(tmp != NULL)
					{
						node* tmp2 = tmp;
						tmp = tmp->next;
						ary[tmp2->sensor].ptr = NULL;
						delete tmp2;
					}
				}
				delete [] vector[i];
			}
			delete [] ary;
			delete [] vector;
			p--;
			continue;
		}*/
		/*************************************************/
		avg_num += set_count;
        cout<<RandomNode<<" "<<setw(4)<<set_count<<" "<<setw(4)<<(float)(RandomNode-set_count)/RandomNode*100<<"%";
//cout<<endl;continue;
        X_Y* sets_vector = new X_Y [set_count];
        float** dist_vector = new float* [set_count];
        for(int i=0;i<set_count;i++)
        {
            dist_vector[i] = new float [set_count];
            dist_vector[i][i] = 0;
        }
        calculate_set_dist(sets_vector,set_count,dist_vector,RandomNode,ary);//計算set之間距離, dist_vector
        //=============================
        //---------------
        double** adap_degree;
        double threshold;
        double low_bound; 
        double* up_bound;//find the largest except 1
        if(method1 == 6)//AFDX crossover
        {
        	adap_degree = new double* [set_count];
        	up_bound = new double [set_count];//find the largest except 1
        	for(int i=0;i<set_count;i++) adap_degree[i] = new double [set_count];
        	AFDX_initial(set_count,dist_vector,adap_degree,threshold,up_bound,low_bound);
        }
        //---------------
        node_list* cache = new node_list [set_count];//index為set編號
        calculate_set_node(cache,ary,sets_vector,set_count);
        
        chromosome first_chromo;
        first_chromo.chromo = new int [set_count];
        float first_cost = greedy(set_count,dist_vector,first_chromo.chromo);//first chromosome is decided by greedy algorithm
        cout<<"  greedy::"<<setw(8)<<first_cost<<" ";
		original += first_cost;
        chromosome* all_chromosome;
        /*chromosome* all_chromosome2;
        chromosome* all_chromosome3;
        chromosome* all_chromosome4;
        
        chromosome first_chromo2;
        first_chromo2.chromo = new int [set_count];
        memcpy(first_chromo2.chromo,first_chromo.chromo,set_count*sizeof(int));
        chromosome first_chromo3;
        first_chromo3.chromo = new int [set_count];
        memcpy(first_chromo3.chromo,first_chromo.chromo,set_count*sizeof(int));
        chromosome first_chromo4;
        first_chromo4.chromo = new int [set_count];
        memcpy(first_chromo4.chromo,first_chromo.chromo,set_count*sizeof(int));*/
        
        //initialization
        int population = population_sz;
        int remain = migration_interval;
        chromosome migration_original;
        migration_original.fitness_val = INT_MAX;
        migration_original.chromo = new int [set_count];
         
        switch(method3)
        {
        	case 1:
        	    all_chromosome = CGA(first_chromo.chromo,set_count,population);//------1
        	    break;
        	case 2:
        	    all_chromosome = random_initial(first_chromo.chromo,set_count,population);//------2
        	    break;
        	case 3:
        	    all_chromosome = balance_SD(set_count,population);//------3
        	    break;
        	case 4:
        	    all_chromosome = balance_SD_2(set_count,population);//------3
        	    break;
        }
        
        fitness_CP(population,all_chromosome,set_count,sets_vector,ary,cache);
        //for(int v=0;v<population;v++) cout<<all_chromosome[v].fitness_val<<" ";getchar();
        
		X_Y path[set_count];//存沒有LLA的waypoint
        fitness_CP(1,&first_chromo,set_count,sets_vector,ary,cache,path);//印出沒有LLA的waypoint
        greedy_avg += first_chromo.fitness_val;//after CP
        greedy_rotate_degree += first_chromo.rotate_degree;//after CP
        greedy_max_step_len += first_chromo.max_step_len;//after CP
        cout<<"  CP::"<<setw(8)<<first_chromo.fitness_val; 
        
        float LLA_hit = final_fitness_CP(0,&first_chromo,set_count,sets_vector,ary,path);//印出有LLA的waypoint
        greedy_LLA_success_count += LLA_hit;
		LLA_hit_ratio_greedy += LLA_hit/set_count;
		
		greedy_avg_LLA += first_chromo.fitness_val;//after CP LLA
        greedy_rotate_degree_LLA += first_chromo.rotate_degree;//after CP LLA
        greedy_max_step_len_LLA += first_chromo.max_step_len;//after CP LLA
        
        float SD_ = evaluate_population_diversity_SD(population,set_count,all_chromosome);
        SD_avg += SD_;
        cout<<"  SD::"<<setw(8)<<SD_<<"";
        //cout<<endl<<"=============="<<endl;
        
        	//int min_fitness_val = min_value_chromosome(population,all_chromosome,"fitness_val");
        	//initial_population_best += all_chromosome[min_fitness_val].fitness_val;
        		/*cout<<all_chromosome[min_fitness_val].fitness_val<<"  "
        		<<all_chromosome[min_fitness_val].rotate_degree<<"  "
        		<<all_chromosome[min_fitness_val].max_step_len<<endl;*/
        		
        	//int min_rotate = min_value_chromosome(population,all_chromosome,"rotate_degree");
        	//initial_rotate += all_chromosome[min_rotate].rotate_degree;
        		/*cout<<all_chromosome[min_rotate].fitness_val<<"  "
        		<<all_chromosome[min_rotate].rotate_degree<<"  "
        		<<all_chromosome[min_rotate].max_step_len<<endl;*/
        		
        	//int min_max_step = min_value_chromosome(population,all_chromosome,"max_step_len");
        	//initial_max_step += all_chromosome[min_max_step].max_step_len;
        		/*cout<<all_chromosome[min_max_step].fitness_val<<"  "
        		<<all_chromosome[min_max_step].rotate_degree<<"  "
        		<<all_chromosome[min_max_step].max_step_len<<endl;*/
        		
        //cout<<"=============="<<endl;
        float converge_value = 0;
        int converge_cnt = convergence;
        float old_value = INT_MAX;//for AFDX
	    int best_chro;
        for(int kk=0;kk<generation;kk++)
        {
            all_chromosome = GA_loop(population,all_chromosome,set_count,dist_vector,sets_vector,ary,cache,method1,method2,
        	adap_degree,threshold,up_bound);
            fitness_CP(population,all_chromosome,set_count,sets_vector,ary,cache);
            best_chro = min_value_chromosome(population,all_chromosome,"fitness_val");
            float best_v1 = all_chromosome[best_chro].fitness_val;
            //cout<<kk<<" ("<<(int)best_v1<<") ";
        	
        		if(best_v1 < old_value) threshold += 0.05;
        		else threshold -= 0.05;
        		if(threshold > 1) threshold = 1;
        		else if(threshold < low_bound) threshold = low_bound;
        		old_value = best_v1;
        		//cout<<"   threshold::"<<setw(7)<<threshold<<" ";getchar();/*//AFDX*/
               
	    	if(converge_value == best_v1)
            {
                converge_cnt--;
                if(converge_cnt == 0)
                {
                    int total_evolution = kk-convergence+2;
	    			cout<<"   cvg::"<<setw(4)<<total_evolution<<" ";
                    conv += total_evolution;
                    break;
                }
            }
            else
            {
                converge_cnt = convergence;//reset
                converge_value = best_v1;
            }/*//calculate convergence*/
	    	
            /*remain--; 
            if(remain == 0)//migration, share best
            {
                remain = migration_interval;
                if((double)rand()/(double)RAND_MAX <= migration_rate)
                {
                    int give, receive;
                    give = rand()%4;
                    do{receive = rand()%4;}while(give == receive);
                    chromosome ptr;
                    switch(give)
                    {
                        case 0: ptr = all_chromosome[best_chro]; break;
                        case 1: ptr = all_chromosome2[best_chro2]; break;
                        case 2: ptr = all_chromosome3[best_chro3]; break;
                        case 3: ptr = all_chromosome4[best_chro4]; break;
                        default: ;
                    }
                    if(ptr.fitness_val < migration_original.fitness_val)
                    {
                        migration_original.fitness_val = ptr.fitness_val;
                        memcpy(migration_original.chromo,ptr.chromo,set_count*sizeof(int));
                    }
                    switch(receive)
                    {
                        case 0: memcpy(all_chromosome[w1].chromo,migration_original.chromo,set_count*sizeof(int)); break;
                        case 1: memcpy(all_chromosome2[w2].chromo,migration_original.chromo,set_count*sizeof(int)); break;
                        case 2: memcpy(all_chromosome3[w3].chromo,migration_original.chromo,set_count*sizeof(int)); break;
                        case 3: memcpy(all_chromosome4[w4].chromo,migration_original.chromo,set_count*sizeof(int)); break;
                        default: ;
                    }
                }
            }*/
        }
        float best_v1 = all_chromosome[best_chro].fitness_val;
        rotate_degree += all_chromosome[best_chro].rotate_degree;
        max_step_len += all_chromosome[best_chro].max_step_len;
        cout<<"  "<<best_v1<<"  ";
        average += best_v1;
        if(best_v1 < best_ite) best_ite = best_v1;
        
		fitness_CP(1,all_chromosome+best_chro,set_count,sets_vector,ary,cache,path);//印出沒有LLA的waypoint
        /*for(int z=0;z<set_count;z++)//印出原本CA給每個set的(x,y)
		{
            int ss = all_chromosome[best_chro].chromo[z];
            cout<<"("<<cache[ss].set_node_x[0]<<", "<<cache[ss].set_node_y[0]<<") ";
        }*/
		
		LLA_hit = final_fitness_CP(best_chro,all_chromosome,set_count,sets_vector,ary,path);//印出有LLA的waypoint
        LLA_success_count += LLA_hit;
		LLA_hit_ratio += LLA_hit/set_count;
		
		avg_after_lookahead += all_chromosome[best_chro].fitness_val;
        rotate_after_lookahead += all_chromosome[best_chro].rotate_degree;
        max_step_after_lookahead += all_chromosome[best_chro].max_step_len;

        if(method1 == 6)//AFDX crossover
        {
        	for(int i=0;i<set_count;i++) delete [] adap_degree[i];
        	delete [] adap_degree;
        	delete [] up_bound;
        }
        delete [] migration_original.chromo;
        for(int m=0;m<population;m++)
        {
            delete [] all_chromosome[m].chromo;
        }
        delete [] all_chromosome;
        delete [] first_chromo.chromo;
        //=============================

		/*X_Y* ext_set = new X_Y [max_ext_set];
		int ext_count = extend_set(ext_set,set_count,sets_vector,ary,Radius,vector);
		cout<<""<<setw(5)<<ext_count<<"";
	
		//print_coordinates(sets_vector,set_count,ext_set,ext_count);
	
		float** dist_vector_joint = new float* [set_count+ext_count];
		for(int i=0;i<set_count+ext_count;i++)
		{
			dist_vector_joint[i] = new float [set_count+ext_count];
			dist_vector_joint[i][i] = 1;//initialize to 1
		}
		X_Y* joint_vector = new X_Y [set_count+ext_count];
		merge(sets_vector,set_count,ext_set,ext_count,dist_vector_joint,joint_vector);//算set間距離,合併、delete兩vector
		Spanning_Set_Covering(RandomNode,set_count+ext_count,joint_vector,dist_vector_joint);
		free_memory2(dist_vector_joint,joint_vector,set_count+ext_count);*/

        free_memory(ary,vector,RandomNode,dist_vector,set_count,cache);
        //RandomNode+=50;
        cout<<endl; 
    }
	cout<<"greedy_nCP:          "<<original/loop<<endl;//greedy without CP
	cout<<"greedy_CP::          "<<greedy_avg/loop<<endl;//greedy with CP
	cout<<"greedy_LLA::         "<<greedy_avg_LLA/loop<<endl;
	cout<<"greedy_improved::    "<<(greedy_avg-greedy_avg_LLA)/greedy_avg*100<<"%"<<endl;//compare to greedy_CP, LLA
	
	cout<<"max_step_len::       "<<greedy_max_step_len/loop<<endl;
	cout<<"rotate_degree::      "<<greedy_rotate_degree/loop<<endl;
	cout<<"max_step_len_opt::   "<<greedy_max_step_len_LLA/loop<<endl;
	cout<<"rotate_degree_opt::  "<<greedy_rotate_degree_LLA/loop<<endl;
	
	cout<<"LLA_hit_count::      "<<greedy_LLA_success_count/loop<<endl
		<<"LLA_hit_ratio::      "<<LLA_hit_ratio_greedy/loop*100<<"%"<<endl;	

	cout<<"greedy_rotate::      "<<(greedy_rotate_degree-greedy_rotate_degree_LLA)/greedy_rotate_degree*100<<"%"<<endl;
	cout<<"greedy_max_step::    "<<(greedy_max_step_len_LLA-greedy_max_step_len)/greedy_max_step_len*100<<"%"<<endl;	
	cout<<endl;
	//======================================================================above is LLA in greedy
	//cout<<"initial_best::       "<<initial_population_best/loop<<endl;
    //cout<<"initial_max_step::   "<<initial_max_step/loop<<endl;
	//cout<<"initial_rotate::     "<<initial_rotate/loop<<endl;
	cout<<"GA_avg::             "<<average/loop<<endl;
	cout<<"GA_avg_opt::         "<<avg_after_lookahead/loop<<endl;
    
	cout<<"best::               "<<best_ite<<endl;
    cout<<"cluster::            "<<avg_num/loop<<endl;
	cout<<"improved::           "<<(greedy_avg-average)/greedy_avg*100<<"%"<<endl;//compare to greedy_CP, GA
	cout<<"improved_opt::       "<<(greedy_avg-avg_after_lookahead)/greedy_avg*100<<"%"<<endl;//compare to greedy_CP, GA, LLA
	cout<<"improve_diff::       "<<(average-avg_after_lookahead)/greedy_avg*100<<"%"<<endl;
	
	cout<<"converge::           "<<conv/loop<<endl;
	cout<<"SD_avg::             "<<SD_avg/loop<<endl;
	
	cout<<"max_step_len::       "<<max_step_len/loop<<endl;
	cout<<"rotate_degree::      "<<rotate_degree/loop<<endl;
	cout<<"max_step_len_opt::   "<<max_step_after_lookahead/loop<<endl;
	cout<<"rotate_degree_opt::  "<<rotate_after_lookahead/loop<<endl;
	
	cout<<"LLA_hit_count::      "<<LLA_success_count/loop<<endl
		<<"LLA_hit_ratio::      "<<LLA_hit_ratio/loop*100<<"%"<<endl;
	
	cout<<"rotate_improve::     "<<(rotate_degree-rotate_after_lookahead)/rotate_degree*100<<"%"<<endl;
	cout<<"max_increase::       "<<(max_step_after_lookahead-max_step_len)/max_step_len*100<<"%"<<endl;	

	footprint.close();
	print_waypoint.close();
    return 0;
}
