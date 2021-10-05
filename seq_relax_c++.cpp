#include <iostream>
#include <queue>
#include <fstream>
#include <vector>
#include <bitset>
#include <string>
#include <sstream>
#include "gurobi_c++.h"
#include "SetAndTestSchedule.h"
#include "Uheap.h"
#include "Q_ELEMENT.h"
#include "PriorityQueue.h"
#include "lb_compute_param.h"
using namespace std;

// type
typedef struct Q_ELEMENT E;
typedef deque<Q_ELEMENT> Qe;
typedef vector<int> Vi;
typedef vector<double> Vd;
typedef vector<Vi> VVi;
typedef bitset<100> B;
typedef vector<B> Vb;

// VVi D;

// // working time of test jobs , requirement jobs
// Vi P = {0, 21, 22, 14, 5, 19, 24, 13, 10, 30, 9, 8, 18, 7, 9, 26, 7, 20, 25, 29, 13, 4, 4, 17, 20, 15, 14, 11, 18, 6, 9, 5, 3, 12, 18, 5, 28, 8, 1, 20, 27, 6, 18, 25, 24, 8, 1, 27, 20, 13, 4, 22, 27, 14, 17, 9, 5, 13, 6, 16, 18};
// int R = 30;      // number of requirement jobs
// int J = 30;      // number of test jobs
// int N = R + J;   // total number of jobs
const int M = 1;

// Big Constant
const int Mw = 10;
const int My = 100000;

// hyper param
// 解放的constraint不連續，難relax
const int lambda = 50;


// VVi toDependecyMatrix(string s)
// {
//     VVi dep_matrix;
//     ifstream input;
//     input.open(s);
//     string str;
//     bool vec_init = true;
//     while(std::getline(input, str))
//     {        
//         bool i;
//         int set_up_job_idx = 0;        
//         if(vec_init)
//         {
//             stringstream ss(str);            
//             while(ss >> i)
//             {                
//                 Vi v;
//                 v.push_back(i);
//                 dep_matrix.push_back(v);
//             }
//             vec_init = false;
//         }
//         else
//         {
//             stringstream ss(str);
//             while(ss >> i)
//             {
//                 dep_matrix[set_up_job_idx].push_back(i);
//                 set_up_job_idx++;
//             }
//         }
        
//     }
//     input.close();

//     // printf("size: %d * %d\n", dep_matrix.size(), dep_matrix[1].size());
//     // for(auto it: dep_matrix)
//     // {
//     //     for(auto itt: it)
//     //     {
//     //         cout << itt << " ";
//     //     }
//     //     cout << endl;
//     // }

//     VVi matrix;    
//     matrix.assign(dep_matrix.size()+dep_matrix[0].size() + 1, Vi(dep_matrix.size()+dep_matrix[0].size() + 1));
//     for(int i = 1; i <= 30; i++)
//     {
//         for(int j = 31; j <= 60; j++)
//         {
//             matrix[i][j] = dep_matrix[i-1][j-31];
//         }
//     }
//     return matrix;
// }


// 用child在改的時候比較有效率，因為我們固定的是set-up jobs
int relax_lb(Vb &child, Vd &s, Vd &t)
{    
    /*  == Variables decalration == */
    GRBVar **S;
    GRBVar **C;
    GRBVar **W;
    GRBVar **Y;    

    int N = s.size() + t.size();
    int J = t.size();
    int R = s.size();
    Vd P(N+1);
    copy(t.begin()+1, t.end(), P.begin()+1);
    copy(s.begin()+1, s.end(), P.begin()+J+1);


    try {

    // Create an environment
    GRBEnv env = GRBEnv(true);
    env.set("LogFile", "mip1.log");
    env.start();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    // model.set("SolutionLimit", "1");

    /*  ================ Variable Population =============
        // S: start time, C_p: completion time
        // W: work on machine(or not)
        // Y: if Y_p[i][j] = 1, job i completes before j
        ==================================================    
    */
    S = new GRBVar*[M+1];
    for(int i = 1; i <= M; i++)
    {
        S[i] = new GRBVar[N+1];
        for(int j = 1; j <= N; j++)
            S[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, "S");
    }

    C = new GRBVar*[M+1];
    for(int i = 1; i <= M; i++)
    {
        C[i] = new GRBVar[N+1];
        for(int j = 1; j <= N; j++)
            C[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, "C");

    }

    W = new GRBVar*[M+1];
    for(int i = 1; i <= M; i++)
    {
        W[i] = new GRBVar[N+1];
        for(int j = 1; j <= N; j++)
            W[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "W");
    }

    Y = new GRBVar*[N+1];
    for(int i = 1; i <= N; i++)
    {
        Y[i] = new GRBVar[N+1];
        for(int j = 1; j <= N; j++)
            Y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "Y");
    }
       

    /*  ================ Set objective =================
    ================================================  */    
    GRBLinExpr objExpr = 0;
    for(int m = 1; m <= M; m++)
    {
        for(int j = 1; j <= J; j++)
        {
            // Add up all the completion time of test-jobs
            objExpr += C[m][j];
        }
    }
    for(int m = 1; m <= M; m++)
    {
        for(int j = 1; j <= J; j++)
        {
            for(int r = J+1; r <= N; r++)
            {
                if(child[r-J].test(j))
                {                    
                    string order_constr = "t"+ to_string(j) + "_r" + to_string(r-J)+ "_r_than_j";                    
                    objExpr += (1 - Y[r][j]) * lambda;
                    // model.addConstr(Y[r][j] == 1, order_constr);
                }
            }
        }
    }
    model.setObjective(objExpr, GRB_MINIMIZE);


    /*  ================ Set Constraints =================
    ================================================  */ 

    // j 個 test job 由 m 個 machine 共同分工完成 -- 完整性, (cmpltExpr = completion expression)

    GRBLinExpr cmpltExpr = 0;
    for(int j = 1; j <= J; j++)
    {
        for(int m = 1; m <= M; m++)
        {
            cmpltExpr += W[m][j];
        }
    }
    model.addConstr(cmpltExpr == J, "c1a");


    // 同一個 test job j 在整次工作中只能在一個machine上執行，且僅執行一次 -- 一次性
    for(int j = 1; j <= J; j++)
    {
        GRBLinExpr onceExpr = 0;
        for(int m = 1; m <= M; m++)
        {
            onceExpr += W[m][j];
        }
        string name = "onceExpr";
        name += to_string(j);
        model.addConstr(onceExpr == 1, name);
    }
    
    // 相依性
    for(int m = 1; m <= M; m++)
    {
        for(int j = 1; j <= J; j++)
        {
            for(int r = J+1; r <= N; r++)
            {
                if(child[r-J].test(j))
                {
                    string rt_constr = "t"+ to_string(j) + "_r" + to_string(r-J)+ "_dep";
                    string order_constr = "t"+ to_string(j) + "_r" + to_string(r-J)+ "_r_than_j";
                    model.addConstr(( 1 - W[m][j] ) * Mw + ( W[m][r] - 1 ) >= 0, rt_constr);
                    // model.addConstr(Y[r][j] == 1, order_constr);
                }
            }
        }
    }

    // 一定要有先後, order_clarity
    for(int i = 1; i <= N; i++)
    {
        for(int j = i+1; j <= N; j++)
        {
            string order_clar_constr = "job" + to_string(i) + "_job" + to_string(j) + "_order_clarity";
            model.addConstr(Y[i][j] + Y[j][i] == 1, order_clar_constr);
        }
    }

    // start time >= 0
    // completion time definitions
    for(int m = 1; m <= M; m++)
    {
        for(int n = 1; n <= N; n++)
        {
            string start_constr = "job" + to_string(n) + "_on_m" + to_string(m) + "_startTime";
            string completion_constr = "job" + to_string(n) + "_on_m" + to_string(m) + "_CompletionTime";
            model.addConstr(S[m][n] >= 0, start_constr);
            model.addQConstr(C[m][n] - ( S[m][n] + P[n] ) * W[m][n] == 0, completion_constr);
        }
    }

    // // 排他性, Exclusivity
    for(int m = 1; m <= M; m++)
    {
        for(int i = 1; i <= N; i++)
        {
            for(int j = i+1; j <= N; j++)
            {
                string ij_constr = "job" + to_string(i) + "_finishes_bf" + "_job" + to_string(j) + "_on_" + to_string(m);
                string ji_constr = "job" + to_string(j) + "_finishes_bf" + "_job" + to_string(i) + "_on_" + to_string(m);
                model.addConstr(S[m][j] - C[m][i] + My * ( 1 - Y[i][j] ) + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) >= 0, ij_constr);
                model.addConstr(S[m][i] - C[m][j] + My * Y[i][j] + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) >= 0);
            }
        }
    } 


    // Optimize model
    model.optimize();    
    return model.get(GRB_DoubleAttr_ObjVal);

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    // return lower bound
    return 0;
}

Lb_Compute_Param lbComputeParamSetUp(const Vb &_child, const Vd &_s, const Vd &_t, const Vi &setup_done)
{
    Vb child(_child);
    Vd s(_s);
    Vd t(_t);
    for(int i = 0; i < setup_done.size(); ++i)
    {        
        child[setup_done[i]].reset();
        s[setup_done[i]] = 0;
    }
    return Lb_Compute_Param(child, s, t);
}


int main()
{
    // all 1-based, init and ref
    int tmp, Sn, Tn;;
    cin >> Sn >> Tn;
    Vd s ,t;
    Vb prec;
    Vb child;
    SetAndTestSchedule st(Sn, Tn);
    child.assign(Sn+1, B(0));
    prec.assign(Tn+1, B(0));
    s.resize(Sn+1);
    t.resize(Tn+1);
    st.j.resize(Sn+Tn+1);
     

    // construct
    for(int i = 1; i <= Sn; i++){cin >> s[i]; } 
    for(int i = 1; i <= Tn; i++){cin >> t[i]; } 
    cin >> tmp;
    if(tmp != -1) cout << "input error" << endl;
    copy(s.begin()+1, s.end(), st.j.begin()+1);
    copy(t.begin()+1, t.end(), st.j.begin()+1+Sn);     



    for(int i = 1; i <= Sn; i++)
    {
        int b;
        for(int j = 1; j <= Tn; j++)
        {
            cin >> b;
            if(b) {
                prec[j].set(i);
                child[i].set(j);
            }
        }
    }

    // done test
    Vi setup_done = {1, 2};    
    Lb_Compute_Param param = lbComputeParamSetUp(child, s, t, setup_done);

    cout << "optimal: " << relax_lb(param.child, param.s, param.t) << endl;
}