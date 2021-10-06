/*  ================== seq_relax_c++ =======================
   
    - The order of set-up and testjob in general 
      when mixed is different with that of in permSet_outtree
    - the permSet_outtree it's (set_up, test_job)
    - while here it's (test_job, set_up)
    - so in this code, test_job has smaller index 
      when using the mixed version index
    - 這裡使用了兩種indexing方式
      跟permSet_outtree裡面的indexing剛好相反

    Indecies:
    There are two kinds of indexing in this code
    1. unmixed
        - In this version, set_up and test jobs has their 
        own indecies, starting from 1 respectively

    2. mixed
        - In this version, set_up and test jobs are being 
        seen as one kind and test jobs have smaller indecies 
        starting from 1
        - range: 
            - test job:     1 ~ J
            - set_up job: J+1 ~ N 

    ==================================================    */

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
const int lambda = 1000;


// 用child在改的時候比較有效率，因為我們固定的是set-up jobs
int relax_lb(const Vb &_child, const Vd &_s, const Vd &_t, const int &Sn, const int &Tn, const Vi &setup_done)
{    
    /*  == Variables decalration == */
    Vb child(_child);
    Vd s(_s);
    Vd t(_t);    

    GRBVar **S;
    GRBVar **C;
    GRBVar **W;
    GRBVar **Y;
    GRBVar** y_ij;
    GRBVar** Kx;

    int N = Sn + Tn;
    int J = Tn;
    int R = Sn;
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
        {
            string tag = "S_" + to_string(j);
            S[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, tag);
        }
    }

    C = new GRBVar*[M+1];
    for(int i = 1; i <= M; i++)
    {
        C[i] = new GRBVar[N+1];
        for(int j = 1; j <= N; j++)
        {
            string tag = "C_" + to_string(j);
            C[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, tag);
        }

    }

    W = new GRBVar*[M+1];
    for(int i = 1; i <= M; i++)
    {
        W[i] = new GRBVar[N+1];
        for(int j = 1; j <= N; j++)
        {
            string tag = to_string(j) + "Work_on" + to_string(i);
            W[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, tag);
        }
    }

    Y = new GRBVar*[N+1];
    for(int i = 1; i <= N; i++)
    {
        Y[i] = new GRBVar[N+1];
        for(int j = 1; j <= N; j++)
        {
            string tag = to_string(i) + "_precede_" + to_string(j);
            Y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, tag);
        }
    }

    // Auxiliary variables for lagrangian relaxation
    y_ij = new GRBVar*[N+1];    
    Kx = new GRBVar*[N+1];
    for(int i = 1; i <= N; i++)
    {            
        y_ij[i] = new GRBVar[N+1];      
        Kx[i] = new GRBVar[N+1];     
    }
    for(int i = 1; i <= N; i++)
    {                       
        for(int j = i+1; j <= N; j++)
        {
            // add auxiliary bool variables
            string aux_ij = to_string(i) + "_" + to_string(j);
            string aux_ji = to_string(j) + "_" + to_string(i);
            y_ij[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, aux_ij+"y_ij");
            y_ij[j][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, aux_ji+"y_ij");
            Kx[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, aux_ij+"_Kx");
            Kx[j][i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, aux_ij+"_Kx");
        }
    }
       
    /*  ================ Set objective =================
    ================================================  */    
    GRBQuadExpr objExpr = 0;
    for(int m = 1; m <= M; m++)
    {
        for(int j = 1; j <= J; j++)
        {
            // Add up all the completion time of test-jobs
            objExpr += C[m][j];
        }
    }

    // Lagrangian relaxation objective tweak
    for(int i = 1; i <= N; i++)
    {            
        for(int j = i+1; j <= N; j++)
        {
            objExpr += Kx[i][j] * lambda * (1 - y_ij[i][j]);
            objExpr += Kx[j][i] * lambda * (1 - y_ij[j][i]);
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
                    model.addConstr(Y[r][j] == 1, order_constr);
                }
            }
        }
    }

    // set-up job間 相依性    
    if(setup_done.size())
    {

        cout << "here!" << endl << endl;
        // last set up job 比餘下任一set up job早做
        vector<bool> not_done_setup(R+1, true);
        for(int i = 0; i < setup_done.size(); i++)  // 這裡0是對的
        {
            not_done_setup[setup_done[i]] = false;
        }
        for(int i = 1; i <= R; i++) // 這裡要1
        {
            if(not_done_setup[i])
            {
                string order_constr = "r"+ to_string(setup_done.back()) + "_r" + to_string(i)+ "_r_than_r";
                cout << order_constr << endl;
                model.addConstr(Y[setup_done.back() + J][i + J] == 1, order_constr);
            }
        }

        // set-up job間相依性，會連續影響
        if(setup_done.size() >= 2)
        {    for(int i = 0; i < setup_done.size()-1; i++)
            {
                // r_idx = setup_done[i] + J;
                string order_constr = "r"+ to_string(setup_done[i]) + "_r" + to_string(setup_done[i+1])+ "_r_than_r";
                cout << order_constr << endl;
                model.addConstr(Y[setup_done[i] + J][setup_done[i+1] + J] == 1, order_constr);
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

    // 試著讓排他性 relax
    /*================ Lagrangian Relaxation =================

        # Tx = S[m][j] - C[m][i] + My * ( 1 - Y[i][j] ) + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) + My(1 - y_ij[i][j])
            or S[m][i] - C[m][j] + My * Y[i][j] + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) + My(1 - y_ij[j][i])

        old constraint:
            Tx >= 0

        new constraints:
            Tx + M(1 - y) >= 0
            Kx - Tx * y >= 0
            Kx - Tx(y - 1) >= 0

        new obj:
            obj += Kx * lambda * (1 - y)

    ==========================================================  */  

    for(int m = 1; m <= M; m++)
    {
        for(int i = 1; i <= N; i++)
        {            
            for(int j = i+1; j <= N; j++)
            {                
                string aux_ij = to_string(i) + "_" + to_string(j);
                string aux_ji = to_string(j) + "_" + to_string(i);

                // decide whether the relaxed constraint is violated
                string ij_constr = "job" + to_string(i) + "_finishes_bf" + "_job" + to_string(j) + "_on_" + to_string(m);
                string ji_constr = "job" + to_string(j) + "_finishes_bf" + "_job" + to_string(i) + "_on_" + to_string(m);                
                model.addConstr(S[m][j] - C[m][i] + My * ( 1 - Y[i][j] ) + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) + My * ( 1 - y_ij[i][j] ) >= 0, ij_constr);
                model.addConstr(S[m][i] - C[m][j] + My * Y[i][j] + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) + My * ( 1 - y_ij[j][i] ) >= 0, ji_constr);
            
                // provide Kx as absolute subtitution for Tx
                string ij_abs1 = aux_ij + "_abs1";
                string ij_abs2 = aux_ij + "_abs2";
                string ji_abs1 = aux_ji + "_abs1";
                string ji_abs2 = aux_ji + "_abs1";
                model.addQConstr(Kx[i][j] - ( S[m][j] - C[m][i] + My * ( 1 - Y[i][j] ) + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) ) * y_ij[i][j] >= 0, ij_abs1);
                model.addQConstr(Kx[i][j] - ( S[m][j] - C[m][i] + My * ( 1 - Y[i][j] ) + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) ) * ( y_ij[i][j] - 1 ) >= 0, ij_abs2);
                model.addQConstr(Kx[j][i] - ( S[m][i] - C[m][j] + My * Y[i][j] + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) ) * y_ij[j][i] >= 0, ji_abs1);
                model.addQConstr(Kx[j][i] - ( S[m][i] - C[m][j] + My * Y[i][j] + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) ) * ( y_ij[j][i] - 1 ) >= 0, ji_abs2);
            }
        }
    } 

    // Optimize model
    model.optimize();    

    auto it = model.getVars();
    for(int i = 0; i < model.get(GRB_IntAttr_NumVars); i++)
    {
        cout << it[i].get(GRB_StringAttr_VarName) << ":" << it[i].get(GRB_DoubleAttr_X) << endl;
    }



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
    Vi setup_done = {};    

    cout << "optimal: " << relax_lb(child, s, t, Sn, Tn, setup_done) << endl;
}












// Deprecated Functions


// Lb_Compute_Param lbComputeParamSetUp(const Vb &_child, const Vd &_s, const Vd &_t, const Vi &setup_done)
// {
//     Vb child(_child);
//     Vd s(_s);
//     Vd t(_t);
//     for(int i = 0; i < setup_done.size(); ++i)
//     {        
//         child[setup_done[i]].reset();
//         s[setup_done[i]] = 0;
//     }
//     return Lb_Compute_Param(child, s, t);
// }


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