//
// Created by 林恩衍 on 2021/10/6.
//

#ifndef RELAX_GUROBI_H
#define RELAX_GUROBI_H

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


struct Relax_Gurobi
{
    const static int M = 1;

    // Big Constant
    const static int Mw = 10;
    const static int My = 100000;

    // hyper param
    // 解放的constraint不連續，難relax
    const static int lambda = 50;

    static int relax_lb(const Vb &_child, const Vd &_s, const Vd &_t, const Vi &setup_done);
};



// 回傳objbound
int Relax_Gurobi::relax_lb(const Vb &_child, const Vd &_s, const Vd &_t, const Vi &setup_done)
{    
    /*  == Variables decalration == */
    Vb child(_child);
    Vd s(_s);
    Vd t(_t);    

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
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set("LogFile", "mip1.log");
    env.start();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    model.set("SolutionLimit", "1");

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

                    // If Relax
                    model.addConstr(Y[r][j] == 1, order_constr);
                }
            }
        }
    }

    // set-up job間 相依性    
    if(setup_done.size())
    {

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
                model.addConstr(Y[setup_done.back() + J][i + J] == 1, order_constr);
            }
        }

        // set-up job間相依性，會連續影響
        if(setup_done.size() >= 2)
        {    for(int i = 0; i < setup_done.size()-1; i++)
            {
                // r_idx = setup_done[i] + J;
                string order_constr = "r"+ to_string(setup_done[i]) + "_r" + to_string(setup_done[i+1])+ "_r_than_r";
                // cout << order_constr << endl;
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

    // 排他性, Exclusivity
    // 不是你先就是我先
    for(int m = 1; m <= M; m++)
    {
        for(int i = 1; i <= N; i++)
        {
            for(int j = i+1; j <= N; j++)
            {
                string ij_constr = "job" + to_string(i) + "_finishes_bf" + "_job" + to_string(j) + "_on_" + to_string(m);
                string ji_constr = "job" + to_string(j) + "_finishes_bf" + "_job" + to_string(i) + "_on_" + to_string(m);
                // cout << ji_constr << endl;
                model.addConstr(S[m][j] - C[m][i] + My * ( 1 - Y[i][j] ) + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) >= 0, ij_constr);
                model.addConstr(S[m][i] - C[m][j] + My * Y[i][j] + My * ( 1 - W[m][i] ) + My * ( 1 - W[m][j] ) >= 0);
            }
        }
    } 

    // Optimize model
    model.optimize();    
    return model.get(GRB_DoubleAttr_ObjBound);

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    // return lower bound
    return 0;
}


#endif //RELAX_GUROBI_H
