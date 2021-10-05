#ifndef SETANDTESTSCHEDULE_H
#define SETANDTESTSCHEDULE_H

#include <iostream>
#include <vector>
#include <map>
#include <queue>
#include "Q_ELEMENT.h"
#include "PriorityQueue.h"
#include "Uheap.h"
using namespace std;
typedef struct Q_ELEMENT E;
typedef deque<Q_ELEMENT> Qe;
typedef vector<int> Vi;
typedef vector<double> Vd;
typedef vector<Vi> VVi;
typedef bitset<100> B;
typedef vector<B> Vb;

/* struct SetAndTestSchedule

/ Data structure
*   Basic structure for 1|bipartite|SigmaCj
*   <2d-Boolean vector> child, prec: information of the precendence dependecies 
*   <Double vector>     s, t: the processing time of set-up and test jobs

/ Usage
*   The constructor is not well-designed(I believe this has something to do with Bitset),
    so user should initiate all the attribute in their main function

*/

// p_si, p_ti should be 1-based
struct SetAndTestSchedule
{
private:
    int Sn, Tn, Jn;
    vector<B> child, prec;
    Vd s, t;
public:
    Vd j;
    bool debug = false;
    SetAndTestSchedule(int Sn, int Tn) {
        this->Sn = Sn; this->Tn = Tn;
        this->Jn = Sn+Tn;
        child.assign(Sn+1, B(0));
        prec.assign(Tn+1, B(0));
        s.resize(Sn+1);
        t.resize(Tn+1);
    }
    SetAndTestSchedule(){}
    void inputInit();
    void init(Vb prec, Vb child, Vd s, Vd t){this->prec = prec; this->child = child; this->s = s; this->t = t;}
    inline int stoParentidx(int i){return i;}
    inline int tTojobidx(int i){return i+Sn;}
    inline int originidx(int i) { return (i > Sn) ? i-(Sn) : i; }
    inline bool isT(int i){return i > Sn;}
    Uheap toUheap(const Vi &seq, const Vb &prec, const Vb &child, const Vd& s, const Vd &t);
    void precedenceClearUp(const Vi &seq, Vb prec, Vb child, Vi &parent);
    pair<list<int>, double> solve(Vb prec, Vb child, Vd s, Vd t, bool print = false);
    pair<list<int>, double> solve(){return solve(prec, child, s, t, true);};
    double computeSeq(list<int> v);
    pair<list<int>, double> DFSsolve(Vb prec, Vb child, Vd s, Vd t, bool print = false);
    void DFS(double &min, list<int> &min_seq, B visited, Vi seq, int lv);
    pair<list<int>, double> BFSsolve(Vb prec, Vb child, Vd s, Vd t, bool print);
    pair<list<int>, double> BFSBBsolve(Vb prec, Vb child, Vd s, Vd t, bool print);
    pair<list<int>, double> BFSBBsolve(){return BFSBBsolve(prec, child, s, t, true);};
    inline QQjr qqjrGen(E &e, bool debug = false){return qqjrGen_v0(e, debug);}
    QQjr qqjrGen_v0(E &e, bool debug);  // consider sequenced & non_sequenced set up job
    QQjr qqjrGen_v1(E &e, bool debug);  // take released test_job into consideration
    pair<list<int>, double> BFSCBBsolve(Vb prec, Vb child, Vd s, Vd t, bool print);
    pair<list<int>, double> localSearch(Vb prec, Vb child, Vd s, Vd t, bool print);
};

#endif