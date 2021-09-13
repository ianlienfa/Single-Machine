//
// Created by 林恩衍 on 2021/6/21.
//

#ifndef TP_SETANDTESTBB_H
#define TP_SETANDTESTBB_H
#include "Q_ELEMENT.h"
struct SetAndTestBB
{
    // data members
    int Sn, Tn, Jn;
    vector<B> parent;
    Vi s, t, j;

    // inline functions
    inline int tTojobidx(int i){return i + Sn;}
    inline int originidx(int i) { return (i >= Sn) ? i - Sn : i; }
    inline bool isT(int i){return i >= Sn;}

    SetAndTestBB(int Sn, int Tn) {
        this->Sn = Sn; this->Tn = Tn; this->Jn = Sn + Tn;
        parent.assign(Sn+Tn, B(0));
        j.resize(Sn+Tn);
        s.resize(Sn);
        t.resize(Tn);
    }
    int calSigmaCj(const Vi &seq, const Vi &j);
    void BSTsolve(Vi &ans, const Vi &j, const Vb &parent);
    void solve(const Vi &j, const Vb &parent);
    bool isRunnable(const B &visited, const B &parent);
    B toUnvisited(const B &visited);
    void printAns(const Vi &ans, const Vi &j);
    void printAns(const Vi &ans, const Vi &j, int sigmaCj);
//    int calSigmaWjCj(const Vi &seq, const Vi &j);
    Uheap toUheap(const Vi &seq)
//    void precedenceClearUp(const Vi &seq, vector<double> &prec, vector<double> &pi);

};

int SetAndTestBB::calSigmaCj(const Vi &seq, const Vi &j)
{
    int cum = 0;
    int ans = 0;
    for(int i = 0; i < seq.size(); i++)
    {
        cum += j[i];
        if(isT(seq[i])) ans += cum;
    }
    return ans;
}

void SetAndTestBB::printAns(const Vi &ans, const Vi &j)
{
    // print ans seq
    for(int i = 0; i < ans.size(); i++)
    {
        (isT(i)) ? printf("t") : printf("s");
        printf("%d", originidx(ans[i]));
        if(i != ans.size()-1) printf(" ,");
    }
    printf("\n");

    // print ans
    printf("Sigma Cj: %d\n", calSigmaCj(ans, j));
}

void SetAndTestBB::printAns(const Vi &ans, const Vi &j, int sigmaCj)
{
    // print ans seq
    for(int i = 0; i < ans.size(); i++)
    {
        (isT(i)) ? printf("t") : printf("s");
        printf("%d", originidx(ans[i]));
        if(i != ans.size()-1) printf(" ,");
    }
    printf("\n");

    // print ans
    printf("Sigma Cj: %d\n", sigmaCj);
}

void SetAndTestBB::solve(const Vi &j, const Vb &parent)
{
    // ans seq
    Vi ans;

    // BST solve ans
    BSTsolve(ans, j, parent);

    // printAns
    printAns(ans, j);
}

B SetAndTestBB::toUnvisited(const B &visited)
{
    B mask; mask.set();
    mask >>= (mask.size() - (Sn+Tn));
    B unvisited(visited);
    unvisited.flip();
    unvisited &= mask;
    return unvisited;
}

bool SetAndTestBB::isRunnable(const B &visited, const B &parent)
{
    B unvisited;
    unvisited = toUnvisited(visited);
    unvisited &= parent;
    return !unvisited.any();
}

//void precedenceClearUp(const Vi &seq, vector<double> &prec, vector<double> &pi)
//{
//    // iterate through seq and
//    for(int i = 0; i < seq.size(); i++)
//    {
//
//    }
//}



void SetAndTestBB::BSTsolve(Vi &ans, const Vi &j, const Vb &parent)
{
    ans.resize(Jn);

    // init state and q
    Qe q;
    q.push_front(E(B(0), Vi(), 0));
    int sz = 1, lv = 0;
    int best = 0x3FFFFFFF;
    Vi bestSeq;
    while(q.size())
    {
        // traverse every node on one level
        int next_sz = 0;
        while(sz--)
        {
            // get the state
            E &e = q.back();

            // calculate at the last lv
            if(lv == Jn)
            {
                // 直接把sigma Cj算出來
                int sigmaCj = calSigmaCj(e.seq, j);
                if(sigmaCj < best) {
                    best = sigmaCj;
                    bestSeq = e.seq;
                }
                printAns(e.seq, j, sigmaCj);
            }

            // retrieve its unvisited state
            B unvisited;
            unvisited = toUnvisited(e.visited);

            // start branching
            for(int i = 0; i < Jn; i++)
            {
                // prune if its parent is not yet visited
                if(unvisited.test(i) && isRunnable(e.visited, parent[i]))
                {
                    Vi v_in(e.seq);
                    v_in.push_back(j[i]);
                    q.push_front(E((B(e.visited).set(i)), v_in, e.c + j[i]));
                    next_sz++;
                }
            }

            // pop from queue
            q.pop_back();
        }
        sz = next_sz;
        lv++;
    }
    printf("Best Ans: ");
    printAns(bestSeq, j, best);
}



#endif //TP_SETANDTESTBB_H
