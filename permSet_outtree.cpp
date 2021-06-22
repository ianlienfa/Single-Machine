#include <iostream>
#include <bitset>
#include <vector>
#include <queue>
#include "Uheap.h"
#include "Q_ELEMENT.h"
using namespace std;

typedef struct Q_ELEMENT E;
typedef vector<int> Vi;
typedef vector<Vi> VVi;
typedef bitset<100> B;
typedef vector<B> Vb;
typedef vector<double> Vd;
typedef deque<Q_ELEMENT> Qe;
#define SWAP(i, j){int tmp = i; i = j; j = tmp;}

typedef struct BFSElement{
        Vi seq; bool visited; 
        BFSElement(Vi seq, bool visited = false) {this->seq = seq; this->visited = visited;}
}Eb;

// p_si, p_ti should be 1-based
struct SetAndTestSchedule
{
private:
    int Sn, Tn, Jn;
    vector<B> child, prec;
    Vd s, t;
public:
    bool debug = false;
    SetAndTestSchedule(int Sn, int Tn) {
        this->Sn = Sn; this->Tn = Tn;
        this->Jn = Sn+Tn;
        child.assign(Sn+1, B(0));
        prec.assign(Tn+1, B(0));
        s.resize(Sn+1);
        t.resize(Tn+1);
    }
    void init(Vb prec, Vb child, Vd s, Vd t){this->prec = prec; this->child = child; this->s = s; this->t = t;}
    inline int stoParentidx(int i){return i;}
    inline int tTojobidx(int i){return i+Sn;}
    inline int originidx(int i) { return (i > Sn) ? i-(Sn) : i; }
    inline bool isT(int i){return i > Sn;}   
    Uheap toUheap(const Vi &seq, const Vb &prec, const Vb &child, const Vd& s, const Vd &t);
    void precedenceClearUp(const Vi &seq, Vb prec, Vb child, Vi &parent);
    pair<list<int>, double> solve(Vb prec, Vb child, Vd s, Vd t, bool print = false);
    double computeSeq(list<int> v);
    pair<list<int>, double> DFSsolve(Vb prec, Vb child, Vd s, Vd t, bool print = false);
    void DFS(double &min, list<int> &min_seq, B visited, Vi seq, int lv);
    pair<list<int>, double> BFSsolve(Vb prec, Vb child, Vd s, Vd t, bool print);
};

void SetAndTestSchedule::precedenceClearUp(const Vi &seq, Vb prec, Vb child, Vi &parent)
{
    // seq: 存有set up job idx, 0-based
    // prec: 大小為Tn+1, prec[test job idx]為存有map到 set up job的 bitset array，兩個都是1-based
    // parent: 大小為 Tn+Sn+1, idx 要經過轉換
    // parse seq and set up the precedencies of set up jobs
    parent.resize(Jn+1);
    for(int i = seq.size() - 1; i >= 0; i--) {
        if (i == 0)
            parent[seq[i]] = seq[i];
        else
            parent[seq[i]] = seq[i - 1];

        // for all test job of set job i clean precedence if test job can be process later
        for (int j = 1; j <= Tn; j++) {
            if (child[seq[i]].test(j)) {
                parent[tTojobidx(j)] = seq[i];
                for (int p = 1; p <= Sn; p++) {
                    if (p == seq[i])
                        continue;
                    if (prec[j].test(p)) {
                        child[p].reset(j);
                        prec[j].reset(p);
                    }
                }
            }
        }
    }
}

Uheap SetAndTestSchedule::toUheap(const Vi &seq, const Vb &prec, const Vb &child, const Vd& s, const Vd &t)
{
    // info init and some checking
    vector<NodeInfo> info;
    info.resize(Tn+Sn+1); // 0 is forbidden
    info[0] = NodeInfo();
    if(child.size() != Sn+1 || prec.size() != Tn+1 || seq.size() != Sn || s.size() != Sn+1 || t.size() != Tn+1)
        printf("child or prec size error! child.size: %d, prec: %d, seq: %d\n", child.size(), prec.size(), seq.size());

    // use seq & prec & child to construct new 1-1 parent dependencies
    vector<int> parent;

    // precedence clear up
    precedenceClearUp(seq, prec, child, parent);

    // import set up job first
    for(int i = 1; i <= Sn; i++)
    {
        if(seq[0] == i)
            info[i] = (NodeInfo(-1.0e+10, s[i], parent[i]));
        else
            info[i] = (NodeInfo(0, s[i], parent[i]));
    }
    for(int j = 1; j <= Tn; j++)
    {
        info[tTojobidx(j)] = (NodeInfo(1, t[j], parent[tTojobidx(j)]));
    }

    if(info.size() != Sn + Tn + 1)
        printf("toUheap: info size error!\n");

    return Uheap(info);
}

double SetAndTestSchedule::computeSeq(list<int> seq)
{
    if(seq.size() != Sn + Tn)
        printf("ComputeSeq: seq size error!\n");
    int cum = 0;
    int ans = 0;
    if(debug) printf("whole seq: ");
    for(auto i: seq)
    {
        // if(debug) printf("%d, ", i);
        // i is the number in seq, which are job indecies
        cum += (isT(i)) ? t[originidx(i)] : s[i];
        ans += (isT(i)) ? cum : 0;
        if(debug) printf("i: %d, cum: %d, ans: %d\n", i, cum, ans);
    }
    if(debug) printf("\n");
    return ans;
}


pair<list<int>, double> SetAndTestSchedule::BFSsolve(Vb prec, Vb child, Vd s, Vd t, bool print)
{
    //init
    this->s = s;
    this->t = t;
    this->child = child;
    this->prec = prec;

    deque<BFSElement> q; 
    Vi v;
    for(int i = 1; i <= Sn; i++)
        v.push_back(i);
    q.push_front(Eb(v));

    double min = 0x3FFFFFFF;
    list<int> min_seq;

    int sz = 1;
    int l = 0;
    while(q.size())
    {        
        int next_sz = 0;
        while(sz--)
        {
            Eb &now = q.back();            
            if(!now.visited)
            {
                // run outtree 
                // for(auto it: now.seq) printf("%d, ", it); printf("\n");
                Uheap hp = toUheap(now.seq, prec, child, s, t);
//                hp.debug = true;
                list<int> v = hp.find_seq();
                double ans = computeSeq(v);
                if(ans < min)
                {
                    if(debug) printf("ans: %.2f, min: %.2f\n", ans, min);
                    min = ans;
                    min_seq = v;
                    for(auto it: min_seq) printf("%d, ", it); printf("\n");
                    printf("now min: %.1f\n", min);
                }
            }

            for(int i = l; i < now.seq.size(); i++)
            {
                Eb e(now);
                SWAP(now.seq[l], now.seq[i]);
                if(i == l)                
                    e.visited = true;
                q.push_front(e);
                next_sz++;
            }

            q.pop_back();
        }
        sz = next_sz;
        l++;
    }
    return make_pair(min_seq, min);
}


void SetAndTestSchedule::DFS(double &min, list<int> &min_seq, B visited, Vi seq, int lv)
{
    if(lv == Sn)
    {
        Uheap hp = toUheap(seq, prec, child, s, t);
        list<int> v = hp.find_seq();
        double ans = computeSeq(v);
        if(ans < min)
        {
            if(debug) printf("ans: %.2f, min: %.2f\n", ans, min);
            min = ans;
            min_seq = v;
            printf("now min: %.1f\n", min);
            printf("seq: "); for(auto it: min_seq) printf("%d, ", it); printf("\n");
        }
        return;
    }

    // for(int i = 1; i <= Sn; i++)
    // {
    //     if(!visited.test(i))
    //     {
    //         Vi v(seq);
    //         B b(visited);
    //         b.set(i);
    //         v.push_back(i);
    //         DFS(min, min_seq, b, v, lv+1);
    //     }
    // }

    for(int i = Sn; i >= 1; i--)
    {
        if(!visited.test(i))
        {
            Vi v(seq);
            B b(visited);
            b.set(i);
            v.push_back(i);
            DFS(min, min_seq, b, v, lv+1);
        }
    }
}

pair<list<int>, double> SetAndTestSchedule::DFSsolve(Vb prec, Vb child, Vd s, Vd t, bool print)
{
    //init
    this->s = s;
    this->t = t;
    this->child = child;
    this->prec = prec;

    double min = 0x3FFFFFFF;
    list<int> min_seq;
    DFS(min, min_seq, B(0), Vi(0), 0);
    return make_pair(min_seq, min);
}

pair<list<int>, double> SetAndTestSchedule::solve(Vb prec, Vb child, Vd s, Vd t, bool print)
{
    //init
    this->s = s;
    this->t = t;

    // BST generate permutations, here we don't need c
    Qe q;
    q.push_front(E(B(0), Vi(), 0));
    int sz = 1, lv = 0;
    double min = 0x3FFFFFFF;
    list<int> min_seq;

    while(q.size())
    {
        int next_sz = 0;
        while(sz--)
        {
            E &e = q.back();

            // for each permutation, solve it by Uheap
            if(lv == Sn)
            {

                if(debug) {
                    printf("\n\nnew case starts.\n");
                    printf("set job seq: "); for(auto it: e.seq) printf("%d, ", it); printf("\n");
                }
                Uheap hp = toUheap(e.seq, prec, child, s, t);
//                hp.debug = true;
                list<int> v = hp.find_seq();
                double ans = computeSeq(v);
                if(ans < min)
                {
                    if(debug) printf("ans: %.2f, min: %.2f\n", ans, min);
                    min = ans;
                    min_seq = v;
                    printf("now min: %.1f\n", min);
                }
            }

            // branch
            for(int i = 1; i <= Sn; i++)
            {
                if(!e.visited.test(i))
                {
                    Vi v_in(e.seq);
                    v_in.push_back(i);
                    q.push_front(E((B(e.visited).set(i)), v_in, 0));
                    next_sz++;
                }
            }
            q.pop_back();
        }
        sz = next_sz;
        lv++;
    }

    if(print)
    {
        printf("optimal seq:");
        for(auto it: min_seq) printf("%d, ", it); printf("\n");
        printf("min obj: %.1f\n", min);
    }

    return make_pair(min_seq, min);
}


int main()
{
    // all 1-based, init and ref
    int tmp, Sn, Tn;;
    cin >> Sn >> Tn;
    SetAndTestSchedule st(Sn, Tn);
    Vd s ,t;
    Vb prec;
    Vb child;
    child.assign(Sn+1, B(0));
    prec.assign(Tn+1, B(0));
    s.resize(Sn+1);
    t.resize(Tn+1);

    // construct
    for(int i = 1; i <= Sn; i++){cin >> s[i]; }
    for(int i = 1; i <= Tn; i++){cin >> t[i]; }
    cin >> tmp;
    if(tmp != -1) cout << "input error" << endl;

    // input precedencies
    // i means set up job i is test job j's parent
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
    
    // prec 1 ~ Tn, child 1~Sn 對應的bitset
    // vector<d> t, s has the pi of jobs
    // st.debug = true;
    pair<list<int>, double> pp = st.DFSsolve(prec, child, s, t, true);
    for(auto it : pp.first) printf("%d, ", it); printf("\n");
    // st.init(prec, child, s, t);
    // list<int> v = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 34, 40, 38, 37, 55, 20, 21, 44, 22, 25, 41, 29, 43, 46, 24, 52, 26, 33, 30, 51, 60, 56, 27, 59, 53, 23, 57, 50, 28, 42, 58, 35, 47, 54, 31, 32, 36, 48, 45, 49, 39};
    // cout << st.computeSeq(v) << endl;

}
