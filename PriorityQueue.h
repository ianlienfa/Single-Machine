//
// Created by 林恩衍 on 2021/6/28.
//

#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
using namespace std;
struct Jr;
typedef deque<Jr> Qjr;
typedef deque<Qjr> QQjr;

template <class T>
class PriorityQueue
{
private:
    int sz;
    bool debug = false;
    vector<T> arr;
    bool (*cmpr)(const T& t1, const T &t2);
    void swap(T &t1, T &t2){T tmp = t1; t1 = t2; t2 = tmp;}
    void heapify(int i);
    inline int left(int i){return i * 2;}
    inline int right(int i){return i * 2 + 1;}
    inline int parent(int i){return i / 2;}
    void sizeInit(vector<T> v){
        arr.resize(v.size()+1);
        copy(v.begin(), v.end(), arr.begin()+1);    // make arr 0-based
        sz = v.size();
    }
//    void modeInit(bool MinMode)
//    {
//        if(MinMode == true)
//            cmpr = [](const T& t1, const T &t2){return t1 < t2;};
//        else
//            cmpr = [](const T& t1, const T &t2){return t1 > t2;};
//    }
public:
    // Default: Max heap
    PriorityQueue()
    {
        sz = 0;
        arr.clear();
    }
    PriorityQueue(bool (*cmpr)(const T& t1, const T &t2)) : PriorityQueue()
    {
        this->cmpr = cmpr;
    }
    PriorityQueue(vector<T> v, bool (*cmpr)(const T& t1, const T &t2) = [](const T& t1, const T &t2){return t1 > t2;}) : PriorityQueue()
    {
        sizeInit(v);
        this->cmpr = cmpr;
        construct_heap();
    }
//    PriorityQueue(bool MinMode)
//    {
//        PriorityQueue();
//        modeInit(MinMode);
//    }
//    PriorityQueue(vector<T> v, bool MinMode)
//    {
//        PriorityQueue();
//        modeInit(MinMode);
//        sizeInit(v);
//        construct_heap();
//    }
    void construct_heap();
    void push(T t);
    T extract();
    void bst_print();
    int size(){return sz;}
    T top(){return arr[1];}
    bool empty(){return (sz < 1);}
};

template<class T>
void PriorityQueue<T>::heapify(int i) {
    int l = left(i);
    int r = right(i);
    int largest = i;
    if(l <= sz && cmpr(arr[l], arr[i]))
        largest = l;
    if(r <= sz && cmpr(arr[r], arr[largest]))
        largest = r;
    if(i != largest) {
        swap(arr[largest], arr[i]);
        heapify(largest);
    }
}
// waiting for check
template<class T>
void PriorityQueue<T>::push(T t) {
    sz++;
    if(sz >= arr.size())
        arr.resize(sz+1);
    int cur = sz;
    arr[cur] = t;
    while(cur >= 1 && parent(cur) && cmpr(arr[cur], arr[parent(cur)])){
        int p = parent(cur);
        swap(arr[cur], arr[p]);
        cur = p;
    }

    if(debug)
    {
        cout << "arr: ";
        for(auto it: arr)
            cout << it << ", ";
        cout << endl;
    }
}

template<class T>
T PriorityQueue<T>::extract() {
    if(sz > 0) {
        T best = arr[1];
        swap(arr[1], arr[sz]);
        sz--;
        heapify(1);
        return best;
    }
    else {
    //    printf("error extracting.\n");
        return T();
    }
}

template<class T>
void PriorityQueue<T>::construct_heap() {
    for(int i = sz / 2; i >= 1; i--)
    {
        heapify(i);
    }
}

template<class T>
void PriorityQueue<T>::bst_print() {
    deque<int> q;
    if(sz < 1)
        return;
    q.push_front(1);
    int s = 1;
    while(q.size())
    {
        int next_s = 0;
        for(int ct = 0; ct < s; ct++)
        {
            int i = q.back();
            q.pop_back();
            if(left(i) <= sz) {
                q.push_front(left(i));
                next_s++;
            }
            if(right(i) <= sz) {
                q.push_front(right(i));
                next_s++;
            }
            cout << arr[i] << ", ";
        }
        s = next_s;
        cout << endl;
    }
}

// Jr means jobs with release time
struct Jr
{
    int r;
    int idx;
    int p;
    Jr(){}
    Jr(int idx, int p, int r = 0){this->r = r; this->idx = idx; this->p = p;}
    friend ostream& operator<<(ostream &out, const Jr &j){
        return out << "(idx: " << j.idx << ", p: " << j.p << ", r: " << j.r << ")";
    }
    friend bool operator<(const Jr &jr1, const Jr& jr2){return jr1.p < jr2.p;}
    friend bool operator>(const Jr &jr1, const Jr& jr2){return jr1.p > jr2.p;}
};

// qu is a 2d queue with release date, the job with same release date "should be save in the same row"
int SRPT(QQjr qu, bool debug)
{
    struct seqJ{
        int p;
        int ci;
        int idx;
        seqJ(int idx, int p, int ci){this->idx = idx; this->p = p; this->ci = ci;}
    };
    vector<seqJ> seq;

    // priority queue init
    bool (*cmpr)(const Jr &jr1, const Jr& jr2) = [](const Jr &jr1, const Jr& jr2){return jr1.p < jr2.p;};
    PriorityQueue<Jr> pq(cmpr);
    if(!qu.size()) return 0;
    int t = qu[0][0].r;
    int crtcl = qu[0][0].r;

    // iterations
    do
    {
        // push job at current time
        Qjr jobs;
        if(qu.size())
        {
            jobs = qu.front();  // get the earliest releasing jobs
            if (jobs.size()) {
                if (debug) cout << "add: " << endl;
                if (t == jobs[0].r) // if current time matches the release time, push these jobs into pq
                {
                    qu.pop_front();
                    for (Qjr::iterator it = jobs.begin(); it != jobs.end(); it++) {
                        pq.push(*it);    
                        if(debug) cout << *it << endl;                   
                    }
                }
            }
        }

        // priority queue inspection
        if (debug) 
        {
            cout << "pq snapshot: " << endl;
            pq.bst_print();
            cout << "extract: " << pq.top() << endl;
            cout << endl;
        }

        Jr leading;
        if(pq.size())
        {
            // pop the leading job from priority queue
            leading = pq.extract();
            // if(debug) cout << "leading: " << leading << endl;

            // update critical time point
            crtcl = (qu.size()) ? min(qu.front()[0].r, t + leading.p) : t + leading.p;

            // update vec and heap(if needed)
            int ci = crtcl;
            int pi = crtcl - t;
            int less = leading.p + t - crtcl;
            if(less)
                pq.push(Jr(leading.idx, less));
            else
                seq.push_back(seqJ(leading.idx, pi, ci));
        }
        else
        {
            // update critical time point
            if(qu.size())
                crtcl = qu.front()[0].r;
        }
        
        if(debug) cout << "crtrl: " << crtcl << endl;

        // update current time
        t = crtcl;

    }while(qu.size() || pq.size());

    // after calculation, compute sigma Cj
    int sigmaCj = 0;
    for(int i = 0; i < seq.size(); i++)
    {
        sigmaCj += seq[i].ci;

        // print for debug
        if(debug)
        printf("(idx: %d, Pi: %d, Ci: %d), ", seq[i].idx, seq[i].p, seq[i].ci);
    }
    if(debug) printf("\n");

    if(debug) cout << "sigmaCj: " << sigmaCj << endl;
    return sigmaCj;
}

int SRPT(QQjr &qu)
{
    return SRPT(qu, false);
}

#endif //PRIORITYQUEUE_H
