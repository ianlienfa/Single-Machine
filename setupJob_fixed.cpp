#include <iostream>
#include <vector>
#include <limits>
#include <list>
using namespace std;
#define SWAP(i, j, tmp){tmp = i; i = j; j = tmp;}

struct NodeInfo
{
    double w, p;
    int p_init;
    NodeInfo(double w, double p, int p_init){this->w = w; this->p = p; this->p_init = p_init;}
    NodeInfo(){w = 0, p = -1, p_init = -1;};
};

struct HpNode
{
    double q;
    int idx;
};

struct UsetNode
{
    list<int> vi;
    int parent;
    int pos;
    bool connect;
};

struct Uheap
{
    vector<HpNode> arr;
    vector<UsetNode> st;
    vector<NodeInfo> info;
    int hpsize;

    Uheap(vector<NodeInfo> info);

    // Heap functions
    double extract_max();
    void max_heapify(int i);
    void update(double newq, int i);     // node i in heap has a new q value
    void upfloat(int i);                 // node i in heap has larger q then its parents, floats up
    void heap_node_swap(int i, int j);   // swaps two node in heap, should be used max_heapify or upfloat to maintain the heap properties

    // Set functions
    void merge(int i);      // node i in st is being chosen to merge with its parent
    int find_parent(int i); // set instructions, finds the representative of the set

    // find seq
    list<int> find_seq()
    {
        int ct = 0;
        int last_nd_p;
        while(hpsize > 0)
        {
            if(hpsize == 1)
                last_nd_p = find_parent(arr[0].idx);
            int maxq = extract_max();
            printf("round %d: maxq = %d\n", ct++, maxq);
        }

        // 2 ways to find the last node

        // method 1: find for the node with connect = 0
        for(int i = 1; i < st.size(); i++)
        {
            if(!st[i].connect)
                for(auto it: st[i].vi)
                    printf("%d, ", it);
            printf("\n");
        }

        // method 2: the last node
        for(auto it: st[last_nd_p].vi)
            printf("%d, ", it);
        printf("\n");

        return st[last_nd_p].vi;
    }
};


Uheap::Uheap(vector<NodeInfo> in) // info should be 1-based
{
    // initialize set, heap
    // there should be a dummy job infront of info, 除了idx 0是不能用的, idx 1 留給 w = -infinite, parent 為自己, pi = 0的job
    vector<NodeInfo> info = in;
    if(info.size() >= 1)
    {
        info[0] = (NodeInfo());
        info[1].w = -1.0e+20;
        info[1].p_init = 1;
    }
    else{printf("info size error!\n");}
    arr.assign(info.size(), HpNode());
    st.assign(info.size(), UsetNode());
    hpsize = info.size()-1;
    this->info = info;
    for(int i = 1; i < info.size(); i++)
    {
        arr[i].q = info[i].w/info[i].p;
        arr[i].idx = i;
        st[i].vi = list<int>(); st[i].vi.push_back(i);
        st[i].parent = info[i].p_init;
        st[i].connect = false;
        st[i].pos = i;
    }

    // build heap
    for(int i = hpsize/2; i > 0; i--)
    {
        max_heapify(i);
    }
}

void Uheap::heap_node_swap(int i, int j)
{
    int stnode_a = arr[i].idx;
    int stnode_b = arr[j].idx;
    int temp;
    HpNode tempnode;

    // swap pos->heap index
    SWAP(st[stnode_a].pos, st[stnode_b].pos, temp);

//    // swap heap->pos index
//    SWAP(arr[i].idx, arr[j].idx, temp);

    // swap heap element
    SWAP(arr[i], arr[j], tempnode);
}

void Uheap::max_heapify(int i)
{
    int l = i * 2;
    int r = i * 2 + 1;
    int largest;
    if(l <= hpsize && arr[l].q >  arr[i].q)
        largest = l;
    else
        largest = i;
    if(r <= hpsize && arr[r].q > arr[largest].q)
        largest = r;
    if(i != largest) {
        heap_node_swap(i, largest);
        max_heapify(largest);
    }
}

double Uheap::extract_max()
{
    if(hpsize < 1)
        cout << "extract_max: heap_size error." << endl;
    // debug
    for (auto it: arr)
        printf("%.2f ", it.q);
        printf(" ,heap size: %d", hpsize);
        printf("\n");
    
    HpNode maxnode = arr[1];
    double maxq = arr[1].q;
    heap_node_swap(1, hpsize);
    hpsize--;
    max_heapify(1);
    merge(maxnode.idx);
    return maxq;
}

int Uheap::find_parent(int i)
{
    return st[i].parent = (st[st[i].parent].connect) ? find_parent(st[i].parent) : st[i].parent;    // set: find_union
}

void Uheap::merge(int i)
{
    // union i and its parent
    int pt = find_parent(i);
    if(pt == i) return;

    // debug
    printf("concat: parent -- ");
    for(auto it: st[pt].vi)
        printf("%d ", it);
    printf(", child -- ");
    for(auto it: st[i].vi)
        printf("%d ", it);
    printf("\n");

    st[pt].vi.splice(st[pt].vi.end(), st[i].vi);    // concatenate two lists
    st[i].connect = true;
    st[i].pos = pt;

    // q recalculate
    info[pt].w += info[i].w;
    info[pt].p += info[i].p;
    double newq = info[pt].w/info[pt].p;

    if(st[pt].connect)
        printf("merge conflict: set[%d] as a representative has its connect on.\n", pt);

    // update heap
    update(newq, st[pt].pos);
}

void Uheap::update(double newq, int i)
{
    if(i <= 0 || i > hpsize)
    {
        printf("update: heap index problem, might be problem in merge.  insight: newq: %f, i: %d\n", newq, i);
    }
    if(newq == arr[i].q)
        return;
    if(newq > arr[i].q)
    {
        arr[i].q = newq;
        upfloat(i);
    }
    else
    {
        arr[i].q = newq;
        max_heapify(i);
    }
}

void Uheap::upfloat(int i)
{
    int p = i/2;
    if(p >= 1 && arr[i].q > arr[p].q)
    {
        heap_node_swap(i, p);
    }
}

int main()
{
    vector<NodeInfo> info;
    double w, p;
    int p_init;
    info.push_back(NodeInfo());
    // warning ! the parent of the first set up job should be its self
    // warning ! the weight of the first set up job should be -infinite
    while(cin >> w >> p >> p_init)
    {
        if(w == p && p == p_init && p_init == 0)
            break;
        // printf("%f %f %d\n", w, p, p_init);
        info.push_back(NodeInfo(w, p, p_init));
    }
    Uheap uheap(info);
    uheap.find_seq();
}

