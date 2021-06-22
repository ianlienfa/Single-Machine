#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
#define SWAP(i, j, tmp) {tmp = i; i = j; j = tmp;}


// Heap implementation
// Max Heap with Update support
// arr 自己多開就好
// arr 的for loop reference之後 都是
// for(int i = 0; i < n; i++) {
//            cout << arr[i+1] << endl;
//        }
// bitset 的lowest bit視為1號位置
// 想得知bitset的內容使用印出比較安全

template <size_t N>
struct heapNode
{
    int qi;
    bitset<N> b;    // 01101 代表 job0, job2, job3 在此集合
    friend bool operator>(const heapNode<N> &nd1, const heapNode<N> &nd2)
    {
        return nd1.qi > nd2.qi;
    }
    friend ostream& operator<<(ostream& os, const heapNode<N>& node)
    {
        os << "[qi: " << node.qi << ", index: ";
        for(int i = 0; i < N; i++)
        {
            if((bitset<N>().set(i) & node.b).any())
                cout << i+1 << " ";
        }
        cout << "]" << endl;
        return os;
    }
    void sync(bitset<N> a)
    {
        b |= a;
    }
    heapNode(int qi, bitset<N> b)
    {
        this->qi = qi;
        this->b = b;
    }
};

// requires that n is fixed, n == N

template <size_t N>
struct MaxHeap{

    vector<heapNode<N>> arr;
    int n;
    int heapsize;
    inline int left(int i){return i*2;}
    inline int right(int i){return i*2+1;}
    inline int parent(int i){return i/2;}

    MaxHeap(int n, vector<double> qi)
    {
        arr.assign(n+1, heapNode<N>(0, bitset<N>(0)));
        this->n = n;
        heapsize = n;
        for(int i = 0; i < n; i++)
        {
            arr[i+1] = (heapNode<N>(qi[i], bitset<N>().set(i)));
        }
        for(int i = 0; i < n; i++) {
            cout << arr[i+1];
        }
        cout << endl;
        build_max_heap();
        for(int i = 0; i < n; i++) {
            cout << arr[i+1];
        }
        cout << endl;
    };

    void max_heapify(int i)
    {
        int l = left(i);
        int r = right(i);
        int largest;
        if(l <= heapsize && arr[l] > arr[i])
            largest = l;
        else
            largest = i;
        if(r <= heapsize && arr[r] > arr[largest])
            largest = r;
        if(largest != i)    // means that i is being changed
        {
            heapNode<N> temp(0, bitset<N>(0));
            SWAP(arr[largest], arr[i], temp);
            max_heapify(largest);
        }
    }

    heapNode<N> extract_max()
    {
        heapNode<N> max = arr[1];
        arr[1] = arr[heapsize];
        heapsize--;
        max_heapify(1);
        return max;
    }

    void update(heapNode<N> target, heapNode<N> newval)
    {

    }

    void build_max_heap()
    {
        for(int i = n/2; i >= 1; i--)
        {
            max_heapify(i);
        }
    };
};


int main()
{
    MaxHeap<5> heap(5, vector<double>({3, 5, 7, 89, 0}));
    for(int i = 0; i < 5; i++)
    {
        cout << heap.extract_max();
    }
}