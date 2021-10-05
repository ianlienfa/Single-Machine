//
// Created by 林恩衍 on 2021/6/21.
//

#ifndef Q_ELEMENT_H
#define Q_ELEMENT_H

#include <iostream>
#include <vector>
#include <bitset>
#include <limits>
#include <list>
using namespace std;
typedef vector<int> Vi;
typedef vector<Vi> VVi;
typedef bitset<100> B;
typedef vector<B> Vb;
typedef struct Q_ELEMENT E;

struct Q_ELEMENT
{
    B visited;
    Vi seq;
    double lb;
    Q_ELEMENT(){}
    Q_ELEMENT(B visited, Vi seq, double lb){this->visited = visited; this->seq = seq; this->lb = lb;}
    friend ostream& operator<<(ostream &out, const E &e)
    {
        out << "[(";
        for(auto it: e.seq)
            out << it << ", ";
        out << "), ";
        out << "lb:" << e.lb << "]";
        return out;
    }
};



#endif //TP_Q_ELEMENT_H
