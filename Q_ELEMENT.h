//
// Created by 林恩衍 on 2021/6/21.
//

#ifndef Q_ELEMENT_H
#define Q_ELEMENT_H

typedef struct Q_ELEMENT E;
typedef vector<int> Vi;
typedef vector<Vi> VVi;
typedef bitset<100> B;
typedef vector<B> Vb;
typedef deque<Q_ELEMENT> Qe;

struct Q_ELEMENT
{
    B visited;
    Vi seq;
    int c;
    Q_ELEMENT(){}
    Q_ELEMENT(B visited, Vi seq, int c){this->visited = visited; this->seq = seq; this->c = c;}
};



#endif //TP_Q_ELEMENT_H
