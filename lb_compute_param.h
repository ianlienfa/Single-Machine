#ifndef LBCOMPUTEPARAM_H
#define LBCOMPUTEPARAM_H

#include <iostream>
#include <vector>
#include <bitset>
#include "Q_ELEMENT.h"
using namespace std;

typedef vector<double> Vd;
typedef vector<B> Vb;

struct Lb_Compute_Param
{
    Vb child;
    Vd s;
    Vd t;
    Lb_Compute_Param(){}
    Lb_Compute_Param(Vb _child, Vd _s, Vd _t){child = _child; s = _s; t = _t;}
};

#endif //QTEST_LBCOMPUTEPARAM_H
