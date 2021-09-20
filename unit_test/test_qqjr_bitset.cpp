#include <iostream>
#include <bitset>
#include <vector>
using namespace std;

typedef bitset<30> B;

int main()
{
    B b;
    B visited(110);


    // extract the sequenced_setup_jobs from test job j's parents
    B sequenced_setup_jobs = (visited) & prec;

    // extract the parents in non_sequenced_setup_jobs 
    B non_sequenced_setup_jobs = (~sequenced_setup_jobs) & prec;
}