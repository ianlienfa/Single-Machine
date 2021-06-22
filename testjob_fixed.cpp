#include <iostream>
#include <vector>
using namespace std;



void fixed_test_job_solve(vector<int> testj, vector<vector<bool>> vl, int setj_ct)
{
    // testj is sequenced
    // vl are the setupjob sets of each testjob
    vector<string> seq; 
    for(int i = 0; i < testj.size(); i++)
    {
        for(int j = 0; j < setj_ct; j++)
        {
            if(vl[testj[i]-1][j])
            {
                seq.push_back("s"+to_string(j));
                for(int k = i+1; k < testj.size(); k++)
                {
                    vl[testj[k]-1]][j] = 0;                    
                }
            }
        }
        seq.push_back("t"+to_string(testj[i]));
    }
    for(auto it: seq)
    {
        cout << it << " ";
    }
}

int main()
{
    vector<vector<bool>> vl;    // the set up job set of each job
    vector<int> testj;    
    testj = {1, 2, 0};  // index of testj
    vl.assign(3, vector<bool>(5));
    vl[0] = {1, 0, 0, 0, 1};
    vl[1] = {0, 0, 1, 1, 0};
    vl[2] = {0, 1, 0, 0, 0};
}


