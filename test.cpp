#include <iostream>
#include <bitset>
#include <vector>
using namespace std;

typedef bitset<30> B;

int main()
{
    vector<int> s = {1, 5, 6, 7, 8};
    int it;
    for(it = s.size()-1; it >= 0; --it)
    {
        if(s[it] == 2)
            break; 
    }
    for(it = it+1; it < s.size(); ++it)
    {
        cout << s[it] << endl;
    }

}