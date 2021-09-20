// qqjr means 2d queue of job-with-release-time, this function generates the sequence combining fixed and SRPT
QQjr SetAndTestSchedule::qqjrGen(E &e, int start)
{
    // for all the setup jobs not yet done, set their release time to "start"
    QQjr qqjr;
    map<int, vector<int>> mp;
    for(int i = 1; i <= Sn; i++)
    {
        if(!e.visited.test(i))
        {
            map<int, vector<int>>::iterator it = mp.find(start);
            if(it == mp.end())
            {
                mp.insert(make_pair(start, vector<int>({i})));
            }
            else
            {
                it->second.push_back(i);
            }
        }
    }

    for(int i = 1; i <= Tn; i++)
    {
        B parents_to_consider = (~e.visited) & prec[i];
        int r = start;
        for(int j = 1; j <= Sn; j++)
        {
            if(parents_to_consider.test(j))
                r += s[i];
        }

        map<int, vector<int>>::iterator it = mp.find(r);
        if(it == mp.end())
        {
            mp.insert(make_pair(r, vector<int>({tTojobidx(i)})));
        }
        else
        {
            it->second.push_back(tTojobidx(i));
        }
    }

    // turn back to QQjr
    for(auto it: mp)
    {
        deque<Jr> qjr;
        for(int i = 0; i < it.second.size(); i++)
        {
            // 用了新的j的pi值，可能會錯
            qjr.push_back(Jr(it.second[i], j[it.second[i]], it.first));
        }
        qqjr.push_back(qjr);
    }
    return qqjr;
}




pair<list<int>, double> SetAndTestSchedule::BFSCBBsolve(Vb prec, Vb child, Vd s, Vd t, bool print)
{
//init
this->s = s;
this->t = t;
this->prec = prec;
this->child = child;

// BST generate permutations, here we don't need c
vector<PriorityQueue<E>> contour;
// each level has one contour(priority queue)
contour.assign(Tn+Sn+1, PriorityQueue<E>([](const E &e1, const E &e2){return e1.lb < e2.lb;}));
// PriorityQueue<E> q([](const E &e1, const E &e2){return e1.lb < e2.lb;});
contour[0].push(E(B(0), Vi(), 0));
long long sz = 1, lv = 0;
double min = 0x3FFFFFFF;
list<int> min_seq;
int srptBest = 0x3FFFFFFF;
B levelVisited; levelVisited.set();     // levelvisited = 1111....1111, 用levelvisited來記錄每個contour是否為空
B mask(0); mask.flip(); mask >>= (mask.size() - (Tn + Sn+ 1));    // mask用來去掉不合理的level

while((levelVisited & mask).any())  // if there exists node unvisited...
{
for(int level = 0; level <= Tn+Sn; level++)
{
if(contour[level].size())
{
// cout << "round: " << sz++ << endl;
// q.bst_print();
// long long next_sz = 0;
E e = contour[level].top();
contour[level].extract();   // 會不會害有人沒branch到？

// for each permutation, solve it by Uheap
if(e.seq.size() == Sn)
{
if(debug) {
printf("\n\nnew complete seq starts.\n");
printf("set job seq: "); for(auto it: e.seq) printf("%d, ", it); printf("\n");
}
Uheap hp = toUheap(e.seq, prec, child, s, t);
// hp.debug = true;
list<int> v = hp.find_seq();
double ans = computeSeq(v);
if(ans < min)
{
if(debug) printf("ans: %.2f, min: %.2f\n", ans, min);
min = ans;
min_seq = v;
printf("now min: %.1f\n", min);
for(auto it: min_seq)
cout << it << ", " ; cout << endl << endl;
}
continue;
}

// branch -> wrong
// double fixed = 0;
// double sigmaCj_fixed = 0;
// for(int i = 0; i < e.seq.size(); i++)
// {
//     fixed += s[e.seq[i]];
//     sigmaCj_fixed += fixed;
// }

// compute out-tree generating left sequence
//         Vi outTreeVec(e.seq);
//         for(int i = 1; i <= Sn; i++)
//         {
//             if(!e.visited.test(i))
//             {
//                 outTreeVec.push_back(i);
//             }
//         }
//         cout << "outtree seq: ";
//         for(auto it: outTreeVec)
//             cout << it << " "; cout << endl;
//         Uheap hp = toUheap(outTreeVec, prec, child, s, t);
// //                hp.debug = true;
//         list<int> v = hp.find_seq();
//         double ans = computeSeq(v);
//         // cout << "ans: " << ans << endl;
//         // cout << "min: " << min << endl;
//         if(ans < min)
//         {
//             if(debug) printf("ans: %.2f, min: %.2f\n", ans, min);
//             min = ans;
//             min_seq = v;
//             for(auto it: min_seq)
//             cout << it << ", " ; cout << endl;
//             printf("now min: %.1f\n", min);
//         }

// branch
for(int i = 1; i <= Sn; i++)
{
// 現在node還沒去過的地方
// branching: move forward one more job
// compute the lower bound of the sequence and then push the node into the next level contour(priority queue)
if(!e.visited.test(i))  // if this set up job is not yet sequenced
{
Vi v_in(e.seq);
v_in.push_back(i);
E topush = E((B(e.visited).set(i)), v_in, 0);

// 1. tidy up the seq till this time point
// 2. compute lb using SRPT
QQjr job_with_rj = qqjrGen(topush);
// for(auto it: job_with_rj){
//     cout << "(";
//     for(auto itt: it)
//         cout << itt << ", ";
//     cout << ")";
//     }
// cout << endl;
int srpt = SRPT(job_with_rj);
topush.lb = srpt;
// cout << topush << " SigmaCjFixed: " << sigmaCj_fixed << " SRPT: " << srpt << endl;
contour[level+1].push(topush);
// if(srpt < srptBest)
// {
//     srptBest = srpt;

//     // // add computation
//     // for(int i = 1; i <= Sn; i++)
//     // {
//     //     if(!e.visited.test(i))
// }
// cout << " best srpt: " << srptBest << endl;
// next_sz++;
}
}
}
else
levelVisited.set(level, 0);
}
//        q.bst_print();
}

if(print)
{
printf("optimal seq:");
for(auto it: min_seq) printf("%d, ", it); printf("\n");
printf("min obj: %.1f\n", min);
}

return make_pair(min_seq, min);
}
