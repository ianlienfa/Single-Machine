
int main()
{
    SetAndTestSchedule st(3, 4);
    Vi seq = {1, 2, 3};
    B t1; t1.set(1);
    B t2; t2.set(1);
    B t3; t3.set(1).set(2);
    B t4; t4.set(1).set(2);
    B s1; s1.set(1).set(2).set(3).set(4);
    B s2; s2.set(3).set(4);
    B s3(0);
    Vb child = {B(0), s1, s2, s3};
    Vb prec = {B(0), t1, t2, t3, t4};
    Vd s = {0, 3, 2, 1};
    Vd t = {0, 5, 1, 2, 3};
    Uheap up = st.toUheap(seq, prec, child, s, t);
    list<int> v = up.find_seq();

}
