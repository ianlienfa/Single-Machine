
int main()
{
    int Tn = 2, Sn = 3;
    SetAndTestSchedule st(Sn, Tn);
    B t1; t1.set(2);
    B t2; t2.set(1).set(2).set(3);
    B s1; s1.set(2);
    B s2; s2.set(1).set(2);
    B s3; s3.set(2);
    Vi seq = {2, 1, 3}; Vb prec = {B(0), t1, t2}; Vb child = {B(0), s1, s2, s3}; Vi parent;
    st.precedenceClearUp(seq, prec, child, parent);
    for(int i = 1; i < parent.size(); i++)
        (st.isT(i)) ? printf("p[t%d]: %d ", st.originidx(i), parent[i]) : printf("p[s%d]: %d ", st.originidx(i), parent[i]);
    printf("\n");
}