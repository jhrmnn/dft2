#include <cmath>
#include <iostream>
#include <ctime>
#include <libint2.h>
#include <boys.h>
#include <mex.h>

using namespace std;
using namespace libint2;

struct Pair_data_t {
    double zeta;
    double oo2z;
    double AB[3];
    double P[3];
    double PA[3];
    double K;
    double cacb;
};

struct Basis_func_t {
    double R[3];
    unsigned int l;
    double *zeta;
    double *d;
    unsigned int contr;
};

void print_clock(clock_t &start, clock_t &end, string name);

void prepare_pair_data(Pair_data_t &pair_data,
                       double ca, double zetaa, double A[3],
                       double cb, double zetab, double B[3]);

void prepare_data(Libint_eri_t &erieval, FmEval_Chebyshev3 &fmeval, int max_m,
                  Pair_data_t &pair_ab, Pair_data_t &pair_cd);

void calc2ints(double *eris, int nbasis, vector<Basis_func_t> &basis)
{
    libint2_static_init();
    int max_contr = 0;
    int max_l = 0;
    int nshell = basis.size();
    vector<int> idx(nshell);
    vector<int> idx_prim(nshell);
    for (int i = 0; i < nshell; i++) {
        if (basis[i].contr > max_contr)
            max_contr = basis[i].contr;
        if (basis[i].l > max_l)
            max_l = basis[i].l;
        idx[i] = (i == 0) ? 0 : idx[i-1]+(basis[i-1].l+1)*(basis[i-1].l+2)/2;
        idx_prim[i] = (i == 0) ? 0 : idx_prim[i-1]+basis[i-1].contr;
    }
    int nprims = idx_prim[nshell-1]+basis[nshell-1].contr;

    vector<Pair_data_t> pair_data(nprims*nprims);
    for (int s0 = 0; s0 < nshell; s0++) {
    for (int s1 = 0; s1 < nshell; s1++) {
        if (basis[s0].l < basis[s1].l)
            continue;
        for (int p0 = 0; p0 < basis[s0].contr; p0++) {
        for (int p1 = 0; p1 < basis[s1].contr; p1++) {
            prepare_pair_data(pair_data[(idx_prim[s0]+p0)+(idx_prim[s1]+p1)*nprims],
                              basis[s0].d[p0], basis[s0].zeta[p0], basis[s0].R,
                              basis[s1].d[p1], basis[s1].zeta[p1], basis[s1].R);
        }}
    }}

    int max_contr4 = max_contr*max_contr*max_contr*max_contr;
    vector<Libint_eri_t> erieval(max_contr4);
    libint2_init_eri(&erieval[0], max_l, 0);
    int max_m = 4*max_l;
    FmEval_Chebyshev3 fmeval(max_m);

    int p0123, max_shell_m;
    int ni, nj, nk, nl;
    clock_t timer = 0, tmp, zero = 0;
    for (int s0 = 0; s0 < nshell; s0++) {
    for (int s1 = 0; s1 < nshell; s1++) {
    for (int s2 = 0; s2 < nshell; s2++) {
    for (int s3 = 0; s3 < nshell; s3++) {
        if (basis[s0].l < basis[s1].l
            || basis[s2].l < basis[s3].l
            || basis[s2].l+basis[s3].l < basis[s0].l+basis[s1].l)
            continue;
        p0123 = 0;
        max_shell_m = basis[s0].l+basis[s1].l+basis[s2].l+basis[s3].l;
        tmp = clock();
        for (int p0 = 0; p0 < basis[s0].contr; p0++) {
        for (int p1 = 0; p1 < basis[s1].contr; p1++) {
        for (int p2 = 0; p2 < basis[s2].contr; p2++) {
        for (int p3 = 0; p3 < basis[s3].contr; p3++) {
            prepare_data(erieval[p0123], fmeval, max_m,
                         pair_data[(idx_prim[s0]+p0)+(idx_prim[s1]+p1)*nprims],
                         pair_data[(idx_prim[s2]+p2)+(idx_prim[s3]+p3)*nprims]);
            p0123++;
        }}}}
        timer += clock()-tmp;
        erieval[0].contrdepth = p0123;
        libint2_build_eri[basis[s0].l][basis[s1].l][basis[s2].l][basis[s3].l](&erieval[0]);
        int iub = idx[s0]+(basis[s0].l+1)*(basis[s0].l+2)/2;
        int jub = idx[s1]+(basis[s1].l+1)*(basis[s1].l+2)/2;
        int kub = idx[s2]+(basis[s2].l+1)*(basis[s2].l+2)/2;
        int lub = idx[s3]+(basis[s3].l+1)*(basis[s3].l+2)/2;
        double *eri = erieval[0].targets[0];
        for (int i = idx[s0]; i < iub; i++) {
        for (int j = idx[s1]; j < jub; j++) {
        for (int k = idx[s2]; k < kub; k++) {
        for (int l = idx[s3]; l < lub; l++) {
            eris[i+nbasis*(j+nbasis*(k+l*nbasis))] = *eri;
            eris[j+nbasis*(i+nbasis*(k+l*nbasis))] = *eri;
            eris[i+nbasis*(j+nbasis*(l+k*nbasis))] = *eri;
            eris[j+nbasis*(i+nbasis*(l+k*nbasis))] = *eri;
            eris[k+nbasis*(l+nbasis*(i+j*nbasis))] = *eri;
            eris[l+nbasis*(k+nbasis*(i+j*nbasis))] = *eri;
            eris[k+nbasis*(l+nbasis*(j+i*nbasis))] = *eri;
            eris[l+nbasis*(k+nbasis*(j+i*nbasis))] = *eri;
            eri++;
        }}}}
    }}}}
    print_clock(zero, timer, "timer");

    libint2_cleanup_eri(&erieval[0]);
    libint2_static_cleanup();
}

void prepare_pair_data(Pair_data_t &pair_data,
                       double ca, double zetaa, double A[3],
                       double cb, double zetab, double B[3])
{
    pair_data.zeta = zetaa+zetab;
    pair_data.oo2z = 0.5/pair_data.zeta;
    pair_data.AB[0] = A[0]-B[0];
    pair_data.AB[1] = A[1]-B[1];
    pair_data.AB[2] = A[2]-B[2];
    pair_data.P[0] = (zetaa*A[0]+zetab*B[0])/pair_data.zeta;
    pair_data.P[1] = (zetaa*A[1]+zetab*B[1])/pair_data.zeta;
    pair_data.P[2] = (zetaa*A[2]+zetab*B[2])/pair_data.zeta;
    pair_data.PA[0] = pair_data.P[0]-A[0];
    pair_data.PA[1] = pair_data.P[1]-A[1];
    pair_data.PA[2] = pair_data.P[2]-A[2];
    double AB2 = pair_data.AB[0]*pair_data.AB[0]+pair_data.AB[1]*pair_data.AB[1]
                 +pair_data.AB[2]*pair_data.AB[2];
    pair_data.K = exp(-zetaa*zetab*AB2/pair_data.zeta)/pair_data.zeta;
    pair_data.cacb = ca*cb;
}

void prepare_data(Libint_eri_t &erieval, FmEval_Chebyshev3 &fmeval, int max_m,
                  Pair_data_t &pair_ab, Pair_data_t &pair_cd)
{
    erieval.AB_x[0] = pair_ab.AB[0];
    erieval.AB_y[0] = pair_ab.AB[1];
    erieval.AB_z[0] = pair_ab.AB[2];
    erieval.PA_x[0] = pair_ab.PA[0];
    erieval.PA_y[0] = pair_ab.PA[1];
    erieval.PA_z[0] = pair_ab.PA[2];
    erieval.oo2z[0] = pair_ab.oo2z;
    erieval.CD_x[0] = pair_cd.AB[0];
    erieval.CD_y[0] = pair_cd.AB[1];
    erieval.CD_z[0] = pair_cd.AB[2];
    erieval.QC_x[0] = pair_cd.PA[0];
    erieval.QC_y[0] = pair_cd.PA[1];
    erieval.QC_z[0] = pair_cd.PA[2];
    erieval.oo2e[0] = pair_cd.oo2z;
    double PQ[3];
    PQ[0] = pair_ab.P[0]-pair_cd.P[0];
    PQ[1] = pair_ab.P[1]-pair_cd.P[1];
    PQ[2] = pair_ab.P[2]-pair_cd.P[2];
    double zetapeta = pair_ab.zeta+pair_cd.zeta;
    double W[3];
    W[0] = (pair_ab.zeta*pair_ab.P[0]+pair_cd.zeta*pair_cd.P[0])/zetapeta;
    W[1] = (pair_ab.zeta*pair_ab.P[1]+pair_cd.zeta*pair_cd.P[1])/zetapeta;
    W[2] = (pair_ab.zeta*pair_ab.P[2]+pair_cd.zeta*pair_cd.P[2])/zetapeta;
    double rho = pair_ab.zeta*pair_cd.zeta/zetapeta;
    erieval.WP_x[0] = W[0]-pair_ab.P[0];
    erieval.WP_y[0] = W[1]-pair_ab.P[1];
    erieval.WP_z[0] = W[2]-pair_ab.P[2];
    erieval.WQ_x[0] = W[0]-pair_cd.P[0];
    erieval.WQ_y[0] = W[1]-pair_cd.P[1];
    erieval.WQ_z[0] = W[2]-pair_cd.P[2];
    erieval.oo2ze[0] = 0.5/zetapeta;
    erieval.roz[0] = rho/pair_ab.zeta;
    erieval.roe[0] = rho/pair_cd.zeta;
    double PQ2 = PQ[0]*PQ[0]+PQ[1]*PQ[1]+PQ[2]*PQ[2];
    const double twopi25 = 34.98683665524972497; // 2*pi^2.5
    double prefactor = twopi25*pair_ab.K*pair_cd.K/sqrt(zetapeta)*pair_ab.cacb*pair_cd.cacb;
    double* F = new double[max_m+1];
    fmeval.eval(F, rho*PQ2, max_m);
    switch (max_m)
    {
        case 20: erieval.LIBINT_T_SS_EREP_SS(20)[0] = prefactor*F[20];
        case 19: erieval.LIBINT_T_SS_EREP_SS(19)[0] = prefactor*F[19];
        case 18: erieval.LIBINT_T_SS_EREP_SS(18)[0] = prefactor*F[18];
        case 17: erieval.LIBINT_T_SS_EREP_SS(17)[0] = prefactor*F[17];
        case 16: erieval.LIBINT_T_SS_EREP_SS(16)[0] = prefactor*F[16];
        case 15: erieval.LIBINT_T_SS_EREP_SS(15)[0] = prefactor*F[15];
        case 14: erieval.LIBINT_T_SS_EREP_SS(14)[0] = prefactor*F[14];
        case 13: erieval.LIBINT_T_SS_EREP_SS(13)[0] = prefactor*F[13];
        case 12: erieval.LIBINT_T_SS_EREP_SS(12)[0] = prefactor*F[12];
        case 11: erieval.LIBINT_T_SS_EREP_SS(11)[0] = prefactor*F[11];
        case 10: erieval.LIBINT_T_SS_EREP_SS(10)[0] = prefactor*F[10];
        case 9: erieval.LIBINT_T_SS_EREP_SS(9)[0] = prefactor*F[9];
        case 8: erieval.LIBINT_T_SS_EREP_SS(8)[0] = prefactor*F[8];
        case 7: erieval.LIBINT_T_SS_EREP_SS(7)[0] = prefactor*F[7];
        case 6: erieval.LIBINT_T_SS_EREP_SS(6)[0] = prefactor*F[6];
        case 5: erieval.LIBINT_T_SS_EREP_SS(5)[0] = prefactor*F[5];
        case 4: erieval.LIBINT_T_SS_EREP_SS(4)[0] = prefactor*F[4];
        case 3: erieval.LIBINT_T_SS_EREP_SS(3)[0] = prefactor*F[3];
        case 2: erieval.LIBINT_T_SS_EREP_SS(2)[0] = prefactor*F[2];
        case 1: erieval.LIBINT_T_SS_EREP_SS(1)[0] = prefactor*F[1];
        case 0: erieval.LIBINT_T_SS_EREP_SS(0)[0] = prefactor*F[0];
    }
    delete[] F;
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    int nshell = mxGetN(prhs[0]);
    vector<Basis_func_t> basis(nshell);
    int nbasis = 0;
    for (int i = 0; i < nshell; i++) {
        double *R = mxGetPr(mxGetField(prhs[0], i, "R"));
        copy(R, R+3, basis[i].R);
        int l = mxGetScalar(mxGetField(prhs[0], i, "l"));
        basis[i].l = l;
        basis[i].zeta = mxGetPr(mxGetField(prhs[0], i, "zeta"));
        basis[i].d = mxGetPr(mxGetField(prhs[0], i, "d"));
        basis[i].contr = mxGetM(mxGetField(prhs[0], i, "d"));
        nbasis += (l+1)*(l+2)/2;
    }
    mwSize dims[4] = {nbasis, nbasis, nbasis, nbasis};
    plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    double *eris = mxGetPr(plhs[0]);

    calc2ints(eris, nbasis, basis);
}

void print_clock(clock_t &start, clock_t &end, string name)
{
    cout << name << ": " << (double)(end-start)/CLOCKS_PER_SEC << endl;
}

