#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include <libint2.h>
#include <boys.h>
#include <mex.h>

using namespace std;
using namespace libint2;

struct Basis_func_t {
    double R[3];
    unsigned int l;
    double *zeta;
    double *d;
    unsigned int contr;
};

struct Pair_data_t {
    double zeta;
    double oo2z;
    double AB[3];
    double P[3];
    double PA[3];
    double K;
    double cacb;
};

struct Ssss_data_t {
    double zetapeta;
    double rho;
    double PQ2;
    double prefactor;
};

void print_clock(clock_t start, clock_t end, string name);

void prepare_pair(Pair_data_t &pair,
                  double ca, double zetaa, double A[3],
                  double cb, double zetab, double B[3]);

void prepare_ssss(Ssss_data_t &out,
                  Pair_data_t &pair_ab, Pair_data_t &pair_cd);

void prepare_erival(Libint_eri_t &erieval, FmEval_Chebyshev3<double> &fmeval, int max_m,
                    double *F, Ssss_data_t &ssss,
                    Pair_data_t &ab, Pair_data_t &cd);

int l_to_n(int l)
{
    return (l+1)*(l+2)/2;
}

void calc2ints(double *eris, int nbasis, vector<Basis_func_t> &basis)
{
    libint2_static_init();
    int nshell = basis.size();
    vector<int> ibas(nshell);
    vector<int> iprim(nshell);
    int max_contr = 0, max_l = 0;
    for (int i = 0; i < nshell; i++) {
        if (basis[i].contr > max_contr)
            max_contr = basis[i].contr;
        if (basis[i].l > max_l)
            max_l = basis[i].l;
        ibas[i] = (i == 0) ? 0 : ibas[i-1]+l_to_n(basis[i-1].l);
        iprim[i] = (i == 0) ? 0 : iprim[i-1]+basis[i-1].contr;
    }
    int nprims = iprim[nshell-1]+basis[nshell-1].contr;

    int max_contr4 = max_contr*max_contr*max_contr*max_contr;
    vector<Libint_eri_t> erieval(max_contr4);
    libint2_init_eri(&erieval[0], max_l, 0);
    int max_m = 4*max_l;
    FmEval_Chebyshev3<double> fmeval(max_m);

    vector<Pair_data_t> pairs(nprims*nprims);
    for (int s0 = 0; s0 < nshell; s0++) {
    for (int s1 = 0; s1 < nshell; s1++) {
        for (int p0 = 0; p0 < basis[s0].contr; p0++) {
        for (int p1 = 0; p1 < basis[s1].contr; p1++) {
            prepare_pair(pairs[(iprim[s0]+p0)+(iprim[s1]+p1)*nprims],
                         basis[s0].d[p0], basis[s0].zeta[p0], basis[s0].R,
                         basis[s1].d[p1], basis[s1].zeta[p1], basis[s1].R);
        }}
    }}

    int p0123, max_shell_m;
    int s0l, s1l, s2l, s3l;
    double F;
    double *F_p = new double[max_m+1];
    double *eri, eri_v;
    Ssss_data_t ssss;
    for (int s0 = 0; s0 < nshell; s0++) {
    for (int s1 = 0; s1 <= s0; s1++) {
    for (int s2 = 0; s2 <= s0; s2++) {
    for (int s3 = 0; s3 <= s2; s3++) {
        s0l = s0; s1l = s1; s2l = s2; s3l = s3;
        if (basis[s0l].l < basis[s1l].l)
            swap(s0l, s1l);
        if (basis[s2l].l < basis[s3l].l)
            swap(s2l, s3l);
        if (basis[s2l].l+basis[s3l].l < basis[s0l].l+basis[s1l].l) {
            swap(s0l, s2l);
            swap(s1l, s3l);
        }
        max_shell_m = basis[s0l].l+basis[s1l].l+basis[s2l].l+basis[s3l].l;
        if (max_shell_m == 0) {
            eri_v = 0;
            for (int p0 = 0; p0 < basis[s0l].contr; p0++) {
            for (int p1 = 0; p1 < basis[s1l].contr; p1++) {
            for (int p2 = 0; p2 < basis[s2l].contr; p2++) {
            for (int p3 = 0; p3 < basis[s3l].contr; p3++) {
                prepare_ssss(ssss,
                             pairs[(iprim[s0l]+p0)+(iprim[s1l]+p1)*nprims],
                             pairs[(iprim[s2l]+p2)+(iprim[s3l]+p3)*nprims]);
                fmeval.eval(&F, ssss.rho*ssss.PQ2, 0);
                eri_v += ssss.prefactor*F;
            }}}}
            eri = &eri_v;
        }
        else {
            p0123 = 0;
            for (int p0 = 0; p0 < basis[s0l].contr; p0++) {
            for (int p1 = 0; p1 < basis[s1l].contr; p1++) {
            for (int p2 = 0; p2 < basis[s2l].contr; p2++) {
            for (int p3 = 0; p3 < basis[s3l].contr; p3++) {
                prepare_erival(erieval[p0123], fmeval, max_m, F_p, ssss,
                               pairs[(iprim[s0l]+p0)+(iprim[s1l]+p1)*nprims],
                               pairs[(iprim[s2l]+p2)+(iprim[s3l]+p3)*nprims]);
                p0123++;
            }}}}
            erieval[0].contrdepth = p0123;
            libint2_build_eri[basis[s0l].l][basis[s1l].l][basis[s2l].l][basis[s3l].l](&erieval[0]);
            eri = erieval[0].targets[0];
        }
        int iub = ibas[s0l]+l_to_n(basis[s0l].l);
        int jub = ibas[s1l]+l_to_n(basis[s1l].l);
        int kub = ibas[s2l]+l_to_n(basis[s2l].l);
        int lub = ibas[s3l]+l_to_n(basis[s3l].l);
        for (int i = ibas[s0l]; i < iub; i++) {
        for (int j = ibas[s1l]; j < jub; j++) {
        for (int k = ibas[s2l]; k < kub; k++) {
        for (int l = ibas[s3l]; l < lub; l++) {
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
    delete[] F_p;

    libint2_cleanup_eri(&erieval[0]);
    libint2_static_cleanup();
}

void prepare_pair(Pair_data_t &pair,
                       double ca, double zetaa, double A[3],
                       double cb, double zetab, double B[3])
{
    pair.zeta = zetaa+zetab;
    pair.oo2z = 0.5/pair.zeta;
    pair.AB[0] = A[0]-B[0];
    pair.AB[1] = A[1]-B[1];
    pair.AB[2] = A[2]-B[2];
    pair.P[0] = (zetaa*A[0]+zetab*B[0])/pair.zeta;
    pair.P[1] = (zetaa*A[1]+zetab*B[1])/pair.zeta;
    pair.P[2] = (zetaa*A[2]+zetab*B[2])/pair.zeta;
    pair.PA[0] = pair.P[0]-A[0];
    pair.PA[1] = pair.P[1]-A[1];
    pair.PA[2] = pair.P[2]-A[2];
    double AB2 = pair.AB[0]*pair.AB[0]
                 +pair.AB[1]*pair.AB[1]
                 +pair.AB[2]*pair.AB[2];
    pair.K = exp(-zetaa*zetab*AB2/pair.zeta)/pair.zeta;
    pair.cacb = ca*cb;
}

void prepare_ssss(Ssss_data_t &out,
                  Pair_data_t &ab, Pair_data_t &cd)
{
    out.zetapeta = ab.zeta+cd.zeta;
    out.rho = ab.zeta*cd.zeta/out.zetapeta;
    double PQ[3];
    PQ[0] = ab.P[0]-cd.P[0];
    PQ[1] = ab.P[1]-cd.P[1];
    PQ[2] = ab.P[2]-cd.P[2];
    out.PQ2 = PQ[0]*PQ[0]+PQ[1]*PQ[1]+PQ[2]*PQ[2];
    const double twopi25 = 34.98683665524972497; // 2*pi^2.5
    out.prefactor = twopi25*ab.K*cd.K/sqrt(out.zetapeta)*ab.cacb*cd.cacb;
}

void prepare_erival(Libint_eri_t &erieval, FmEval_Chebyshev3<double> &fmeval, int max_m,
                    double *F, Ssss_data_t &ssss,
                    Pair_data_t &ab, Pair_data_t &cd)
{
    prepare_ssss(ssss, ab, cd);
    fmeval.eval(F, ssss.rho*ssss.PQ2, max_m);
    switch (max_m)
    {
        case 16: erieval.LIBINT_T_SS_EREP_SS(16)[0] = ssss.prefactor*F[16];
        case 15: erieval.LIBINT_T_SS_EREP_SS(15)[0] = ssss.prefactor*F[15];
        case 14: erieval.LIBINT_T_SS_EREP_SS(14)[0] = ssss.prefactor*F[14];
        case 13: erieval.LIBINT_T_SS_EREP_SS(13)[0] = ssss.prefactor*F[13];
        case 12: erieval.LIBINT_T_SS_EREP_SS(12)[0] = ssss.prefactor*F[12];
        case 11: erieval.LIBINT_T_SS_EREP_SS(11)[0] = ssss.prefactor*F[11];
        case 10: erieval.LIBINT_T_SS_EREP_SS(10)[0] = ssss.prefactor*F[10];
        case  9: erieval.LIBINT_T_SS_EREP_SS( 9)[0] = ssss.prefactor*F[ 9];
        case  8: erieval.LIBINT_T_SS_EREP_SS( 8)[0] = ssss.prefactor*F[ 8];
        case  7: erieval.LIBINT_T_SS_EREP_SS( 7)[0] = ssss.prefactor*F[ 7];
        case  6: erieval.LIBINT_T_SS_EREP_SS( 6)[0] = ssss.prefactor*F[ 6];
        case  5: erieval.LIBINT_T_SS_EREP_SS( 5)[0] = ssss.prefactor*F[ 5];
        case  4: erieval.LIBINT_T_SS_EREP_SS( 4)[0] = ssss.prefactor*F[ 4];
        case  3: erieval.LIBINT_T_SS_EREP_SS( 3)[0] = ssss.prefactor*F[ 3];
        case  2: erieval.LIBINT_T_SS_EREP_SS( 2)[0] = ssss.prefactor*F[ 2];
        case  1: erieval.LIBINT_T_SS_EREP_SS( 1)[0] = ssss.prefactor*F[ 1];
        case  0: erieval.LIBINT_T_SS_EREP_SS( 0)[0] = ssss.prefactor*F[ 0];
    }
    erieval.AB_x[0] = ab.AB[0];
    erieval.AB_y[0] = ab.AB[1];
    erieval.AB_z[0] = ab.AB[2];
    erieval.PA_x[0] = ab.PA[0];
    erieval.PA_y[0] = ab.PA[1];
    erieval.PA_z[0] = ab.PA[2];
    erieval.oo2z[0] = ab.oo2z;
    erieval.CD_x[0] = cd.AB[0];
    erieval.CD_y[0] = cd.AB[1];
    erieval.CD_z[0] = cd.AB[2];
    erieval.QC_x[0] = cd.PA[0];
    erieval.QC_y[0] = cd.PA[1];
    erieval.QC_z[0] = cd.PA[2];
    erieval.oo2e[0] = cd.oo2z;
    double W[3];
    W[0] = (ab.zeta*ab.P[0]+cd.zeta*cd.P[0])/ssss.zetapeta;
    W[1] = (ab.zeta*ab.P[1]+cd.zeta*cd.P[1])/ssss.zetapeta;
    W[2] = (ab.zeta*ab.P[2]+cd.zeta*cd.P[2])/ssss.zetapeta;
    erieval.WP_x[0] = W[0]-ab.P[0];
    erieval.WP_y[0] = W[1]-ab.P[1];
    erieval.WP_z[0] = W[2]-ab.P[2];
    erieval.WQ_x[0] = W[0]-cd.P[0];
    erieval.WQ_y[0] = W[1]-cd.P[1];
    erieval.WQ_z[0] = W[2]-cd.P[2];
    erieval.oo2ze[0] = 0.5/ssss.zetapeta;
    erieval.roz[0] = ssss.rho/ab.zeta;
    erieval.roe[0] = ssss.rho/cd.zeta;
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

void print_clock(clock_t start, clock_t end, string name)
{
    cout << name << ": " << (double)(end-start)/CLOCKS_PER_SEC << endl;
}

