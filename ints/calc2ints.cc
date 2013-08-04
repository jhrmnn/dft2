#include <cmath>
#include <iostream>
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

void print_clock(clock_t &start, clock_t &end, string name);

void prepare_data(Libint_eri_t &erieval, FmEval_Chebyshev3 &fmeval, int max_m,
                  double ca, double zetaa, double A[3], 
                  double cb, double zetab, double B[3],
                  double cc, double zetac, double C[3], 
                  double cd, double zetad, double D[3]);

void calc2ints(double *eris, int nbasis, vector<Basis_func_t> &basis)
{
    libint2_static_init();
    int max_contr = 0;
    int max_l = 0;
    int nshell = basis.size();
    vector<int> idx(nshell);
    for (int i = 0; i < nshell; i++) {
        if (basis[i].contr > max_contr)
            max_contr = basis[i].contr;
        if (basis[i].l > max_l)
            max_l = basis[i].l;
        idx[i] = (i == 0) ? 0 : idx[i-1]+(basis[i-1].l+1)*(basis[i-1].l+2)/2;
    }
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
                         basis[s0].d[p0], basis[s0].zeta[p0], basis[s0].R,
                         basis[s1].d[p1], basis[s1].zeta[p1], basis[s1].R,
                         basis[s2].d[p2], basis[s2].zeta[p2], basis[s2].R,
                         basis[s3].d[p3], basis[s3].zeta[p3], basis[s3].R);
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

void prepare_data(Libint_eri_t &erieval, FmEval_Chebyshev3 &fmeval, int max_m,
                  double ca, double zetaa, double A[3], 
                  double cb, double zetab, double B[3],
                  double cc, double zetac, double C[3], 
                  double cd, double zetad, double D[3])
{
    double zeta = zetaa+zetab;
    double ABx = A[0]-B[0];
    double ABy = A[1]-B[1];
    double ABz = A[2]-B[2];
    double Px = (zetaa*A[0]+zetab*B[0])/zeta;
    double Py = (zetaa*A[1]+zetab*B[1])/zeta;
    double Pz = (zetaa*A[2]+zetab*B[2])/zeta;
    double PAx = Px-A[0];
    double PAy = Py-A[1];
    double PAz = Pz-A[2];
    erieval.AB_x[0] = ABx;
    erieval.AB_y[0] = ABy;
    erieval.AB_z[0] = ABz;
    erieval.PA_x[0] = PAx;
    erieval.PA_y[0] = PAy;
    erieval.PA_z[0] = PAz;
    erieval.oo2z[0] = 0.5/zeta;
    double eta = zetac+zetad;
    double CDx = C[0]-D[0];
    double CDy = C[1]-D[1];
    double CDz = C[2]-D[2];
    double Qx = (zetac*C[0]+zetad*D[0])/eta;
    double Qy = (zetac*C[1]+zetad*D[1])/eta;
    double Qz = (zetac*C[2]+zetad*D[2])/eta;
    double QCx = Qx-C[0];
    double QCy = Qy-C[1];
    double QCz = Qz-C[2];
    erieval.CD_x[0] = CDx;
    erieval.CD_y[0] = CDy;
    erieval.CD_z[0] = CDz;
    erieval.QC_x[0] = QCx;
    erieval.QC_y[0] = QCy;
    erieval.QC_z[0] = QCz;
    erieval.oo2e[0] = 0.5/eta;
    double PQx = Px-Qx;
    double PQy = Py-Qy;
    double PQz = Pz-Qz;
    double Wx = (zeta*Px+eta*Qx)/(zeta+eta);
    double Wy = (zeta*Py+eta*Qy)/(zeta+eta);
    double Wz = (zeta*Pz+eta*Qz)/(zeta+eta);
    double rho =zeta*eta/(zeta+eta);
    erieval.WP_x[0] = Wx-Px;
    erieval.WP_y[0] = Wy-Py;
    erieval.WP_z[0] = Wz-Pz;
    erieval.WQ_x[0] = Wx-Qx;
    erieval.WQ_y[0] = Wy-Qy;
    erieval.WQ_z[0] = Wz-Qz;
    erieval.oo2ze[0] = 0.5/(zeta+eta);
    erieval.roz[0] = rho/zeta;
    erieval.roe[0] = rho/eta;
    double AB2 = ABx*ABx+ABy*ABy+ABz*ABz;
    double CD2 = CDx*CDx+CDy*CDy+CDz*CDz;
    double PQ2 = PQx*PQx+PQy*PQy+PQz*PQz;
    double K1 = exp(-zetaa*zetab*AB2/zeta);
    double K2 = exp(-zetac*zetad*CD2/eta);
    double prefactor = 2*pow(M_PI, 2.5)*K1*K2/(zeta*eta*sqrt(zeta+eta))*ca*cb*cc*cd;
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

