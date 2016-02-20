read "basics.mpl";
read "pbec.df";

epsilon[Ctpps] := proc(epsiloncrvepkzb, z)
    epsiloncrvepkzb*(1+d*epsiloncrvepkzb*z^3) end;
d := 2.8;

epsilon[Cpkzb] := proc(nu, nd, gradu, gradd, grad, Czetaxi, z)
    epsilon[Cpbe](nu+nd, (nu-nd)/(nu+nd), grad)
    *(1+Czetaxi*z^2)
    -(1+Czetaxi)*z^2
    *(nu/(nu+nd)*epsilon[Cpbe](nu, 1, gradu)
     +nd/(nu+nd)*epsilon[Cpbe](nd, 1, gradd)) end;
epsilon[CrevpkzbU] := proc(nu, nd, gradd, grad, Czetaxi, z)
    epsilon[Cpbe](nu+nd, (nu-nd)/(nu+nd), grad)
    *(1+Czetaxi*z^2-(1+Czetaxi)*z^2*nu/(nu+nd))
    -(1+Czetaxi)*z^2
    *nd/(nu+nd)*epsilon[Cpbe](nd, 1, gradd) end;
epsilon[CrevpkzbUD] := proc(n, zeta, grad, z)
    epsilon[Cpbe](n, zeta, grad)*(1-z^2) end;

CC := proc(zeta, xi)
    CC0(zeta)/(1+xi^2*((1+zeta)^(-4/3)+(1-zeta)^(-4/3))/2)^4 end;
CC0 := proc(zeta)
    0.53+0.87*zeta^2+0.50*zeta^4+2.26*zeta^6 end;
