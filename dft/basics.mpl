k[F] := proc(n) (3*pi^2*n)^(1/3) end;
S := proc(n, grad) grad/(2*k[F](n)*n) end;
r[s] := proc(n) (3/(4*pi*n))^(1/3) end;
k[s] := proc(n) sqrt(4*k[F](n)/pi) end;
T := proc(n, zeta, grad) grad/(2*phi(zeta)*k[s](n)*n) end;
phi := proc(zeta) ((1+zeta)^(2/3)+(1-zeta)^(2/3))/2 end;
P := proc(n, grad) S(n, grad)^2 end;
Z := proc(n, grad, kin) tW(n, grad)/kin end;
Alpha := proc(p, z) (5*p/3)*(1/z-1) end;
tW := proc(n, grad) 1/8*grad^2/n end;
Xi := proc(n, gradzeta) gradzeta/(2*(3*pi^2*n))^(1/3) end;
Gradzeta := proc(n, zeta, sigmauu, sigmadd, sigmaud)
    sqrt((1-zeta)^2*sigmauu+(1+zeta)^2*sigmadd-2*(1-zeta)*(1+zeta)*sigmaud)/n end;

