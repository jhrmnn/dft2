function E = scf(S, H, eri, N, func, grid, basis)
    if nargin == 7
        dft = true;
    else
        eri = eri-1/2*permute(eri, [1 4 3 2]);
        dft = false;
    end
	X = ortho(S);
    diis_data = [];
    F = H;
    for i = 1:60
        P = solve_fock(F, X, N);
        if i > 1 && norm(P-P_old) < 1e-8
            break
        end
        F = H+tprod(eri, [1 2 -1 -2], P, [-1 -2]);
        if dft
            [Exc, Fxc] = eval_xc(P, upper(func), grid, basis);
            F = F+Fxc;
        end
        [F, diis_data] = diis(diis_data, P, S, X, F);
        P_old = P;
    end
    if dft
        E = 1/2*tprod(P, [-1 -2], H+F-Fxc, [-1 -2])+Exc;
    else
        E = 1/2*tprod(P, [-1 -2], H+F, [-1 -2]);
    end
end

function P = solve_fock(F, X, N)
	C = X*eigen(X'*F*X);
    P = 2*C(:, 1:N/2)*C(:, 1:N/2)';
end

function X = ortho(S)
	[U, S] = eigen(S);
	X = U*diag(1./sqrt(diag(S)));
end