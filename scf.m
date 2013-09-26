function E = scf(S, H, eri, N, funcs, grid, basis)
    xc_funcs = rmfield(funcs, 'HF');
	X = ortho(S);
    diis_data = [];
    F = H;
    for i = 1:60
        P = solve_fock(F, X, N);
        if i > 1 && norm(P-P_old) < 1e-6
            break
        end
        J = tprod(eri, [1 2 -1 -2], P, [-1 -2]);
        F = H+J;
        if funcs.HF > 0
            K = -1/2*tprod(eri, [1 -1 2 -2], P, [-1 -2]);
            F = F+funcs.HF*K;
        end
        if ~isempty(fieldnames(xc_funcs))
            [Exc, Fxc] = eval_xc(P, xc_funcs, grid, basis);
            F = F+Fxc;
        end
        [F, diis_data] = diis(diis_data, P, S, X, F);
        P_old = P;
    end
    if isempty(fieldnames(xc_funcs))
        E = 1/2*tprod(P, [-1 -2], H+F, [-1 -2]);
    else
        E = 1/2*tprod(P, [-1 -2], H+F-Fxc, [-1 -2])+sum(Exc);
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