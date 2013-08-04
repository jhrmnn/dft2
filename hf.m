% Hartree-Fock implementation as described in Szabo & Ostlund

function E = hf(S, H, eri, N)
	X = ortho(S);
    eri_JK = eri-1/2*permute(eri, [1 4 3 2]);
    diis_data = [];
    F = H;
    for i = 1:60
        P = solve_fock(F, X, N);
        if i > 1 && norm(P-P_old) < 1e-6
            break
        end
        F = H+tprod(eri_JK, [1 2 -1 -2], P, [-1 -2]);
        [P, diis_data] = diis(diis_data, P, S, X, F);
        F = H+tprod(eri_JK, [1 2 -1 -2], P, [-1 -2]);
        P_old = P;
    end
    E = 1/2*tprod(P, [-1 -2], H+F, [-1 -2]);
end

function P = solve_fock(F, X, N)
	C = X*eigen(X'*F*X);
    P = 2*C(:, 1:N/2)*C(:, 1:N/2)';
end

function X = ortho(S)
	[U, S] = eigen(S);
	X = U*diag(1./sqrt(diag(S)));
end