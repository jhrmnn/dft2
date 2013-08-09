function main(xyzfile, basisfile, func)
	atoms = read_xyz(xyzfile);
	N = sum([atoms.Z]);
    assert(mod(N, 2) == 0);
	basisdef = read_basis(basisfile);
    basis = build_basis(basisdef, atoms);
    [S, H, basis] = calc1ints(basis, atoms); % also normalizes contractions
    if strcmpi(func, 'hf')
        dft = false;
    else
        grid = build_grid(atoms, 50, 302);
        dft = true;
    end
    eri = calc2ints(basis);
    if dft
        E = scf(S, H, eri, N, func, grid, basis);
    else
        E = scf(S, H, eri, N);
    end
    E = E+nuclear(atoms);
    fprintf('Energy: %.15f\n', E);
end