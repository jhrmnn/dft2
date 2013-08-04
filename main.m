function main(xyzfile, basisfile)
	atoms = read_xyz(xyzfile);
	N = sum([atoms.Z]);
    assert(mod(N, 2) == 0);
	basisdef = read_basis(basisfile);
    basis = build_basis(basisdef, atoms);
    [S, H, basis] = calc1ints(basis, atoms); % also normalizes contractions
    eri = calc2ints(basis);
    E = hf(S, H, eri, N);
    E = E+nuclear(atoms);
    fprintf('Energy: %.15f\n', E);
end