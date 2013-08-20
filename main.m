function main(xyzfile, basisfile, funcs)
    if ~isfield(funcs, 'HF')
        funcs.HF = 0;
    end
    atoms = read_xyz(xyzfile);
    N = sum([atoms.Z]);
    assert(mod(N, 2) == 0);
    basisdef = read_basis(basisfile);
    basis = build_basis(basisdef, atoms);
    assert(~isempty(basis));
    [S, H, basis] = calc1ints(basis, atoms); % also normalizes contractions
    if length(fieldnames(funcs)) > 1 % if any xc functional
        grid = build_grid(atoms, 50, 300);
    else
        grid = [];
    end
    eri = calc2ints(basis);
    E_el = scf(S, H, eri, N, funcs, grid, basis);
    E = E_el+nuclear(atoms);
    fprintf('Energy: %.15f\n', E);
end
