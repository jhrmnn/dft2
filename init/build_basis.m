% assign basis functions to atoms

function basis = build_basis(basisdef, atoms)
    idx = bsxfun(@eq, [atoms.Z]', [basisdef.Z]);
    n = sum(idx(:));
    basis = struct('R', cell(1, n), 'l', [], 'zeta', [], 'd', []);
    k = 1;
    for i = 1:length(atoms)
        for j = 1:length(basisdef)
            if ~idx(i, j)
                continue
            end
            b = basisdef(j);
            basis(k) = struct('R', atoms(i).R, 'l', b.l,...
                              'zeta', b.zeta, 'd', b.d);
            k = k + 1;
        end
    end
end