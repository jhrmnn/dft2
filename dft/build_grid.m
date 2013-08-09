function grid = build_grid(atoms, nrad, nang)
    nang = get_degree(nang);
    grid = make_grid_points(atoms, nrad, nang);
    w = eval_global_weights(grid, atoms);
    grid.w = grid.w.*w;
    grid = rmfield(grid, 'idx');
end

function grid = make_grid_points(atoms, nrad, nang)
    N = length(atoms);
    atoms = vertcat(atoms.R);
    grid_1 = make_one_grid(nang, nrad);
    ngrid = length(grid_1.w);
    xyz = bsxfun(@plus, grid_1.xyz, permute(atoms, [3 2 1]));
    w = repmat(grid_1.w, 1, N);
    idx = repmat(1:N, ngrid, 1);
    ngrid = ngrid*N;
    xyz = reshape(permute(xyz, [1 3 2]), ngrid, 3);
    w = reshape(w, ngrid, 1);
    idx = reshape(idx, ngrid, 1);
    grid = struct('xyz', xyz, 'w', w, 'idx', idx);
end

function grid = make_one_grid(nang, nrad)
    sphere = make_sphere(nang);
    radial = make_radial_grid(nrad);
    xyz = bsxfun(@times, sphere.xyz, permute(radial.r, [1 3 2]));
    w = sphere.w*radial.w;
    ngrid = nang*nrad;
    xyz = reshape(permute(xyz, [1 3 2]), ngrid, 3);
    w = reshape(w, ngrid, 1);
    grid = struct('xyz', xyz, 'w', w);
end

function grid = make_radial_grid(nrad)
    alpha = 5; % this value is for H-He, B-Ne, Ar
    m = 3; % this value is for molecules
    x = (1:nrad)/(nrad+1);
    r = -alpha*log(1-x.^m);
    rp = alpha*m*x.^(m-1)./(1-x.^m);
    w = r.^2.*rp/(nrad+1);
    grid = struct('r', r, 'w', w);
end

function grid = make_sphere(nang)
    sphere = getLebedevSphere(nang);
    grid = struct('xyz', [sphere.x sphere.y sphere.z], 'w', sphere.w);
end

function n =  get_degree(n)
    lebedev_degree = [6 14 26 38 50 74 86 110 146 170 194 230 266 302 ...
                      350 434 590 770 974 1202 1454 1730 2030 2354 2702 ...
                      3074 3470 3890 4334 4802 5294 5810];
	n = lebedev_degree(find(lebedev_degree>=n, 1));
end

function p = eval_global_weights(grid, atoms)
    ngrid = length(grid.w);
    N = length(atoms);
    atoms = vertcat(atoms.R);
    R = dist(atoms, atoms);
    r = dist(grid.xyz, atoms);
    mu = bsxfun(@rdivide, bsxfun(@minus, r, permute(r, [1 3 2])),...
                          permute(R, [3 1 2]));
    s = 1/2*(1-hhh(mu));
    s(logical(repmat(permute(eye(N), [3 1 2]), ngrid, 1))) = 1;
    w = prod(s, 3);
    p = w(sub2ind([ngrid N], (1:ngrid)', grid.idx))./sum(w, 2);
end

function g = hhh(mu)
    t1 = mu.^2.*(-3+mu.^2).^2;
    g = -mu.*(-3+mu.^2).*(-12+t1).*(-768+t1.*(-12+t1).^2)/8192;
end

