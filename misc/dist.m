function R = dist(xyz1, xyz2)
    R = sqrt(sum(...
        bsxfun(@minus, permute(xyz1, [1 3 2]), permute(xyz2, [3 1 2])).^2, 3));
end
