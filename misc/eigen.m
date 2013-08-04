% gives orthonormal eigenvectors and eigenvalues in ascending order

function [V, D] = eigen(A)
	[V, D] = eig((A+A')/2);
end
