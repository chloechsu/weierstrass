function [L] = compute_laplacian(F, V, E, cotan)
% the discrete Laplacian operator according to cotan formula
sz = [size(V,1) size(V,1)];
ind1 = sub2ind(sz, E(:,1), E(:,1));
ind2 = sub2ind(sz, E(:,1), E(:,2));
% set issparse as true for accumarray
L = accumarray([ind1; ind2], [0.5 * cotan; -0.5 * cotan], [prod(sz) 1], [], [], true);
L = reshape(L, sz);
L = sparse(L);
end