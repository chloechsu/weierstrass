function [x] = solve_laplacian(L, b, forced_zeros)
L_nonsingular = L;
L_nonsingular(forced_zeros,:) = [];
L_nonsingular(:,forced_zeros) = [];
b(forced_zeros) = [];
x = -L_nonsingular \ b;
x = insertrows(x, zeros(size(forced_zeros)), forced_zeros-1);
%disp(strcat(['Laplacian condition number: ' num2str(condest(L_nonsingular))]));
end