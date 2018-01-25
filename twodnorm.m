function [v] = twodnorm(A)
v = sqrt(A(:,1) .^ 2 + A(:,2) .^ 2);
end