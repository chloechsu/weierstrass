function [F_, V_, h_, branch_points, is_covering_disconnected] = double_branch_square_root(F, V, g)
% If there is a branch cut, returns a double branch cover F_branched
% and a holomorphic function h on the double branch cover such that h^2 = g
% Otherwise, if there is no branch cut, return F, V, and sqrt(g)

% we number the vertices on the double branch cover in the following way:
% 1 ... |V| are "upper vertices" with positive square root, 
% and |V|+1 ... 2|V| are "lower vertices" with negative square root

flag_no_branch_cut = true;
nV = size(V,1);
nF = size(F,1);
F_branched = [F; F + nV];
branch_points = [];
is_covering_disconnected = false;

% if a face in F is across the branch cut, modify the corresponding edge in
% F_Branched to cross the two layers
for a=1:nF
    for i=1:3
        j = mod(i,3)+1;
        k = mod(j,3)+1;
        gi = g(F(a,i)); gj = g(F(a,j)); gk = g(F(a,k));
        % If vertex F(a,i) is on the other side of the branch cut from 
        % vertex F(a,j) and F(a,k), change the face to cross the cut
        if is_across_branch_cut(gi, gj) && ...
                is_across_branch_cut(gi, gk)
            F_branched(a,i) = F(a,i) + nV;
            F_branched(a+nF,i) = F(a,i);
            flag_no_branch_cut = false;
        end   
        % If only one edge of the face crosses the branch cut (this happens
        % when g has a zero inside the face), the branched version of this
        % face would be a hexagon. 
        % In this case, we add a new vertex at the zero of g. This face is
        % represented as 6 new triangles with the new vertex p:
        % p vj vk', p vj' vk, p vi vj, p vi' vj', p vk vi, p vk' vi' 
        if is_across_branch_cut(gj, gk) && ...
                (~is_across_branch_cut(gi, gj)) && ...
                (~is_across_branch_cut(gi, gk))
            p = locate_zero(gi, gj, gk, ...
                V(F(a,i),:), V(F(a,j),:), V(F(a,k),:));
            branch_points = [branch_points; p];
            % the index of the newly added vertex
            p_ind = size(branch_points, 1) + 2 * nV;
            F_branched(a,:) = [p_ind F(a,j) F(a,k)+nV];
            F_branched(a+nF,:) = [p_ind F(a,j)+nV F(a,k)];
            F_branched = [F_branched; ...
                p_ind F(a,i) F(a,j);
                p_ind F(a,i)+nV F(a,j)+nV;
                p_ind F(a,k) F(a,i);
                p_ind F(a,k)+nV F(a,i)+nV];
            flag_no_branch_cut = false;
        end
    end
end

if flag_no_branch_cut
    F_ = F; V_ = V; h_ = sqrt(g);
else
    F_ = F_branched; V_ = [V; V; branch_points]; 
    h_ = [sqrt(g); -sqrt(g); zeros(size(branch_points,1))];
    if size(branch_points,1) == 0
        is_covering_disconnected = true;
    end
end
end

function [ret] = is_across_branch_cut(z1, z2)
% Returns whether there is a branch cut between two vertices with function
% value z1 and z2
ret = false;
if (sign(imag(z1)) ~= sign(imag(z2)))
    % x_intercept is where line segment between z1 and z2 crosses x axis
    x_intercept = (real(z1)*imag(z2) - real(z2)*imag(z1)) / imag(z2-z1);
    if x_intercept <= 0
        ret = true;
    end
end
end

function [p] = locate_zero(g1, g2, g3, v1, v2, v3)
% given the function value at three vertices v1, v2, v3, interpolate a
% piecewise linear complex-valued function and return the position of its
% zero

% represent the function as a 3 x 2 matrix 
% [Ar Ai; Br Bi; Cr Ci]
% such that g(x,y) = (Ar x + Br y + Cr) + i (Ai x + Bi y + Ci)
M = [v1 1; v2 1; v3 1] \ [real(g1) imag(g1); real(g2) imag(g2); real(g3) imag(g3)];
p = (-M(1:2,:)' \ M(3,:)')';
end
