function [E] = expand_faces_to_edges(F)
% circshift the vertices in each face to have a row for each edge
F_shifted = [F; circshift(F,[0,1]); circshift(F,[0,2])];
F_shifted_flipped = [fliplr(F_shifted(:,1:2)) F_shifted(:,3)];
E = [F_shifted; F_shifted_flipped];
end