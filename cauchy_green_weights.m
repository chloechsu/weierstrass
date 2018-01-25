function W = cauchy_green_weights(V,C)
  % Inputs:
  %  V  list of vertex positions
  %  C  list of control vertex positions
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  
  n = size(V, 1);
  m = size(C, 1);
  V = repmat(complex(V(:,1), V(:,2)), 1, m);
  C = repmat(transpose(complex(C(:,1), C(:,2))), n, 1);

  B = C - V;
  A = C - circshift(C, [0,1]);
  term1 = log(circshift(B,[0,-1]) ./ B);
  term1 = term1 .* circshift(B,[0,-1]) ./ circshift(A, [0,-1]);
  term2 = log(B ./ circshift(B,[0,1])) .* circshift(B,[0,1]) ./ A;
  W = (term1 - term2) ./ complex(0,2*pi);
  
end
