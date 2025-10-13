function E = E_diff(E,diff,kmax,N,dk,negk)
% E = E_diff(E,diff,kmax,N,negk)
%
%    Function to build E operator with diffusion effects for standard EPG
%    3 states per k-value). 
%   
%       E = relaxation matrix (diag(E2 E2 E1))
%       diff = structure with fields:
%              G    - Gradient amplitude(s)
%              tau  - Gradient durations(s)
%              D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%       Attention: Periods of G = 0 must also be listed! In other words, sum(tau) = TR.
%       This is because diffusion happens even if gradients are switched off.
%
% Shaihan Malik July 2017

if ~exist('dk','var')
    dk = [];
end

[bDL, bDT] = EPG_diffusion_weights(diff.G,diff.tau,diff.D,0:kmax,dk);

% E is a simple diagonal matrix - just need to compute this diagonal
Ed = diag(E);
Ed = reshape(Ed,3,[]);
EdT = Ed(1:2,:);
EdL = Ed(3  ,:);

EdT = full(EdT).*reshape(bDT,1,1,[]);
EdL = full(EdL).*reshape(bDL,1,1,[]);

% also include negative k
if exist('negk','var') && negk
    EdT = cat(3,EdT(:,:,end:-1:2),EdT);
    EdL = cat(3,EdL(:,:,end:-1:2),EdL);
end

% Combine them
Ed = cat(1,EdT,EdL);
Ed=Ed(:);

%%% Now use sparse diagonal function to define overall matrix
E = spdiags(Ed,0,N,N);

end
