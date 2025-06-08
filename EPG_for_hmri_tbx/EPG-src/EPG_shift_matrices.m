function S = EPG_shift_matrices(Nmax)
% Generates shift matrices up to order Nmax for 3 standard EPG
% Layout of state vector is [F0 F0* Z0 F1 F-1 Z1 F2 F-2 Z2 ... etc]
% i.e. 3 states per n-value
%
%   Shaihan Malik 2017-07-20
% adapted by Luke J. Edwards to make indexing more explicit

N = (Nmax+1) * 3;
M = [3,Nmax+1];

% indices of non-zero elements
kidx = [];
sidx = [];

% F(k>=1)
kidx = [kidx, sub2ind(M, 1, 2:Nmax+1)]; % first index (F0+) omitted as treated explicitly below
sidx = [sidx, sub2ind(M, 1, 1:Nmax)];

% F(k<1)
kidx = [kidx, sub2ind(M, 2, 1:Nmax)]; % only to Nmax as most negative state relates to nothing
sidx = [sidx, sub2ind(M, 2, 2:Nmax+1)]; % negative states come from more negative states

% Z states don't shift
kall = sub2ind(M, 3, 1:Nmax+1);
kidx = [kidx, kall];
sidx = [sidx, kall];

% finally F0+ relates to F-1
kidx = [kidx, sub2ind(M, 1, 1)];
sidx = [sidx, sub2ind(M, 2, 2)];

% build matrix
S = sparse(kidx, sidx, 1, N, N);

end