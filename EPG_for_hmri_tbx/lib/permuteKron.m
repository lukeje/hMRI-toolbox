% Provides matrix P such that
%   P*kron(A,B)=kron(B,A) 
% where A and B are vectors.
% Uses Equation 20 of
%   https://doi.org/10.1080/03081088108817379
%
% luke.edwards@ucl.ac.uk

function P = permuteKron(dimA,dimB)

P=sparse([],[],[],dimA*dimB,dimA*dimB,dimA*dimB);

idA=speye(dimA);
idB=speye(dimB);
for vecNum=1:dimA
    tempVec=idA(:,vecNum);
    P=P+kron(tempVec.',kron(idB,tempVec));
end

end
    