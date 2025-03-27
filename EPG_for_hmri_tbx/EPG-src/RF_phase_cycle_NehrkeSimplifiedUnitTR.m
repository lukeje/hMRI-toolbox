% Attempt at implementing AFI phase cycling scheme from
%   Nehrke, K. (2009), On the steady-state properties of actual flip angle imaging (AFI). 
%   Magn. Reson. Med., 61: 84-92. https://doi.org/10.1002/mrm.21592
%
% Simplified by choosing offset so that constant phase offsets in each TR are zero
% and assuming that TR2 = n*TR1 where n is a positive integer
function phi = RF_phase_cycle_NehrkeSimplifiedUnitTR(npulse,phi0,n)

phi0 = deg2rad(phi0);

% Different phase increments for even and odd TR
% Original paper assumes zero-based indexing, so here [1,3,5,...] are "even"
% and [2,4,6,...] are "odd"
dphi_oddeven = [1, n]*phi0*(n+1)/2; % pseudo-phase increment
increments = [0,0];
phi = zeros(npulse,1);
for k=2:npulse
    % both increments needs to end up scaled by current k, so update them both
    increments = wrapTo2Pi(increments + dphi_oddeven);

    oddevenidx = mod(k,2)+1;
    phi(k) = wrapTo2Pi(phi(k-1) + increments(oddevenidx));
end

end