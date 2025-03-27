% Attempt at implementing AFI phase cycling scheme from
%   Nehrke, K. (2009), On the steady-state properties of actual flip angle imaging (AFI). 
%   Magn. Reson. Med., 61: 84-92. https://doi.org/10.1002/mrm.21592
function phi = RF_phase_cycle_Nehrke(npulse,phi0,N1,N2,offset)

phi0 = deg2rad(phi0);

if ~exist('offset','var')
    offset = 0;
end

% Different phase increments for even and odd TR
% Original paper assumes zero-based indexing, so here [1,3,5,...] are "even"
% and [2,4,6,...] are "odd"
phi_n = phi0*(N1+N2)/2;
scale_oddeven = [N1, N2];
const_oddeven = scale_oddeven.*(phi0*(1-N2)/2 + offset);

increments = [0,0];
phi = zeros(npulse,1);
for k=2:npulse
    % both increments needs to end up scaled by current k, so update them both
    increments = wrapTo2Pi(increments + phi_n*scale_oddeven);

    oddeven = mod(k,2)+1;
    phi(k) = wrapTo2Pi(phi(k-1) + increments(oddeven) + const_oddeven(oddeven));
end

end