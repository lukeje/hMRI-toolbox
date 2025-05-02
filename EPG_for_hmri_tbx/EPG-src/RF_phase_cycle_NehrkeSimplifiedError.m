% Attempt at implementing AFI phase cycling scheme from
%   Nehrke, K. (2009), On the steady-state properties of actual flip angle imaging (AFI). 
%   Magn. Reson. Med., 61: 84-92. https://doi.org/10.1002/mrm.21592
%
% Simplified by choosing offset so that constant phase offsets in each TR are zero
% This implementation intentionally contains an error in order to fit correction 
% factors for data collected with this scheme rather than the intended scheme
function phi = RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,N1,N2)

phi0 = deg2rad(phi0);

% Different phase increments for even and odd TR
% Original paper assumes zero-based indexing, so here [1,3,5,...] are "even"
% and [2,4,6,...] are "odd"
scale_oddeven = [1, N2/N1];

increment = phi0*scale_oddeven(2);
phi = zeros(npulse,1);
phi(1) = increment;
for k=2:npulse
    increment = mod(increment + phi0*scale_oddeven(mod(k,2)+1), 20*pi);

    phi(k) = mod(phi(k-1) + increment, 20*pi);
end

end