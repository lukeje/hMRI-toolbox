% Attempt at implementing AFI phase cycling scheme from
%   Nehrke, K. (2009), On the steady-state properties of actual flip angle imaging (AFI). 
%   Magn. Reson. Med., 61: 84-92. https://doi.org/10.1002/mrm.21592
%
% Simplified by choosing offset so that constant phase offsets in each TR are zero
% This implementation intentionally contains an error in order to fit correction 
% factors for data collected with this scheme rather than the intended scheme
function phi = RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,scale_firsttr_secondtr)

phi0 = deg2rad(phi0);

% initialise return vector
phi = zeros(npulse,1);

% different phase increments for every other TR
increment = mod(phi0*scale_firsttr_secondtr(1), 20*pi);
phi(1) = increment;
for k=2:npulse
    increment = mod(increment + phi0*scale_firsttr_secondtr(mod(k-1,2)+1), 20*pi);

    phi(k) = mod(phi(k-1) + increment, 20*pi);
end

end