% Attempt at implementing AFI phase cycling scheme from
%   Nehrke, K. (2009), On the steady-state properties of actual flip angle imaging (AFI). 
%   Magn. Reson. Med., 61: 84-92. https://doi.org/10.1002/mrm.21592
%
% Simplified by choosing offset so that constant phase offsets in each TR are zero
% This implementation intentionally contains an error in order to fit correction 
% factors for data collected with this scheme rather than the intended scheme
function phi = RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,TR1,TR2)

phi0 = deg2rad(phi0);

% initialise return vector
phi = zeros(npulse,1);

RFSpoilIncrement = 0;
RFSpoilPhase = 0;
for n=1:npulse
    if mod(n,2)
	    % First (shorter) TR, spoil less ala Nehrke
		RFSpoilIncrement = RFSpoilIncrement + phi0*TR1/TR2;
    else
	    RFSpoilIncrement = RFSpoilIncrement + phi0;
    end
	RFSpoilPhase = mod(RFSpoilPhase + RFSpoilIncrement, 20*pi);

	RFSpoilIncrement = mod(RFSpoilIncrement, 20*pi);

    phi(n) = RFSpoilPhase;
end

end