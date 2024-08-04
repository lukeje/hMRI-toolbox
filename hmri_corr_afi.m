function hmri_corr_afi()

%% Input parameters
% Get sequence parameters
FA      = [60, 60];        % Flip angles [deg]
TR      = [50, 150];       % [ms]
Phi0    = 137;             % [deg]
B1range = (30:10:130)'/100; % convert such that 100% = 1
Gdur{1} = [5];             %#ok<*NBRAK2> % [ms]
Gamp{1} = [40];            % [mT/m]
Gdur{2} = [50];            % [ms]
Gamp{2} = [40];            % [mT/m]

% Get tissue parameters
T1range = fliplr([500 700 1000 1500 3000]);     % [ms]
T2range = [50];      % [ms]
D       = [0.7];     % [um^2/ms]

% Max signal will be scaled to this value and rounded before calculation
% Set to inf to disable discretisation of signal
% Typically WM has largest signal values because of fast relaxation
maxS_WM = 2^8;

%% Numerical simulations with EPG
% Build structure "diff" to account for diffusion effect
assert(length(Gamp)==length(Gdur))
for gIdx=1:length(Gamp)
    assert(length(Gdur{gIdx})==length(Gamp{gIdx}),'The vectors of gradient durations and amplitudes must have the same length!')
    diff(gIdx) = struct('D', D*1e-9, 'G', Gamp{gIdx}, 'tau', Gdur{gIdx}); %#ok<AGROW>
end

assert(length(Gamp)==length(TR),'Each TR must have an associated set of gradients')
assert(FA(1)==FA(2),'AFI equation assumes both flip angles are equal')

% Run EPG simulation
nB1 = length(B1range);
nT1 = length(T1range);
nT2 = length(T2range);
S1  = zeros([nB1 nT1 nT2]);
S2  = zeros([nB1 nT1 nT2]);
for T1idx = 1:nT1 % loop over T1 values, can use parfor for speed
    
    T1 = T1range(T1idx);
    npulse = floor(15*T1/min(TR));   % ensure steady state signal
    npulse = npulse + mod(npulse,2); % ensure number of pulses even
    
    for T2idx = 1:nT2
        T2 = T2range(T2idx);
        
        for B1idx = 1:nB1  % loop over B1+ values
            B1eff = B1range(B1idx);
            
            % make train of flip angles and their phases
            alpha_train = repmat(deg2rad(FA*B1eff), 1, npulse/length(FA)); % flip angles
            phi_train   = RF_phase_cycle(npulse,Phi0);          % phases
            
            % Calculate signals via EPG
            F0 = EPG_GRE_nTR(alpha_train, phi_train, TR, T1, T2, 'diff', diff);
            S1(B1idx,T1idx,T2idx) = abs(F0(end-1));
            S2(B1idx,T1idx,T2idx) = abs(F0(end));
            
        end
    end
end

%% Account for discretisation error during DICOM conversion
if isfinite(maxS_WM)
    sc = maxS_WM/max([S1(:);S2(:)]);
    S1 = round(S1*sc);
    S2 = round(S2*sc);
end

%% Calculate relative B1 map
B1app = calc_AFI(S1,S2,TR(1),TR(2),FA(1));

plot(100*B1range,100*(B1app-B1range))
legend("T1 = "+T1range+" ms",'Location',"Best")
xlabel("B1 (p.u.)")
ylabel("B1est - B1 (p.u.)")

end

function B1map = calc_AFI(Y1,Y2,TR1,TR2,nomFA)

% flip angle map in degrees
r=Y2./Y1;
n=TR2/TR1;
FAmap = acosd((r*n-1)./(n-r)); % Eq. (6) in Yarnykh, MRM (2007)

% relative B1 map
B1map = FAmap/nomFA;

end