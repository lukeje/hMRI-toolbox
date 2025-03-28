function hmri_corr_afi()

%% Input parameters
% Get sequence and tissue parameters
protocol = "ADPCA";
switch protocol
    case "ADPCA"
        FA      = [60, 60];        % Flip angles [deg]
        TR      = [20, 100];       % [ms]
        Phi0    = 36.0/3;            % [deg]
        B1range = (20:10:150)'/100; % convert such that 100% = 1
        dur1 = 11; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; %#ok<*NBRAK2> % [ms]
        Gamp{1} = [26,30,-30,30];           % [mT/m]
        dur2 = 55; % ms
        Gdur{2} = [3,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = [26,30,-30,30];           % [mT/m]
        
        % Get tissue parameters
        T1range = fliplr([1000, 1220, 1500, 3000]);     % [ms]
        T2range = 50;      % [ms]
        D       = 0.7;     % [µm^2/ms]
        
    case "PVPphantom"
        n = 3;
        FA      = [60, 60]; % Flip angles [deg]
        TR      = [1,n]*50; % [ms]
        Phi0    = 129.3;    % [deg]
        B1range = (50:5:150)'/100; % convert such that 100% = 1
        dur1 = 42; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,26,-26,26];           % [mT/m]
        dur2 = 42; % ms
        Gdur{2} = [n,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = [26,26,-26,26];           % [mT/m]
        
        % Get tissue parameters
        T1range = 1000;     % [ms]
        T2range = 196;      % [ms]
        D       = 0.6;     % [µm^2/ms]
end

% Max signal will be scaled to this value and rounded before calculation
% Set to inf to disable discretisation of signal
% Typically WM has largest signal values because of fast relaxation
maxS_WM = inf; 2^8;

%% Numerical simulations with EPG
% Build structure "diff" to account for diffusion effect
assert(length(Gamp)==length(Gdur))
for gIdx=1:length(Gamp)
    assert(length(Gdur{gIdx})==length(Gamp{gIdx}),'The vectors of gradient durations and amplitudes must have the same length!')
end
diff = struct('D', D*1e-9, 'G', Gamp, 'tau', Gdur); % struct assigns cell elements to separate struct array elements

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
            %phi_train   = RF_phase_cycle(npulse,Phi0);          % phases
            phi_train   = RF_phase_cycle_NehrkeSimplifiedUnitTR(npulse,Phi0,TR(2)/TR(1)); % phases
            
            % Calculate signals via EPG
            F0 = EPG_GRE_nTR(alpha_train, phi_train, TR, T1, T2, 'diff',diff, 'kmax',inf);
            S1(B1idx,T1idx,T2idx) = abs(F0(end-1));
            S2(B1idx,T1idx,T2idx) = abs(F0(end));
            
        end
    end
end


%% Simulate using exact result assuming perfect spoiling
S1e = abs(hmri_test_utils.dualTRernstd(B1range*FA(1),TR(1),TR(2),1./T1range));
S2e = abs(hmri_test_utils.dualTRernstd(B1range*FA(1),TR(2),TR(1),1./T1range));

%% Account for discretisation error during DICOM conversion
if isfinite(maxS_WM)
    sc = maxS_WM/max([S1(:);S2(:)]);
    S1 = round(S1*sc);
    S2 = round(S2*sc);

    S1e = round(S1e*sc);
    S2e = round(S2e*sc);
end

%% Calculate relative B1 map
B1app_grsp   = calc_AFI(S1, S2, TR(1),TR(2),FA(1));
B1app_compsp = calc_AFI(S1e,S2e,TR(1),TR(2),FA(1));

p = polyfit(100*mean(B1app_grsp,2),100*B1range,3);
disp(sprintf("%.7f ",p)) %#ok<DSPSP>
B1app_corr = polyval(p,100*B1app_grsp)/100;

figure
plot(100*B1range,100*(B1app_grsp-B1range),'-x')
hold on
plot(100*B1range,100*(B1app_compsp-B1range),'-o')
plot(100*B1range,100*(B1app_corr-B1range),'-s')
legend("T1 = "+T1range(:)+" ms"+[" grad only", " perfect", " corrected"],'Location',"Best")
xlabel("B1 (p.u.)")
ylabel("B1est - B1 (p.u.)")
hold off

end

function B1map = calc_AFI(Y1,Y2,TR1,TR2,nomFA)

% flip angle map in degrees
r=Y2./Y1;
n=TR2/TR1;
FAmap = acosd((r*n-1)./(n-r)); % Eq. (6) in Yarnykh, MRM (2007)

% relative B1 map
B1map = FAmap/nomFA;

end