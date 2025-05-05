function hmri_corr_afi()

%% Input parameters
% Get sequence and tissue parameters
protocol = "KRK";

switch protocol
    case "ADPCA"
        FA      = [60, 60];        % Flip angles [deg]
        TR      = [100, 20];       % [ms]
        Phi0    = 36.0;            % [deg]
        B1range = (20:10:150)'/100; % convert such that 100% = 1

        dur1 = 55; % ms
        Gdur{1} = [3,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,30,-30,30];           % [mT/m]
        dur2 = 11; % ms
        Gdur{2} = [1,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = [26,30,-30,30];           % [mT/m]
        
        % Get tissue parameters
        T1range = 1200; %fliplr([1000, 1220, 1500, 3000]);     % [ms]
        T2range = 50;      % [ms]
        D       = 0.7;     % [µm^2/ms]

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0*TR1/TR2,TR1,TR2);
        
    case "KRK"
        FA      = [55,  55]; % Flip angles [deg]
        TR      = [25, 125]; % [ms]
        Phi0    = 36;        % [deg]
        B1range = (30:5:140)'/100; % convert such that 100% = 1
        dur1 = 7.2; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,26,-26,26];           % [mT/m]
        dur2 = 36;  % ms
        Gdur{2} = [3,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = [26,26,-26,26];           % [mT/m]
        
        % Get tissue parameters
        T1range = [1000,1200,3000]; % [ms]
        T2range = 50;      % [ms]
        D       = 0.7;     % [µm^2/ms]

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0*TR1/TR2,TR1,TR2);
        
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

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0*TR1/TR2,TR1,TR2);
end

%% Numerical simulations with EPG
% Build structure "diff" to account for diffusion effect
assert(length(Gamp)==length(Gdur))
for gIdx=1:length(Gamp)
    assert(length(Gdur{gIdx})==length(Gamp{gIdx}),'The vectors of gradient durations and amplitudes must have the same length!')
end
Gdiff = struct('D', D*1e-9, 'G', Gamp, 'tau', Gdur); % struct assigns cell elements to separate struct array elements

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
            phi_train   = phase_cycle(npulse,Phi0,TR(1),TR(2));            % phases
            
            % Calculate signals via EPG
            F0 = EPG_GRE_nTR(alpha_train, phi_train, TR, T1, T2, 'diff',Gdiff, 'kmax',inf);
            S1(B1idx,T1idx,T2idx) = abs(F0(end-1));
            S2(B1idx,T1idx,T2idx) = abs(F0(end));
            
        end
    end
end


%% Simulate using exact result assuming perfect spoiling
S1e = abs(hmri_test_utils.dualTRernstd(B1range*FA(1),TR(1),TR(2),1./T1range));
S2e = abs(hmri_test_utils.dualTRernstd(B1range*FA(1),TR(2),TR(1),1./T1range));

%% Calculate relative B1 map
B1app_grsp   = calc_AFI(S1, S2, TR(1),TR(2),FA(1));
B1app_compsp = calc_AFI(S1e,S2e,TR(1),TR(2),FA(1));

p = polyfit(100*mean(B1app_grsp,2),100*B1range,2);
disp(sprintf("%.7f ",p)) %#ok<DSPSP>
B1app_corr = polyval(p,100*B1app_grsp)/100;

figure

subplot(2,1,1)
plot(100*B1range,100*B1app_grsp,'-x')
hold on
plot(100*B1range,100*B1app_compsp,'-o')
plot(100*B1range,100*B1app_corr,'-s')
legend("T1 = "+T1range(:)+" ms"+[" grad only", " perfect", " corrected"],'Location',"Best")
xlabel("B1 (p.u.)")
ylabel("B1est (p.u.)")
hold off

subplot(2,1,2)
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