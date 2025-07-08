function hmri_afi_reperror()

%% Input parameters
% Get sequence and tissue parameters
protocol = "ADPCA";
nreps=100;
switch protocol
    case "Lutti"
        FA      = [60, 60];        % Flip angles [deg]
        TR      = [100, 20];       % [ms]
        phi0    = 36.0;            % [deg]
        B1range = (30:10:130)'/100; % convert such that 100% = 1

        Gdur{1} = 55; % [ms]
        Gamp{1} = 26; % [mT/m]
        Gdur{2} = 11; % [ms]
        Gamp{2} = 26; % [mT/m]

        % Get tissue parameters
        [T1range,T2range,D] = tissueparams("invivo7T");

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle(npulse,phi0);

    case "ADPCA"
        FA      = [60, 60];        % Flip angles [deg]
        TR      = [100, 20];       % [ms]
        phi0    = 36.0;            % [deg]
        B1range = [0.3,0.5,1,1.3]; % convert such that 100% = 1

        dur1 = 55; % ms
        Gdur{1} = [3,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,30,-30,30];           % [mT/m]
        dur2 = 11; % ms
        Gdur{2} = [1,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = Gamp{1};           % [mT/m]

        % Get tissue parameters
        [T1range,T2range,D] = tissueparams("invivo7T");

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,[1,TR2/TR1]);

    case "KRK"
        FA      = [55,  55]; % Flip angles [deg]
        TR      = [25, 125]; % [ms]

        phi0    = 36;        % [deg]

        B1range = (30:5:140)'/100; % convert such that 100% = 1
        dur1 = 7.2; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,30,-30,30];           % [mT/m]
        dur2 = 36;  % ms
        Gdur{2} = [3,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = Gamp{1};           % [mT/m]

        % Get tissue parameters
        [T1range,T2range,D] = tissueparams("invivo7T");

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,[TR1/TR2,1]);

    case "JS"
        FA      = [60,  60]; % Flip angles [deg]
        TR      = [20, 100]; % [ms]

        phi0    = 36;        % [deg]

        B1range = (30:5:140)'/100; % convert such that 100% = 1
        dur1 = 7.2; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,30,-30,30];           % [mT/m]
        dur2 = 36;  % ms
        Gdur{2} = [3,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = Gamp{1};           % [mT/m]

        % Get tissue parameters
        [T1range,T2range,D] = tissueparams("invivo7T");

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,[TR1/TR2,1]);

    case "BigBrain"
        FA      = [55,  55]; % Flip angles [deg]
        TR      = [25, 125]; % [ms]

        phi0 = 50;        % [deg]

        B1range = (30:5:140)'/100; % convert such that 100% = 1

        dur1 = 7.2; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,30,-30,30];           % [mT/m]
        dur2 = 36;  % ms
        Gdur{2} = [3,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = Gamp{1};           % [mT/m]

        % Get tissue parameters
        [T1range,T2range,D] = tissueparams("postmortem7T");

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,[TR1/TR2,1]);

    case "PVPphantom"
        n = 3;
        FA      = [60, 60]; % Flip angles [deg]
        TR      = [1,n]*50; % [ms]

        phi0    = 50;    % [deg]

        B1range = (50:5:120)'/100; % convert such that 100% = 1
        dur1 = 42; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,26,-26,26];           % [mT/m]
        dur2 = 42; % ms
        Gdur{2} = [n,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = [26,26,-26,26];           % [mT/m]

        % Get tissue parameters
        [T1range,T2range,D] = tissueparams("PVPphantom3T");

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,[TR1/TR2,1]);

    case "IronSleep3T"
        FA      = [60, 60]; % Flip angles [deg]
        TR      = [50,150]; % [ms]

        phi0    = 129.3;    % [deg]

        B1range = (20:5:120)'/100; % convert such that 100% = 1
        dur1 = 42; % ms
        Gdur{1} = [1,dur1/4,dur1/2,dur1/4]; % [ms]
        Gamp{1} = [26,26,-26,26];           % [mT/m]
        dur2 = 42; % ms
        Gdur{2} = [3,dur2/4,dur2/2,dur2/4]; % [ms]
        Gamp{2} = [26,26,-26,26];           % [mT/m]

        % Get tissue parameters
        [T1range,T2range,D] = tissueparams("invivo3T");

        phase_cycle = @(npulse,phi0,TR1,TR2) RF_phase_cycle_NehrkeSimplifiedError(npulse,phi0,[TR1/TR2,1]);
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
B1app_grsp  = zeros([nreps nB1 nT1 nT2]);
for T1idx = 1:nT1 % loop over T1 values, can use parfor for speed

    T1 = T1range(T1idx);
    npulse = 2*nreps; %2*ceil(6*T1/sum(TR)); % ensure steady state signal

    for T2idx = 1:nT2
        T2 = T2range(T2idx);

        for B1idx = 1:nB1  % loop over B1+ values
            B1eff = B1range(B1idx);

            % make train of flip angles and their phases
            alpha_train = repmat(deg2rad(FA*B1eff), 1, npulse/length(FA)); % flip angles
            phi_train   = phase_cycle(npulse,phi0,TR(1),TR(2));            % phases

            % Calculate signals via EPG
            F0 = EPG_GRE_nTR(alpha_train, phi_train, TR, T1, T2, 'diff',Gdiff);

            B1app_grsp(:,B1idx,T1idx,T2idx) = calc_AFI(abs(F0(1:2:end)), abs(F0(2:2:end)), TR(1),TR(2),FA(1));

        end
    end
end


%% Simulate using exact result assuming perfect spoiling
%S1e = abs(hmri_test_utils.dualTRernstd(B1range*FA(1),TR(1),TR(2),1./T1range));
%S2e = abs(hmri_test_utils.dualTRernstd(B1range*FA(1),TR(2),TR(1),1./T1range));
%B1app_compsp = calc_AFI(S1e,S2e,TR(1),TR(2),FA(1));

plot(1:nreps,100*real(B1app_grsp-B1range))
legend("B1 = "+string(100*B1range)+" (p.u.)","Location","best")
xlabel("repetition")
ylabel("B1 estimation error (p.u.)")

end

function B1map = calc_AFI(Y1,Y2,TR1,TR2,nomFA)

% flip angle map in degrees
r=Y2./Y1;
n=TR2/TR1;
FAmap = acosd((r*n-1)./(n-r)); % Eq. (6) in Yarnykh, MRM (2007)

% relative B1 map
B1map = FAmap/nomFA;

end

function [T1,T2,D] = tissueparams(tissuetype)

switch tissuetype
    case "invivo7T"
        T1 = 1200;             % [ms]
        T2 = 30;               % [ms]
        D  = 0.7;              % [µm^2/ms]
    case "invivo3T"
        T1 = [800,1000,1200];             % [ms]
        T2 = 70;               % [ms]
        D  = 0.7;              % [µm^2/ms]
    case "postmortem7T"
        T1 = [500,1000,2000];  % [ms]
        T2 = 30;               % [ms]
        D  = 0.2;              % [µm^2/ms]
    case "PVPphantom3T"
        T1 = 775; % [ms]
        T2 = 250;  % [ms]
        D  = 0.8;  % [µm^2/ms]
    otherwise
        error("unrecognised tissue type %s", tissuetype)
end

end
