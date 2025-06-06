%% Test phase cycling implementation
% test Nehrke AFI RF phase cycling matches standard phase cycling with pseudo-pulses
N = [1,3];
phi0 = 137;
offset = 0;

npulse = 100;
pcN = RF_phase_cycle_Nehrke(npulse,phi0,N(1),N(2),offset);

npulse_pseudo = sum(N*npulse/2);
pcO = wrapTo2Pi(RF_phase_cycle(npulse_pseudo,phi0));
idx = cumsum([1,repmat(N,1,npulse/2)]); idx = idx(1:end-1);
pcO = pcO(idx);

assert(all(pcN(:)-pcO(:)<1e-10))

% test Nehrke AFI RF phase cycling with offset to zero constant phase term
N = [1,3];
phi0 = 137;

npulse = 100;

offset = -phi0*(1-N(2)/N(1))/2;
pcN0 = RF_phase_cycle_Nehrke(npulse,phi0,N(1),N(2),offset);
pcNS = RF_phase_cycle_NehrkeSimplified(npulse,phi0,N(1),N(2));

assert(all(pcN0(:)-pcNS(:)<1e-14))

%% Reproduce Fig 9e,f from Yarnykh (2010), MRM. https://doi.org/10.1002/mrm.22394
fa = 60; % deg
TR = [15, 75; 20, 100]; % ms
Gmax = 25; % mT/m
diffspoil = [280; 450]; % mT ms/m
delta = (diffspoil./Gmax).*[ones(length(diffspoil),1), TR(:,2)./TR(:,1)]; % ms
phi0 = [34,39]; % deg

% WM, GM, MS lesion
params.T1 = [1000, 1500, 1300, 3800];   % ms
params.T2 = [  70,  100,  200, 1900];   % ms
params.D  = [   0.7,  0.8,  1.0,  3.0]; % mT/m

B1 = linspace(10,100)./fa;

figure
tiledlayout(1,2)
for p = 1:length(phi0)
    B1est = zeros(length(params.T1),length(B1));
    for n = 1:length(params.T1)
        npulse = 2*ceil(6*params.T1(n)/sum(TR(p,:))); % ensure steady state signal
        Gdiff = struct('D', params.D(n)*1e-9, 'G', Gmax, 'tau', num2cell(delta(p,:))); % struct assigns cell elements to separate struct array elements

        for b = 1:length(B1)
            alpha_train = repmat(deg2rad(fa*B1(b)), 1, npulse); % flip angles
            phi_train   = RF_phase_cycle(npulse,phi0(p));          % phases

            % Calculate signals via EPG
            F0 = EPG_GRE_nTR(alpha_train, phi_train, TR(p,:), params.T1(n), params.T2(n), 'diff',Gdiff);
            S1 = abs(F0(end-1));
            S2 = abs(F0(end));

            B1est(n,b) = calc_AFI(S1,S2,TR(p,1),TR(p,2),fa);
        end
    end

    nexttile
    plot(B1*fa, 100*(B1est-B1)./B1)
    ylim([-12,5])
    xlabel("flip angle (°)")
    ylabel("B1 relative error (%)")

end

%% Reproduce Fig 7a from Yarnykh (2010), MRM. https://doi.org/10.1002/mrm.22394
fa = 53.4; % deg
TR = [10, 50]; % ms
Gmax = 25; % mT/m
diffspoil = 11.7; % mT ms/m
delta = (diffspoil./Gmax).*[1, TR(2)/TR(1)]; % ms
phi0 = 0:3:180; % deg

% phantom parameters
params.T1 = 762;
params.T2 = 644;
params.D = 2.2;

B1 = 1;

npulse = 2*ceil(6*params.T1/sum(TR)); % ensure steady state signal
Gdiff = struct('D', params.D*1e-9, 'G', Gmax, 'tau', num2cell(delta)); % struct assigns cell elements to separate struct array elements

figure
S1 = zeros(length(phi0),length(B1));
S2 = zeros(length(phi0),length(B1));
for p = 1:length(phi0)
    for b = 1:length(B1)
        alpha_train = repmat(deg2rad(fa*B1(b)), 1, npulse); % flip angles
        phi_train   = RF_phase_cycle(npulse,phi0(p));          % phases

        % Calculate signals via EPG
        F0 = EPG_GRE_nTR(alpha_train, phi_train, TR, params.T1, params.T2, 'diff',Gdiff);
        S1(p,b) = abs(F0(end-1));
        S2(p,b) = abs(F0(end));

    end
end
S1spoiled = abs(hmri_test_utils.dualTRernstd(fa,TR(1),TR(2),1/params.T1));
S2spoiled = abs(hmri_test_utils.dualTRernstd(fa,TR(2),TR(1),1/params.T1));
Sref = abs(hmri_test_utils.dualTRernstd(fa,25,125,1/params.T1)); % diffusion spoiled result as reference

B1est = calc_AFI(S1,S2,TR(1),TR(2),fa);

plot(phi0, [S1,S2]/Sref)
hold on
plot(xlim(),[1,1]*S1spoiled/Sref)
plot(xlim(),[1,1]*S2spoiled/Sref)
ylim([0.15,1.4])
xlabel("flip angle (°)")
ylabel("B1 relative error (%)")

%%
function B1map = calc_AFI(S1,S2,TR1,TR2,nomFA)

% flip angle map in degrees
r=S2./S1;
n=TR2/TR1;
FAmap = acosd((r*n-1)./(n-r)); % Eq. (6) in Yarnykh, MRM (2007)

% relative B1 map
B1map = FAmap/nomFA;

end
