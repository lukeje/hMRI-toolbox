function hmri_corr_imperf_spoil_afisim(job)
%==========================================================================
% PURPOSE
% Compute coefficients to correct for the effect of imperfect spoiling on
% T1 estimation as described in (Preibisch & Deichmann, MRM 2009).
%
% Numerical simulations are performed using the EPG formalism described in
% Malik et al., MRM 2017 and available here:
% https://github.com/mriphysics/EPG-X
%
% The parameters used for the simulation and the resulting correction
% factors are written to a JSON file in the output folder specified by the
% user.
%==========================================================================

hmri_log(sprintf('\t--- Calculating Imperfect Spoiling Correction Coefficients ---'));
%% ***********************************************%%
% 1./ Numerical simulations with EPG
%*************************************************%%
%%
% Get sequence parameters
FA      =   job.seq_params.FA_deg;              % Flip angles [deg]
TR      =   job.seq_params.TR_ms;               % [ms]
Phi0    =   job.seq_params.Phi0_deg;            % [deg]
Gdur    =   job.seq_params.Gdur_ms;             % [ms]
Gamp    =   job.seq_params.Gamp_mT_per_m;       % [mT/m]

B1range =   job.B1range_percent/100; % convert such that 100% = 1

assert(length(Gdur) == length(Gamp), 'The vectors of gradient durations and amplitudes must have the same length!')
assert(all(sum(Gdur)<=TR), 'The total duration of the gradients cannot exceed TR!')

% Get AFI parameters
FA_afi      = job.afi_params.FA_deg;             % Flip angles [deg]
TR_afi      = job.afi_params.TR_ms;              % [ms]
Phi0_afi    = job.afi_params.Phi0_deg;           % [deg]
Phi0_type   = job.afi_params.rf_spoiling_type;
Gdur_afi{1} = job.afi_params.Gdur_ms_1;          % [ms]
Gamp_afi{1} = job.afi_params.Gamp_mT_per_m_1;    % [mT/m]
Gdur_afi{2} = job.afi_params.Gdur_ms_2;          % [ms]
Gamp_afi{2} = job.afi_params.Gamp_mT_per_m_2;    % [mT/m]

assert(TR_afi(1)~=TR_afi(2), "AFI TRs cannot be equal!")

%% Get tissue parameters
T1range     = job.tissue_params.T1range_ms;     %[ms]
T2range     = job.tissue_params.T2range_ms;     % [ms]
D           = job.tissue_params.D_um2_per_ms;   % [um^2/ms]

%% Build structure "diff" to account for diffusion effect
% Note we include any deadtime during each TR so that diffusion effects
% are calculated correctly
for n=2:-1:1 % go backwards to avoid matlab warning about preallocation
    diff(n).D   = D*1e-9;
    diff(n).G   = [0; Gamp(:)];
    diff(n).tau = [TR(n)-sum(Gdur); Gdur(:)];
end

% AFI typically has different amounts of spoiling in each TR
for n=2:-1:1 % go backwards to avoid matlab warning about preallocation
    diff_afi(n).D   = D*1e-9;
    diff_afi(n).G   = [0; Gamp_afi{n}(:)];
    diff_afi(n).tau = [TR_afi(n)-sum(Gdur_afi{n}); Gdur_afi{n}(:)];
end

%% Run EPG simulation
nT1  = length(T1range);
nT2  = length(T2range);
nB1  = length(B1range);
S1   = zeros([nT1 nT2 nB1]);
S2   = zeros([nT1 nT2 nB1]);
AFI1 = zeros([nT1 nT2 nB1]);
AFI2 = zeros([nT1 nT2 nB1]);
hmri_log(sprintf('\t-------- Simulating signals'));
for T1val = 1 : nT1 % loop over T1 values, can use parfor for speed

    T1 = T1range(T1val);

    npulse = floor(15*T1/min(TR));   % ensure steady state signal
    phi_train = RF_phase_cycle(npulse,Phi0); % phase of the RF pulses

    npulse_afi = floor(15*T1/sum(TR_afi));   % ensure steady state signal
    npulse_afi = npulse_afi + (mod(npulse_afi,2)); % ensure the number of afi TRs is even
    switch lower(Phi0_type)
        case 'standard'
            phi_train_afi = RF_phase_cycle(npulse_afi,Phi0_afi); % phase of the RF pulses
        case 'nehrke'
            if TR1>TR2
                N1 = TR1/TR2;
                N2 = 1;
            elseif TR1<TR2
                N1 = 1;
                N2 = TR2/TR1;
            else
                error("AFI TRs should not be equal!")
            end

            if mod(N1,1)~=0 || mod(N2,1) ~=0
                warning("This function expects that the larger AFI TR is an integer multiple of the smaller TR!")
            end
            phi_train_afi = RF_phase_cycle_Nehrke(npulse_afi,Phi0_afi,N1,N2); % phase of the RF pulses
    end

    for T2val = 1 : nT2
        T2 = T2range(T2val);

        for B1val = 1 : nB1  % loop over B1+ values
            B1eff = B1range(B1val);

            %% Calculate MPM signals via EPG:
            % make train of flip angles
            alpha_train1 = d2r(FA(1)*B1eff)*ones([1 npulse]); % flip angles of the PDw acquisitions
            alpha_train2 = d2r(FA(2)*B1eff)*ones([1 npulse]); % flip angles of the T1w acquisitions

            %PDw
            F0 = EPG_GRE(alpha_train1, phi_train, TR(1), T1, T2, 'diff', diff(1));
            S1(T1val,T2val,B1val) = abs(F0(end));
            %T1w
            F0 = EPG_GRE(alpha_train2, phi_train, TR(2), T1, T2, 'diff', diff(2));
            S2(T1val,T2val,B1val) = abs(F0(end));

            %% Calculate AFI signals via EPG:
            alpha_train = d2r(FA_afi*B1eff)*ones([1 npulse_afi]); % flip angles of the AFI acquisitions
            F0 = EPG_GRE_nTR(alpha_train, phi_train_afi, TR_afi, T1, T2, 'diff', diff_afi);
            AFI1(T1val,T2val,B1val) = abs(F0(end-1));
            AFI2(T1val,T2val,B1val) = abs(F0(end));
        end
    end
end



%% ***********************************************%%
% 2./ Fitting T1=A(B1eff)+B(B1eff)*T1app
%*************************************************%%
hmri_log(sprintf('\t-------- Determining Coefficients'));
T1app = zeros(nT1, nT2, nB1);
B1app = zeros(nT1, nT2, nB1);
for B1val = 1:nB1
    B1app(:,:,B1val) = 0.01*hmri_calc_AFI_B1map(AFI1(:,:,B1val),AFI2(:,:,B1val),TR_afi(2)/TR_afi(1),FA_afi);

    % Calculate T1app, accounting for B1+
    T1app(:,:,B1val) = 1./hmri_calc_R1(...
        struct('data',S1(:,:,B1val),'fa',d2r(FA(1)),'TR',TR(1),'B1',B1app(:,:,B1val)),...
        struct('data',S2(:,:,B1val),'fa',d2r(FA(2)),'TR',TR(2),'B1',B1app(:,:,B1val)),...
        job.small_angle_approx);
end

%% *********************************************************%%
% 3./ Fitting A=P(B1eff) and B=P(B1eff) with 2nd degree polynomial
%***********************************************************%%
X = [B1app(:).^2, B1app(:), ones(nT1*nT2*nB1,1)]; % quadratic polynomial in B1app argument
X = [X, X.*T1app(:)]; % linear polynomial in T1app argument
coeff = X\repmat(T1range(:),nT2*nB1,1);
polyCoeffA = coeff(1:end/2);
polyCoeffB = coeff(end/2+1:end);

%% *********************************************************%%
% 4./ Compute RMSE on T1app and T1
%***********************************************************%%
T1corr = polyval(polyCoeffA, B1app) + polyval(polyCoeffB, B1app).*T1app;
T1_Corr_Err = 100*(T1corr - T1range(:))./T1range;
T1_App_Err  = 100*(T1app  - T1range(:))./T1range;

RMSE_Corr = rms(T1_Corr_Err(:));
RMSE_App  = rms(T1_App_Err(:));

%% *********************************************************%%
% 5./ Write parameters and correction factors in a json file
% in the selected output directory
%***********************************************************%%
hmri_log(sprintf('\t-------- Writing results\n'));
Results.Input = job;
Results.Output.P2_a = round(polyCoeffA,4);
Results.Output.P2_b = round(polyCoeffB,4);
Results.Output.small_angle_approx = job.small_angle_approx;
Results.Output.RMSE_percent.T1app=round(RMSE_App,3);
Results.Output.RMSE_percent.T1corr=round(RMSE_Corr,3);

Results.ToCopy{1}    =['hmri_def.MPMacq_set.names{NN} = ''' job.prot_name ''';' ];
Results.ToCopy{end+1}=['hmri_def.MPMacq_set.tags{NN}  = ''' strrep(job.prot_name,' ','') ''';'];
Results.ToCopy{end+1}=['hmri_def.MPMacq_set.vals{NN}  = [' num2str([TR FA]) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.tag = ''' strrep(job.prot_name,' ','') ''';' ];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_a = [' num2str(round(polyCoeffA',4)) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.P2_b = [' num2str(round(polyCoeffB',4)) '];'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.small_angle_approx = ' mat2str(job.small_angle_approx) ';'];
Results.ToCopy{end+1}=['hmri_def.imperfectSpoilCorr.' strrep(job.prot_name,' ','') '.enabled = hmri_def.imperfectSpoilCorr.enabled;'];

results_filename = fullfile(job.outdir,[strrep(job.prot_name,' ',''),'.json']);

spm_jsonwrite(results_filename{1},Results,struct('indent','\t'));

end