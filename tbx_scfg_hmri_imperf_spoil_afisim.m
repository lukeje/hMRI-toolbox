function imperf_spoil=tbx_scfg_hmri_imperf_spoil_afisim
% 
% PURPOSE: Compute correction factors for imperfect spoiling based on the 
% method described in (Preibisch & Deichmann, MRM, 2009). 
%
% METHODS: Numerical simulations are performed with the Extended Phase 
% Graph framework(code from https://github.com/mriphysics/EPG-X). Sequence 
% parameters, a range of values for B1+ efficiency and T1, and global 
% parameters T2 and diffusion coefficent D must be provided. 
%
%_______________________________________________________________________
% Wellcome Centre for Human Neuroimaging
% Nadège Corbin - May 2020
% ======================================================================

% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where a json file containing the correction parameters will be written to.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

% ---------------------------------------------------------------------
% Name of the protocol
% ---------------------------------------------------------------------
prot_name        = cfg_entry;
prot_name.tag     = 'prot_name';
prot_name.name    = 'Protocol Name ';
prot_name.val     = {'Unit_Test_Protocol_AFI'};
prot_name.strtype = 's';
prot_name.help    = {'Specify the name of the protocol'};

% ---------------------------------------------------------------------
% T1 range [ms]
% ---------------------------------------------------------------------
T1range        = cfg_entry;
T1range.tag     = 'T1range_ms';
T1range.name    = 'T1 range (ms)';
T1range.val     = {500:100:2000};
T1range.strtype = 'e';
T1range.num     = [1 Inf];
T1range.help    = {['Specify the range of T1 times over which the ',...
    'corrections factors will be computed. A linear fitting ',...
    'T1 = A + B*T1app will be performed to estimate A and B for each B1+ '...
    'efficiency.']};

% ---------------------------------------------------------------------
% B1+ efficiency range [%]
% ---------------------------------------------------------------------
B1range         = cfg_entry;
B1range.tag     = 'B1range_percent';
B1range.name    = 'Expected B1+ range (%)';
B1range.val     = {70:5:130};
B1range.strtype = 'e';
B1range.num     = [1 Inf];
B1range.help    = {['Specify the range of transmit field efficiency (B1+) over which the ',...
    'corrections factors will be computed. After the linear fitting ',...
    'T1 = A + B*T1app , a polynomial fitting will be perfomed to estimate ',...
    'A and B such that: A=P(B1) and B=P(B1) with P a 2nd degree polynomial.']};

% ---------------------------------------------------------------------
% T2 [ms] 
% ---------------------------------------------------------------------
T2        = cfg_entry;
T2.tag     = 'T2range_ms';
T2.name    = 'T2 range (ms)';
T2.val     = {70};
T2.strtype = 'r';
T2.num     = [1 Inf];
T2.help    = {'Specify an estimate of the T2 time in ms, or an array of values'};

% ---------------------------------------------------------------------
% D [um^2/ms] 
% ---------------------------------------------------------------------
D        = cfg_entry;
D.tag     = 'D_um2_per_ms';
D.name    = 'D (um^2/ms)';
D.val     = {0.8};
D.strtype = 'r';
D.num     = [1 1];
D.help    = {'Specify an estimate of the diffusion coeffcient (D) in um^2/ms'};

% ---------------------------------------------------------------------
% Readout gradient amplitude [ms] 
% ---------------------------------------------------------------------
Gdur        = cfg_entry;
Gdur.tag     = 'Gdur_ms';
Gdur.name    = 'Spoiler gradient duration (ms)';
Gdur.val     = {3.3876};
Gdur.strtype = 'e';
Gdur.num     = [1 Inf];
Gdur.help    = {['Specify the duration (in ms) of the spoiler gradient  ',...
                'of the FLASH acquisitions. Note here a vector could be '...
                'included to account for the full effect of the readout, e.g. multiple echoes ']};

% ---------------------------------------------------------------------
% Spoiler gradient amplitude [ms] 
% ---------------------------------------------------------------------
Gamp        = cfg_entry;
Gamp.tag     = 'Gamp_mT_per_m';
Gamp.name    = 'Spoiler gradient amplitude (mT/m)';
Gamp.val     = {26.00};
Gamp.strtype = 'e';
Gamp.num     = [1 Inf];
Gamp.help    = {['Specify the amplitude (in mT/m) of the spoiling gradient ',...
                'of the FLASH acquisitions ']};

% ---------------------------------------------------------------------
% TR [ms] 
% ---------------------------------------------------------------------
TR         = cfg_entry;
TR.tag     = 'TR_ms';
TR.name    = 'TR';
TR.val     = {[25 25]};
TR.strtype = 'r';
TR.num     = [1 2];
TR.help    = {['Specify the TR (in ms) of the PD-weighted ',...
    'and the T1-weighted FLASH acquisitions (in that order)']};

% ---------------------------------------------------------------------
% RF Spoiling increment [deg]
% ---------------------------------------------------------------------
Phi0         = cfg_entry;
Phi0.tag     = 'Phi0_deg';
Phi0.name    = 'RF spoiling increment';
Phi0.val     = {137};
Phi0.strtype = 'w';
Phi0.num     = [1 1];
Phi0.help    = {'Specify the RF Spoiling increment (in deg) of the FLASH acquisitions'};

% ---------------------------------------------------------------------
% Flip angles [deg]
% ---------------------------------------------------------------------
FA        = cfg_entry;
FA.tag     = 'FA_deg';
FA.name    = 'Flip angles';
FA.val     = {[6 21]};
FA.strtype = 'e';
FA.num     = [1 2];
FA.help    = {['Specify the flip angles (in deg) of the PD-weighted ',...
    'and the T1-weighted FLASH acquisitions (in that order)']};

% ---------------------------------------------------------------------
% Adapt parameters for AFI
% ---------------------------------------------------------------------
TR_AFI = TR;
TR_AFI.val = {[20 100]};
TR_AFI.help = {'Specify the first and second TR (in ms) of the AFI acquisition'};

FA_AFI         = cfg_entry;
FA_AFI.tag     = 'FA_deg';
FA_AFI.name    = 'Flip angle';
FA_AFI.val     = {60};
FA_AFI.strtype = 'e';
FA_AFI.num     = [1 1];
FA_AFI.help    = {'Specify the flip angle (in deg) of the AFI acquisition'};

Phi0_AFI = Phi0;
Phi0_AFI.val  = {50};
Phi0_AFI.help = {'Specify the RF Spoiling increment (in deg) of the AFI acquisitions'};

Phi0_type        = cfg_menu;
Phi0_type.tag    = 'rf_spoiling_type';
Phi0_type.name   = 'RF spoiling type';
Phi0_type.help   = {'Whether standard RF spoiling or Nehrke spoiling was used.'};
Phi0_type.labels = { 'standard', 'Nehrke' };
Phi0_type.values = { 'standard', 'Nehrke' };
Phi0_type.val    = { 'standard' };

Gamp1 = Gamp;
Gamp1.tag  = [Gamp.tag '_1'];
Gamp1.name = [Gamp.name ' of first AFI TR'];
Gamp1.val  = {40};
Gamp1.help = {['Specify the amplitude (in mT/m) of the spoiling gradient ',...
                'of the first AFI TR']};
Gamp2 = Gamp;
Gamp2.tag = [Gamp.tag '_2'];
Gamp2.name = [Gamp.name ' of second AFI TR'];
Gamp2.val  = {40};
Gamp2.help = {['Specify the amplitude (in mT/m) of the spoiling gradient ',...
                'of the second AFI TR']};

Gdur1 = Gdur;
Gdur1.tag = [Gdur.tag '_1'];
Gdur1.name = [Gdur.name ' of first AFI TR'];
Gdur1.val  = {1};
Gdur1.help = {['Specify the duration (in ms) of the spoiler gradient  ',...
               'of the first AFI TR. Note here a vector could be '...
               'included to account for the full effect of the readout, e.g. multiple echoes ']};
Gdur2 = Gdur;
Gdur2.tag = [Gdur.tag '_2'];
Gdur2.name = [Gdur.name ' of second AFI TR'];
Gdur2.val  = {5};
Gdur2.help = {['Specify the duration (in ms) of the spoiler gradient  ',...
               'of the second AFI TR. Note here a vector could be '...
               'included to account for the full effect of the readout, e.g. multiple echoes ']};

% ---------------------------------------------------------------------
% All tissue parameters
% ---------------------------------------------------------------------
tissue_params            = cfg_branch;
tissue_params.tag        = 'tissue_params';
tissue_params.name       = 'Tissue parameters';
tissue_params.help       = {'Input all the tissue parameters.'};
tissue_params.val        = {T1range T2 D};

% ---------------------------------------------------------------------
% All MPM sequence parameters
% ---------------------------------------------------------------------
seq_params            = cfg_branch;
seq_params.tag        = 'seq_params';
seq_params.name       = 'MPM sequence parameters';
seq_params.help       = {'Input all the MPM sequence parameters.'};
seq_params.val        = {FA TR Phi0 Gdur Gamp};

% ---------------------------------------------------------------------
% All AFI sequence parameters
% ---------------------------------------------------------------------
afi_params            = cfg_branch;
afi_params.tag        = 'afi_params';
afi_params.name       = 'AFI sequence parameters';
afi_params.help       = {'Input all the AFI sequence parameters.'};
afi_params.val        = {FA_AFI TR_AFI Phi0_AFI Phi0_type Gdur1 Gamp1 Gdur2 Gamp2};

% ---------------------------------------------------------------------
% Approximation parameters
% ---------------------------------------------------------------------
small_angle_approx        = cfg_menu;
small_angle_approx.tag    = 'small_angle_approx';
small_angle_approx.name   = 'Small angle approximation';
small_angle_approx.help   = {'Should the small angle approximation be used to compute T1? This can lead to significant differences in the computed parameters.'};
small_angle_approx.labels = { 'yes', 'no' };
small_angle_approx.values = { true false };
small_angle_approx.val    = { false };

% ---------------------------------------------------------------------
% Optimally plot residuals
% ---------------------------------------------------------------------
plot_residuals        = cfg_menu;
plot_residuals.tag    = 'plot_residuals';
plot_residuals.name   = 'Plot residuals';
plot_residuals.help   = {'Should the residual errors be plotted for QA?'};
plot_residuals.labels = { 'yes', 'no' };
plot_residuals.values = { true false };
plot_residuals.val    = { false };

% ---------------------------------------------------------------------
% Compute correction factors for imperfect spoiling
% ---------------------------------------------------------------------
imperf_spoil         = cfg_exbranch;
imperf_spoil.tag     = 'imperf_spoil_afisim';
imperf_spoil.name    = 'Imperfect Spoiling Calc. with AFI Correction';
imperf_spoil.val     = { outdir prot_name seq_params afi_params tissue_params B1range small_angle_approx plot_residuals };
imperf_spoil.help    = {'Given input info about the sequence settings and expected tissue properties, ' ...
    'this module computes coefficients required to correct for imperfect spoiling in the FLASH volumes ' ...
    'using the method proposed by Preibisch & Deichmann, MRM 2009, 61(1):125'};
imperf_spoil.prog    = @hmri_corr_imperf_spoil_afisim;
