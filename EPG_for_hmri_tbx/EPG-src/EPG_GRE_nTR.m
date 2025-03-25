function [F0,Fn,Zn,F] = EPG_GRE_nTR(theta,phi,TR,T1,T2,varargin)
%   [F0,Fn,Zn,F] = EPG_GRE(theta,phi,TR,T1,T2,varargin)
%
%   Single pool EPG (classic version) for gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               phi:        phase per pulse. This function can hence be
%                           used to simulate RF spoiling or balanced
%                           sequences depending on how phase is cycled
%                           see function RF_phase_cycle()
%               TR:         repetition times, ms
%               T1:         T1, ms
%               T2:         T2, ms
%
%   optional arguments (use string then value as next argument)
%
%               kmax:       maximum EPG order to include. Can be used to
%                           accelerate calculation. 
%                           Setting kmax=inf ensures ALL pathways are
%                           computed
%               diff:       structure with fields:
%                           G    - Gradient amplitude(s) mT/m
%                           tau  - Gradient durations(s) ms
%                           D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%
%   Outputs:                
%               F0:         signal (F0 state) directly after each
%                           excitation
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0 F0* Z0 F1 F-1* Z1 F2 F-2* Z2 ...] etc
%
%
%   Shaihan Malik 2017-07-20


%% Extra variables
for ii=1:length(varargin)
    
    % Kmax = this is the maximum EPG 'order' to consider
    % If this is infinity then don't do any pruning
    if strcmpi(varargin{ii},'kmax')
        kmax = varargin{ii+1};
    end
    
    % Diffusion - structure contains, G, tau, D
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
    
end

% Different TRs might have different amounts of spoiling, which affects how
% far we need to move in k-space   
np = length(theta);
ntr = length(TR);
nShifts = ones(1,ntr); % default to implicitly having the same amount of spoiling every TR
if exist('diff','var')
    G0 = zeros(ntr,1);
    for tridx=1:ntr
        G0(tridx) = dot(diff(tridx).G(:),diff(tridx).tau(:));
    end

    % implementing this for now assuming that gradient moments are all zero or an integer 
    % multiple of the smallest moment
    nShifts(G0~=0) = G0(G0~=0)/min(G0(G0~=0));
    assert(all(mod(nShifts,1)==0), 'gradient moments per TR are not all integer multiples of the shortest non-zero moment')

    kall = sum(repmat(nShifts,1,np/ntr));
else
    kall = np - 1;
end

%%% The maximum order varies through the sequence. This can be used to speed up the calculation 
% if not defined, assume want max
if ~exist('kmax','var')
    kmax = kall;
end

if isinf(kmax)
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = kall; 
else
    allpathways = false;
end

%%% Variable pathways
if allpathways
    kmax_per_pulse = cumsum(repmat(nShifts,1,np/ntr)); %<-- +1 because (0:kmax) is correct after each RF pulse, but we must increase order by one to also deal with subsequent shift
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax; %<-- don't exceed kmax, we break after last RF pulse
else
    error('not implemented')
    %{
    kmax_per_pulse = [1:ceil(np/2) (floor(np/2)):-1:1];
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
     
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
    %}
end

%%% Number of states is 6x(kmax +1) -- +1 for the zero order
N=3*(kmax+1);

%%% Build Shift matrix, S
S0 = sparse(EPG_shift_matrices(kmax));
S = cell(ntr,1);
for tridx=1:ntr
    S{ntr} = S0^nShifts(ntr);
end


%% Set up matrices for Relaxation
ES = cell(ntr,1);
b  = cell(ntr,1);
for tridx=1:ntr
    E1 = exp(-TR(tridx)/T1);
    E2 = exp(-TR(tridx)/T2);
    E = diag([E2 E2 E1]);

    %%% regrowth
    b{tridx} = zeros([N 1]);
    b{tridx}(3) = 1-E1;%<--- just applies to Z0

    %%% Add in diffusion at this point 
    if exist('diff','var')
        E = E_diff(E,diff(tridx),kmax,N);
    else
        % If no diffusion, E is the same for all EPG orders
        E = spdiags(repmat([E2 E2 E1],[1 kmax+1])',0,N,N);
    end
        
    %%% Composite relax-shift
    ES{tridx}=sparse(E*S{ntr});
end

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

% store the indices of the top 3x3 corner, this helps build_T
i1 = [];
for ii=1:3
    i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3)));
end


%% F matrix (many elements zero, not efficient)
F = zeros([N np]); %%<-- records the state after each RF pulse 

%%% Initial State
FF = zeros([N 1]);
FF(3)=1;   % M0 - could be variable


%% Main body of gradient echo sequence, loop over TRs 

for jj=1:np 
    %%% RF transition matrix
    A = RF_rot(theta(jj),phi(jj));
   
    %%% Variable order of EPG, speed up calculation
    kmax_current = kmax_per_pulse(jj);
    kidx = 1:3*(kmax_current+1); %+1 because states start at zero
    
    %%% Replicate A to make large transition matrix
    build_T(A);
    
    %%% Apply flip and store this: splitting these large matrix
    %%% multiplications into smaller ones might help
    F(kidx,jj)=T(kidx,kidx)*FF(kidx);
    
    if jj==np
        break
    end
    
    %%% Now deal with evolution
    tridx = mod(jj-1,ntr)+1;
    FF(kidx) = ES{tridx}(kidx,kidx)*F(kidx,jj)+b{tridx}(kidx);
    
    % Deal with complex conjugate after shift
    FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate 
end


%%% Return signal
F0=F(1,:);

%%% phase demodulate
F0 = F0(:) .* exp(-1i*phi(:)) *1i;


%%% Construct Fn and Zn
idx=[fliplr(5:3:size(F,1)) 1 4:3:size(F,1)]; 
kvals = -kmax:kmax;

%%% Now reorder
Fn = F(idx,:);
%%% Conjugate
Fn(kvals<0,:)=conj(Fn(kvals<0,:));

%%% Similar for Zn
Zn = F(3:3:end,:);



    %%% NORMAL EPG transition matrix as per Weigel et al JMR 2010 276-285 
    function Tap = RF_rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
    end

    function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end
