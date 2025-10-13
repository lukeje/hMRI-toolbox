function F0 = EPG_GRE_nTR(theta,phi,TR,T1,T2,varargin)
%   F0 = EPG_GRE_nTR(theta,phi,TR,T1,T2,varargin)
%
%   Single pool EPG (classic version) for gradient echo sequences with
%   different interleaved TRs
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
%               diff:       cell array of struct array with fields:
%                           G    - Gradient amplitude(s) mT/m
%                           tau  - Gradient durations(s) ms
%                           D    - Diffusion coeff       m^2/s (i.e. expect 10^-9)
%                           each element of the cell array represents a different 
%                           gradient axis, and each element of the struct array
%                           is a different TR
%
%   Outputs:                
%               F0:         signal (F0 state) directly after each
%                           excitation
%
%
%   Adapted from Shaihan Malik's EPG_GRE.m by Luke J. Edwards


%% Extra variables
for ii=1:length(varargin)
    
    % kmax = this is the maximum EPG 'order' to consider
    % If this is infinity then don't do any pruning
    if strcmpi(varargin{ii},'kmax')
        kmaxin = varargin{ii+1};
    end
    
    % Diffusion - structure contains, G, tau, D
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
        if isstruct(diff)
            diff = {diff};
        elseif ~iscell(diff)
            error("diff argument must either be a struct or a cell array of structs (one per gradient axis)")
        end
        assert(length(diff)<=3, "there can only be up to three separate gradient axes!")
    end
    
end

% Different TRs might have different amounts of spoiling, which affects how
% far we need to move in k-space   
ntr = length(TR);
if exist('diff','var')
    ngaxes  = length(diff);
    nshifts = cell(1,length(diff));
    dk      = zeros(1,length(diff));
    for d=1:ngaxes
        diff{d} = filltr(diff{d},TR);
        [nshifts{d},dk(d)] = computeshifts(diff{d});
    end
else
    % default to implicitly having the same amount of spoiling every TR
    ngaxes = 1;
    nshifts = {ones(1,ntr)};
    dk = {[]};
end

% maximum k which can be reached in the simulation
np = length(theta);
kall = zeros(1,3); % up to 3 gradient axes
allshifts = cell(1,ngaxes);
for d=1:ngaxes
    allshifts{d} = repmat(nshifts{d},1,ceil(np/ntr)); % all shifts if we always complete the TR cycle
    allshifts{d} = allshifts{d}(1:np);   % all the shifts actually performed
    kall(d) = sum(allshifts{d}(1:np-1)); % ignore last shift as we break after last pulse
end

%%% The maximum order varies through the sequence. This can be used to speed up the calculation 
% if not defined, assume want max
kmax = zeros(1,3); % up to 3 gradient axes
if ~exist('kmaxin','var')
    kmax = kall;
elseif isscalar(kmaxin)
    kmax(1:ngaxes) = kmaxin;
elseif length(kmaxin)<ngaxes
    error("please specify kmax either as a scalar or one value per gradient axis!")
elseif length(kmaxin)>ngaxes
    error("too many elements in kmax! It cannot be greater than the number of gradient axes!")
else
    kmax(1:ngaxes) = kmaxin;
end

if any(isinf(kmax))
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = kall;
else
    allpathways = false;
end

%%% Variable pathways
kmax_per_pulse = repmat({zeros(np,1)},1,3); % 3 as up to 3 gradient axes
for d=1:ngaxes
    if allpathways
        kmax_per_pulse{d} = cumsum(allshifts{d}); % current max state plus subsequent shift
        kmax_per_pulse{d}(kmax_per_pulse{d}>kmax(d))=kmax(d); % don't exceed kmax as we break after last RF pulse
    else
        % reduce the number of required states by using the fact that states must be refocused
        % last shift is zero as we break after last pulse
        kmax_per_pulse{d} = min(cumsum(allshifts{d}),cumsum([allshifts{d}(1:end-1),0],'reverse'));
        kmax_per_pulse{d}(kmax_per_pulse{d}>kmax(d))=kmax(d);
        kmax(d) = min(max(kmax_per_pulse{d}),kmax(d));
    end
end


%%% Number of states is 3x(2*kmax +1) -- +1 for the zero order
N = 3*prod(2*kmax+1);
iF0  = sub2ind([3,2*kmax+1],1,kmax(1)+1,kmax(2)+1,kmax(3)+1);
iF0z = sub2ind([3,2*kmax+1],3,kmax(1)+1,kmax(2)+1,kmax(3)+1);


%%% Build Shift matrices, S
S = cell(1,ntr);
for tridx=1:ntr
    S0 = EPG_shift_matrices(kmax(1),true);
    S{tridx} = S0^nshifts{1}(tridx);
    for d=2:ngaxes
        % remember to permute vectors so that S0 operates on the correct g2 x spin dimension
        % pre:  g2 x g1 x spin -I(g2)xP(g1,spin)-> g2 x spin x g1
        % post: g2 x spin x g1 -I(g2)xP(spin,g1)-> g2 x g1 x spin
        S0 = EPG_shift_matrices(kmax(d),true);
        Pre   = kron(speye(2*kmax(d)+1), permuteKron(prod(2*kmax(1:(d-1))+1), 3));
        Post  = kron(speye(2*kmax(d)+1), permuteKron(3, prod(2*kmax(1:(d-1))+1)));        
        S{tridx} = Post*kron(S0^nshifts{d}(tridx), speye(prod(2*kmax(1:(d-1))+1)))...
            *Pre*kron(speye(2*kmax(d)+1), S{tridx});
    end
end


%% Set up matrices for Relaxation
SE = cell(1,ntr);
b  = cell(1,ntr);
for tridx=1:ntr
    E1 = exp(-TR(tridx)/T1);
    E2 = exp(-TR(tridx)/T2);
    E = diag([E2 E2 E1]);

    %%% regrowth
    b{tridx} = sparse(iF0z,1,1-E1,N,1); % just applies to Z0

    %%% Add in diffusion at this point 
    if exist('diff','var')
        for d=1:ngaxes
            E = E_diff(E,diff{d}(tridx),kmax(d),3*prod(2*kmax(1:d)+1),dk(d),true);
        end
    else
        % If no diffusion, E is the same for all EPG orders
        E = spdiags(repmat([E2 E2 E1],[1 prod(2*kmax+1)])',0,N,N);
    end
        
    %%% Composite relax-shift
    SE{tridx}=S{tridx}*E;
end


%% F matrix (many elements zero, not efficient)
F = zeros([N np]); %%<-- records the state after each RF pulse 

%%% Initial State
FF = zeros([N 1]);
FF(iF0z)=1; % M0 - could be variable


%% Main body of gradient echo sequence, loop over TRs 

for jj=1:np 
    %%% RF transition matrix
    A = RF_rot(theta(jj),phi(jj));
   
    %%% Variable order of EPG, speed up calculation
    %+1 because states start at zero
    [i,j,k,l] = ndgrid(1:3,-kmax_per_pulse{1}(jj):kmax_per_pulse{1}(jj),-kmax_per_pulse{2}(jj):kmax_per_pulse{2}(jj),-kmax_per_pulse{3}(jj):kmax_per_pulse{3}(jj));
    kidx = sub2ind([3,2*kmax+1],i(:),1+kmax(1)+j(:),1+kmax(2)+k(:),1+kmax(3)+l(:));
    
    %%% Replicate A to make large transition matrix
    T = kron(speye(prod(2*kmax+1)), build_T(sparse(3,3),A,0));
    
    %%% Apply flip and store this: splitting these large matrix
    %%% multiplications into smaller ones might help
    F(kidx,jj)=T(kidx,kidx)*FF(kidx);
    
    if jj==np
        break
    end
    
    %%% Now deal with evolution
    tridx = mod(jj-1,ntr)+1;
    FF(kidx) = SE{tridx}(kidx,kidx)*F(kidx,jj)+b{tridx}(kidx);
end


%%% Return signal
F0=F(iF0,:);

%%% phase demodulate
F0 = F0(:) .* exp(-1i*phi(:)) *1i;


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

    function T = build_T(T,AA,kmax)
        ksft = 3*(3*(kmax+1)+1);
        i1 = 1:9;
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end
    
end

function diff = filltr(diff,TR)
    % fill out TRs to ensure b-values computed correctly
    % assume gradients are spoilers played out at the end of the TR
    ntr = length(TR);
    assert(length(diff) == ntr, "each TR must have an associated set of diffusion parameters!")
    for tridx=1:ntr
        dur = sum(diff(tridx).tau);
        assert(dur<=TR(tridx),'diffusion gradients cannot be on for longer than TR!')
        diff(tridx).tau = [TR(tridx) - dur; diff(tridx).tau(:)];
        diff(tridx).G   = [0;               diff(tridx).G(:)];
    end
end

function [nshifts,dk] = computeshifts(diff)

    ntr = length(diff);

    % compute gradient moment for each TR
    G0 = zeros(1,ntr);
    for tridx=1:ntr
        G0(tridx) = dot(diff(tridx).G(:),diff(tridx).tau(:));
    end

    % confirm that gradient moments are all zero or an integer multiple of the smallest moment
    if any(G0~=0)
        deltaG0 = min(abs(G0));
        nshifts = G0/deltaG0;
        assert(all(nshifts>=0), 'negative gradient moments not implemented')

        % allow for small numerical imprecision
        assert(all(abs(nshifts-round(nshifts))<2*eps(G0)), 'gradient moments per TR are not all integer multiples of the shortest non-zero moment')
        nshifts = round(nshifts);
    else   
        deltaG0 = 0;
        nshifts = zeros(1,ntr);
    end

    % total dephasing between two EPG states
    gmT = 42.58e6 * 1e-3 * 2*pi; % rad s^-1 mT^-1
    dk = gmT*deltaG0*1e-3;
end