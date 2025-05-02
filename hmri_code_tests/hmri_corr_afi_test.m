% test Nehrke AFI RF phase cycling matches standard phase cycling with pseudo-pulses
N = [1,3];
TR = 50*N;
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
TR = 50*N;
phi0 = 137;

npulse = 100;

offset = -phi0*(1-N(2)/N(1))/2;
pcN0 = RF_phase_cycle_Nehrke(npulse,phi0,N(1),N(2),offset);
pcNS = RF_phase_cycle_NehrkeSimplified(npulse,phi0,N(1),N(2));

assert(all(pcN0(:)-pcNS(:)<1e-14))