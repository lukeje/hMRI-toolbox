function B1map = hmri_calc_AFI_B1map(S_TR1,S_TR2,TR2TR1ratio,nomFA)

% flip angle map in degrees
r=S_TR2./S_TR1;
n=TR2TR1ratio;
FAmap = acosd((r*n-1)./(n-r)); % Eq. (6) in Yarnykh, MRM (2007)

% relative B1 map in p.u.
B1map = 100*FAmap/nomFA;

end