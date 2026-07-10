%{
m_dRFSpoilIncrement = 0.0; m_dRFSpoilMult = 0.0; m_dRFSpoilMult1 = 0.0; m_dRFSpoilMult2 = 0.0; m_dRFSpoilPhase = 0.0;

m_dRFSpoilIncrement = m_WIPParamTool.getDoubleValue(rMrProt, dWPT_RFSpoilIncrement)*(static_cast<double>(lTR2TR1ratio)+1.0)/2.0;


// Calculate additional phase needed for RF spoiling and add to FreqPhase
// RFSPOIL_INCREMENTdeg has a predefined value
// RFMAXPHASEdeg has a predefined value of 3600
if (rMrProt.fastImaging().getulEnableRFSpoiling() ) {
	m_dRFSpoilMult1 += static_cast<double>(lTR2TR1ratio);
	m_dRFSpoilMult2 += 1.0;
	if (crep == 0) {
		m_dRFSpoilMult += m_dRFSpoilMult1;
	} else {
		m_dRFSpoilMult += m_dRFSpoilMult2;
	}

	m_dRFSpoilPhase = fmod(m_dRFSpoilMult*m_dRFSpoilIncrement, 360.0);

	for(lEchoCounter=0; lEchoCounter<lMAX_ECHOES; lEchoCounter++) {
		m_sADC01zSet[lEchoCounter].increasePhase(m_dRFSpoilPhase);
		m_sADC01zNeg[lEchoCounter].decreasePhase(m_dRFSpoilPhase);
	}

	m_sSRFSinczSet.increasePhase(m_dRFSpoilPhase);
	m_sSRFSinczNeg.decreasePhase(m_dRFSpoilPhase);

	m_sSRFRectzSet.increasePhase(m_dRFSpoilPhase);
	m_sSRFRectzNeg.decreasePhase(m_dRFSpoilPhase);
}
%}

function phi = RF_phase_cycle_v1h(npulse,phi0,N)

RFSpoilIncrement = deg2rad(phi0*(N+1)/2);
RFSpoilMult = 0;
RFSpoilMult1 = 0;
RFSpoilMult2 = 0;

phi = zeros(npulse,1);
phase0 = 0;
for n=1:npulse
    RFSpoilMult1 = RFSpoilMult1 + N;
	RFSpoilMult2 = RFSpoilMult2 + 1;
    if mod(n,2)
		RFSpoilMult = RFSpoilMult + RFSpoilMult1;
	else
		RFSpoilMult = RFSpoilMult + RFSpoilMult2;
    end

	RFSpoilPhase = wrapTo2Pi(RFSpoilMult*RFSpoilIncrement);

    phase0 = phase0 + RFSpoilPhase;
    phi(n) = phase0;
end

end