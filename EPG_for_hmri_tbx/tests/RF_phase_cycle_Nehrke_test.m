classdef RF_phase_cycle_Nehrke_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.
        ndim = {2,3};
        fraction = {0,0.3,0.5,0.7,1.0};
    end
    
    methods (Test)
        function NehrkeTest(testCase)
            npulse = 20;
            N1 = 3;
            N2 = 5;
            phi0 = 129.3;
            phi_Nehrke = wrapTo2Pi(RF_phase_cycle_Nehrke(npulse, phi0, N1, N2));
            phi_standard = wrapTo2Pi(RF_phase_cycle(ceil(npulse/2)*N1+floor(npulse/2)*N2, phi0))';
            testCase.assertEqual(phi_Nehrke(1:2:end), phi_standard(1:(N1+N2):end), 'AbsTol',1e-11) 
            testCase.assertEqual(phi_Nehrke(2:2:end), phi_standard(1+N1:(N1+N2):end), 'AbsTol',1e-11) 
        end


        function NehrkeSimplifiedTest(testCase)
            npulse = 20;
            N1 = 3;
            N2 = 5;
            phi0 = 129.3;
            offset = -phi0*(1-N2)/2;
            phi_Nehrke = wrapTo2Pi(RF_phase_cycle_Nehrke(npulse, phi0, N1, N2, offset));
            phi_simplified = wrapTo2Pi(RF_phase_cycle_NehrkeSimplified(npulse, phi0, N1, N2));
            testCase.assertEqual(phi_simplified, phi_Nehrke, 'AbsTol',1e-11) 
        end

        %{
        function NehrkeV1hTest(testCase)
            npulse = 20;
            N1 = 1;
            N2 = 5;
            phi0 = 129.3;
            phi_v1h = wrapTo2Pi(RF_phase_cycle_v1h(npulse, phi0, N2/N1));
            phi_simplified = wrapTo2Pi(RF_phase_cycle_NehrkeSimplified(npulse, phi0, N1, N2));
            testCase.assertEqual(phi_v1h, phi_simplified, 'AbsTol',1e-11) 
        end

        function NehrkeErrorTest(testCase)
            npulse = 20;
            N1 = 1;
            N2 = 5;
            phi0 = 129.3;
            phi_error = wrapTo2Pi(RF_phase_cycle_NehrkeSimplifiedError(npulse, phi0, [N1/N2,1]));
            phi_simplified = wrapTo2Pi(phi_error(1)+RF_phase_cycle_NehrkeSimplified(npulse, phi0/N2, N1, N2));
            testCase.assertEqual(phi_error, phi_simplified, 'AbsTol',1e-11) 
        end
        %}
    end
end