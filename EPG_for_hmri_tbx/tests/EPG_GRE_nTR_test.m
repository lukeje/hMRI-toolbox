
classdef EPG_GRE_nTR_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.
        ndim = {2,3};
    end
    
    methods (Test)
        
        %% Test Functions
        function NDgradaxisTest(testCase,ndim)
            theta = deg2rad(repmat(30,100,1));
            phi   = zeros(size(theta));
            TR    = 100; % [ms]
                
            T1 = 1200;             % [ms]
            T2 = 30;               % [ms]
            D  = 0.7;              % [µm^2/ms]

            dur1 = 55; % ms
            Gdur = [3,dur1/4,dur1/2,dur1/4]; % [ms]
            Gamp = [26,30,-30,30];           % [mT/m]
            GdiffRef = struct('D', D*1e-9, 'G', Gamp, 'tau', Gdur);

            GdiffAxes = repmat({struct('D', D*1e-9, 'G', Gamp, 'tau', Gdur)},ndim,1);
            for n=1:ndim
                GdiffAxes{n}.G = GdiffRef.G/sqrt(ndim);
            end

            naxis = EPG_GRE_nTR(theta,phi,TR,T1,T2,'diff', GdiffAxes);
            ref   = EPG_GRE_nTR(theta,phi,TR,T1,T2,'diff', GdiffRef);

            assertEqual(testCase, ref, naxis, 'AbsTol',1e-3);
        end
        
    end
 
end