
classdef EPG_GRE_nTR_test < matlab.unittest.TestCase
    properties (TestParameter)
        % Augment TestParameter with parameters over which tests will run,
        % as well as parameters needed by the test functions.
        ndim = {2,3};
        fraction = {0,0.3,0.5,0.7,1.0};
    end
    
    methods (Test)
        
        %% Test Functions
        function NDgradaxisTest(testCase,ndim)
            theta = deg2rad(repmat(30,100,1));
            phi   = RF_phase_cycle(length(theta),137);
            TR    = 100; % [ms]
                
            T1 = 1200; % [ms]
            T2 = 30;   % [ms]
            D  = 0.7;  % [µm^2/ms]

            dur = 55; % ms
            Gdur = [3,dur/4,dur/2,dur/4]; % [ms]
            Gamp = [26,30,-30,30];        % [mT/m]
            GdiffRef = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);

            GdiffAxes = repmat({struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur)},ndim,1);
            for n=1:ndim
                GdiffAxes{n}.G = GdiffRef.G/sqrt(ndim);
            end

            naxis = EPG_GRE_nTR(theta,phi,TR,T1,T2, 'diff',GdiffAxes);
            ref   = EPG_GRE_nTR(theta,phi,TR,T1,T2, 'diff',GdiffRef);

            assertEqual(testCase, naxis, ref, 'AbsTol',1e-12);
        end

        function twogradaxisfractionsTest(testCase,fraction)
            theta = deg2rad(repmat(30,200,1));
            phi   = RF_phase_cycle(length(theta),137);
            TR    = 100; % [ms]
                
            T1 = 1200; % [ms]
            T2 = 30;   % [ms]
            D  = 0.7;  % [µm^2/ms]

            dur1 = 55; % ms
            Gdur = [3,dur1/4,dur1/2,dur1/4]; % [ms]
            Gamp = [26,30,-30,30];           % [mT/m]
            GdiffRef = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);

            GdiffAxes{1} = GdiffRef;
            GdiffAxes{1}.G = GdiffRef.G*fraction;

            GdiffAxes{2} = GdiffRef;
            GdiffAxes{2}.G = GdiffRef.G*sqrt(1-fraction^2);

            naxis = EPG_GRE_nTR(theta,phi,TR,T1,T2, 'diff',GdiffAxes);
            ref   = EPG_GRE_nTR(theta,phi,TR,T1,T2, 'diff',GdiffRef);

            assertEqual(testCase, naxis, ref, 'AbsTol',1e-12);
        end

        function EPG_GRE_2TRTest(testCase)
            theta = deg2rad(repmat(30,200,1));
            phi   = RF_phase_cycle(length(theta),137);
            TR    = [20,20]; % [ms]
                
            T1 = 1200; % [ms]
            T2 = 30;   % [ms]
            D  = 0.7;  % [µm^2/ms]

            dur = 10; % ms
            Gdur = [3,dur/4,dur/2,dur/4]; % [ms]
            Gamp = [26,30,-30,30];           % [mT/m]
            Gdiff(1) = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);
            Gdiff(2) = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);

            naxis = EPG_GRE_nTR(theta,phi,TR(1),T1,T2, 'diff',Gdiff(1));
            ref   = EPG_GRE_nTR(theta,phi,TR,   T1,T2, 'diff',Gdiff);

            assertEqual(testCase, naxis, ref);
        end

        function EPG_GRE_nTRvsEPG_GRETest(testCase)
            theta = deg2rad(repmat(30,200,1));
            phi   = RF_phase_cycle(length(theta),137);
            TR    = 20; % [ms]
  
            T1 = 1200; % [ms]
            T2 = 30;   % [ms]
            D  = 0.7;  % [µm^2/ms]

            dur = 10; % ms
            Gdur = [3,dur/4,dur/2,dur/4]; % [ms]
            Gamp = [26,30,-30,30];        % [mT/m]
            Gdur(end+1) = TR-sum(Gdur); % ensure that we fill the TR
            Gamp(end+1) = 0;
            Gdiff(1) = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);

            naxis = EPG_GRE_nTR(theta,phi,TR,T1,T2, 'diff',Gdiff);
            ref   = EPG_GRE(    theta,phi,TR,T1,T2, 'diff',Gdiff);

            assertEqual(testCase, naxis, ref);
        end

        function EPG_GRE_2TRvsEPG_GRETest(testCase)
            theta = deg2rad(repmat(30,200,1));
            phi   = RF_phase_cycle(length(theta),137);
            TR    = [20,20]; % [ms]
  
            T1 = 1200; % [ms]
            T2 = 30;   % [ms]
            D  = 0.7;  % [µm^2/ms]

            dur = 10; % ms
            Gdur = [3,dur/4,dur/2,dur/4]; % [ms]
            Gamp = [26,30,-30,30];        % [mT/m]
            Gamp = [Gamp,0];
            Gdiff(1) = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);
            Gdiff(1).tau = [Gdiff(1).tau,TR(1)-sum(Gdiff(1).tau)];
            Gdiff(2) = Gdiff(1);

            naxis = EPG_GRE_nTR(theta,phi,TR,   T1,T2, 'diff',Gdiff);
            ref   = EPG_GRE(    theta,phi,TR(1),T1,T2, 'diff',Gdiff(1));

            assertEqual(testCase, naxis, ref);
        end

        function twoTRtwogradaxisTest(testCase)
            theta = deg2rad(repmat(30,200,1));
            phi   = RF_phase_cycle(length(theta),137);
            TR    = [100,50]; % [ms]
                
            T1 = 1200; % [ms]
            T2 = 30;   % [ms]
            D  = 0.7;  % [µm^2/ms]

            Gdur = 3;  % [ms]
            Gamp = 26; % [mT/m]
            GdiffRef(1) = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);
            GdiffRef(2) = struct('D',D*1e-9, 'G',Gamp, 'tau',Gdur);

            f = 0.1;
            GdiffAxes{1} = GdiffRef;
            GdiffAxes{1}(1).G = GdiffRef(1).G*f;
            GdiffAxes{1}(2).G = GdiffRef(2).G*f;

            GdiffAxes{2} = GdiffRef;
            GdiffAxes{2}(1).G = GdiffRef(1).G*sqrt(1-f^2);
            GdiffAxes{2}(2).G = GdiffRef(2).G*sqrt(1-f^2);

            naxis = EPG_GRE_nTR(theta,phi,TR,T1,T2, 'diff',GdiffAxes);
            ref   = EPG_GRE_nTR(theta,phi,TR,T1,T2, 'diff',GdiffRef);

            assertEqual(testCase, naxis, ref, 'AbsTol',1e-12);
        end
        
    end
 
end