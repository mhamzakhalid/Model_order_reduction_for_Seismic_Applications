%%%%%%%%%%%%

% This file describes the code to reproduce various figures in the manuscript.
% Some codes may run for several days depending on the computational architecture.
% Note: the execution must be done in the increasing order of
% CaseIDs
% of Case IDs. By running this file, the user can select from the available
% Case IDs. 
% -------------------------------------------------------------------------
% Figures and Instructions:
%
% Figure 1:
%   - Functionfile: run('Figure_1_Plot.m')
%   - Input: The precomputed file 'optimals0_Ricker.mat' (optimal parameters for the Laplace inversion via Weeks’ method)
%   - Output: Result_Fig1.png (Illustration of the results)
%   - CaseID: '1'
%
% Figure 4:
%   - Functionfile: Figure_4_Plot(alpha, Ntrain)
%   - Inputs:
%       alpha: A vector [1.0*pi, 1.5*pi, 2.0*pi] or a number (e.g., 1.0*pi)
%       Ntrain: 1024 (training set size)
%   - Outputs: ROM_POD_s.mat, ROM_SPOD_s.mat, ROM_Greedy_s.mat, 
%              Reduction_POD_and_Greedy_Fig4.png
%   - CaseID: '4'
%
% Figure 5:
%   - Functionfile: Figure_5_Plot(alpha, Ntest,ROM_POD,ROM_Greedy,ROM_SPOD)
%   - Inputs:
%       alpha: [1.0*pi, 1.5*pi, 2.0*pi]
%       Ntest: 512 Size of the test set
%       ROMs: ROM_POD,ROM_Greedy,ROM_SPOD from CaseID 4.
%   - Required Files: ROM_POD_s.mat, ROM_Greedy_s.mat

%   - Outputs: % Data files containing error values
%              ROM_POD_Test_FD.mat, 
%              ROM_Greedy_Test_FD.mat, 
%              ROM_Greedy_Test_TD.mat,
%              ROM_POD_Test_TD.mat, 
%              ROM_SPOD_Test_TD.mat, 
%              % Illutration of the results
%              Reduction_POD_and_Greedy_Fig4.png, (Figure 5 in the
%              manuscript)
%              Reduction_POD_and_SPOD_Supplementary.png (Supplementary
%              Figure SM1)
%   - CaseID: '5'
%
% Figure 8 (b and c) - POD over a Coarse Sample:
%   - Use: run('Figure_8_bc_Plot.m')
%   - Output: POD_Coarse.png (Illustration of results)
%   - CaseID: '8bc'
%
% Figure 9:
%   - Functionfile: Figure_9abc_12b_Plot(CaseROM, Pchange, Ntrain, Nkmax)
%   - Inputs:
%       For Case 1: CaseROM = 1, Pchange = 30, Ntrain = 512, Nkmax = 300 
%   - Outputs: ROM_PODGreedy_m.mat, ROM_Greedy_m.mat, Reduction_PODGreedy_and_Greedy_Fig9.png
%   - CaseID: '9'
%
% Figure 11a:
%   - Functionfile: Figure_11a_Plot(CaseROM, Pchange, Ntest,ROM_PODGreedy_m,ROM_Greedy_m, Nbasis)
%   - Inputs:
%       CaseROM = 1, Pchange = 30, Ntest = 128 
%       ROM_PODGreedy_m and ROM_Greedy_m (ROMs from Case 9)
%       Nbasis (maximum number of RB functions)
%   - Outputs: Errors_reduction_m.mat, Reduction_Test_PODGreedy_and_Greedy_Fig11a.png
%   - CaseID: '11a'
%
% Figure 11 (b):
%   - Functionfile: Figure_11b_Plot(ROM_PODGreedy_m,ROM_Greedy_m,Nb)
%   - Input: ROMs from Case 9 and Nb = maximum number of RB functions used
%   - Output: SeismogramHcompTD_PODGreedy_and_Greedy_Fig11b.png
%   - CaseID = '11b'
%
% Figure 12 (b) - Construct ROM for Case 2:
%   - Functionfile: Figure_9abc_12b_Plot(CaseROM, Pchange, Ntrain, Nkmax);
%   - Inputs: For Case 2: CaseROM = 2, Pchange = 5 or 10, Ntrain = 512, Nkmax = 450
%   - Outputs: ROM_PODGreedy_m_C2.mat, (Data ROM)
%              Reduction_PODGreedy_Fig12C2.png (Illustration)
%   - CaseID: '12b'
%
% Figure 12 (c) ~ - Test ROM for Case 2:
%   - Functionfile:Figure_12c_Plot(Pchange,ROM_PODGreedy_m_C2,Ntest);
%   - Input: Pchange = 10 (or 5), ROM_PODGreedy_m_C2 from CaseID 12b, Ntest = 128
%   - Outputs: Errors_reduction_m_Case2.mat, Reduction_Test_PODGreedy_Case2_Fig12c.png
%   - CaseID: '12c'
% 
%
% Additional Utilities:
%   - CaseID = '10' To compare inf-sup constants: getinfsupcomparison() outputs Inf_sup_comparison.png
%   - CaseID = 'wrwi' To compute optimal parameters: Compute_Optimal_Parameters() output: [wr, wi]
%   - To find C_W for the a posteriori error estimator: getCW() output: CW (a scalar)
%
% -------------------------------------------------------------------------
% How to Use:
%
%   1. Set the appropriate CaseID from below to choose which figure to reproduce.
%   2. The switch-case block will execute the corresponding MATLAB command.
%   3. Ensure that all required files are available in the working directory.
%

clear all
%parpool(32)

% Sample MATLAB Script to Switch Between Cases
% Possible options 
fprintf('Case IDs : Figure number\n');
fprintf('   1   : Figure 1\n');
fprintf('   4   : Figure 4 ROM construction fixed m\n');
fprintf('   5   : Figure 5 ROM testing fixed m\n');
fprintf('  8bc  : Figure 8 (b and c) (POD over coarse s)\n');
fprintf('   9   : Figure 9abc (ROM construction m change Case 1)\n');
fprintf('  11a  : Figure 11a (ROM test m change Case 1)\n');
fprintf('  11b  : Figure 11b (ROM seismogram comparison with Newmark-beta)\n');
fprintf(' 12b   :  Figure 12b (ROM construct for Case 2)\n');
fprintf(' 12c   : Figure 12c (ROM test for Case 2)\n');

% Prompt user for input; the 's' flag ensures the input is treated as a string
userInput = input('Enter the desired CaseID as e.g, 8bc: ', 's');
CaseID = userInput;  % This converts the input to a string if it isn't already

switch CaseID
    case '1'
        disp('Running Figure 1');
        % Reproduce Figure 1:
        % Precomputed file 'optimals0_Ricker.mat' is available.
        run('Figure_1_Plot.m');
        % Output: Result_Fig1.png

    case '4'
        disp('Running Figure 4');
        %%%% Inputs for Figure 4:
        alpha = [1.0*pi, 1.5*pi, 2.0*pi]; % or a number (e.g., alpha=1.0*pi)
        Ntrain = 1024; % Size of training set
        Figure_4_Plot(alpha, Ntrain) % Respective function file
        %%%% Outputs: For each alpha_id from 1 2 3 (depending on alpha)
        % ROM_POD_s.mat,  POD_s basis func.  as ROM_POD{alpha_id}.V 
        % ROM_SPOD_s.mat, Greedy_s basis func. as ROM_Greedy{alpha_id}.V 
        % ROM_Greedy_s.mat, SPOD basis func. as ROM_SPOD{alpha_id}.V 
        % Reduction_POD_and_Greedy_Fig4.png Figure 4
        % Some Additional Outputs in the respective .mat files include
        % i- Singular values ROM_POD{alpha_id}.Sigma
        % ii- Values of the error estimate ROM_Greedy{alpha_id}.maxDelta
        % iii- Values of the true error ROM_Greedy{alpha_id}.maxTrue

    case '5'
        disp('Running Figure 5');
        % Inputs for Figure 5:
        alpha = [1.0*pi, 1.5*pi, 2.0*pi];
        Ntrain = 512;

        % Ensure the following files are computed from CaseID 4
         load('ROM_POD_s.mat','ROM_POD')
         load('ROM_Greedy_s.mat','ROM_Greedy')
         load('ROM_SPOD_s.mat','ROM_SPOD')
                         
        Figure_5_Plot(alpha, Ntrain,ROM_POD,ROM_Greedy,ROM_SPOD)

        % Outputs: Frequency domain as FD and time domain as TD
        %'ROM_POD_Test_FD.mat' contains 'ErrorX_rel_POD','ErrorX_Omega_POD','Test_set'
        %'ROM_Greedy_Test_FD.mat' contains'ErrorX_rel_Greedy','ErrorX_Omega_Greedy','Test_set' 
        %%% TD Errors
        %'ROM_Greedy_Test_TD.mat' contains'error_seismo_Greedy' (TD error)
        %'ROM_POD_Test_TD.mat' contains 'error_seismo_POD' (TD error)
        %'ROM_SPOD_Test_TD.mat' contains'error_seismo_SPOD' (TD error)
        %%% Illustration of results
        % Reduction_POD_and_Greedy_Fig5.png (Figure 5 in the manuscript)
        % Reduction_POD_and_SPOD_Supplementary.png (Fig. SM1)

    case '8bc'
        disp('Running Figure 8 (b and c)');
        % Reproduce Figure 7 b and c (POD over a coarse sample)
        run('Figure_8_bc_Plot.m');
        % Output: POD_Coarse.png (Illustration of results)

    case '9'
        disp('Running Figure 9abc to construct ROMs for variations in m');
        % Inputs for Figure 9abc:
        CaseROM = 1; % For the First Case
        Pchange = 30; % Percentage change
        Ntrain = 512; % Size of the training set
        Nkmax = 300; % Maximum number of basis functions
        Figure_9abc_12b_Plot(CaseROM, Pchange, Ntrain, Nkmax); 
        
        % Outputs: 
        % ROM_PODGreedy_m.mat (RB as ROM_PODGreedy_m.V)
        % ROM_Greedy_m.mat,  (RB as ROM_Greedy_m.V)
        % Illustration of results
        % Reduction_PODGreedy_and_Greedy_Fig9.png

    case '11a'
        disp('Running Figure 11a testing ROMs');
        % Inputs for Figure 11a:
        CaseROM = 1; Pchange = 30; Ntest = 128; 
        % Ensure the following files exists from Case 9
        load('ROM_PODGreedy_m.mat','ROM_PODGreedy_m')
        load('ROM_Greedy_m.mat','ROM_Greedy_m')

        Nbasis = size(ROM_PODGreedy_m.V,2); % Maximum number of RB functions

        Figure_11a_Plot(CaseROM,Pchange,Ntest,ROM_PODGreedy_m,ROM_Greedy_m,Nbasis)
        % Outputs: 
        % Errors_reduction_m.mat, 
        % Reduction_Test_PODGreedy_and_Greedy_Fig11a.png

    case '11b'
        disp('Running Figure 11 (b) Case 1 seismograms plot');
        % First compute ROM_PODGreedy_m and ROM_Greedy_m by using Case ID 9
        load('ROM_PODGreedy_m.mat','ROM_PODGreedy_m')
        load('ROM_Greedy_m.mat','ROM_Greedy_m')
        
        Nbasis = size(ROM_PODGreedy_m.V,2); % Maximum number of RB functions

        Figure_11b_Plot(ROM_PODGreedy_m,ROM_Greedy_m,Nbasis);
        
        % Output: SeismogramHcompTD_PODGreedy_and_Greedy_Fig11b.png
        % (Horizontal component)
        % SeismogramVcompTD_PODGreedy_and_Greedy_Fig11b.png (Vertical
        % component)

    case '12b'
        disp('Running Figure 12b ROM construction for Case 2...');
        % Inputs for Figure 12b:
        CaseROM = 2; 
        Pchange = 10;  % Alternatively, Pchange can be 5
        Ntrain = 512; 
        Nkmax = 450;
        Figure_9abc_12b_Plot(CaseROM, Pchange, Ntrain, Nkmax);
        % Outputs: ROM_PODGreedy_m_C2.mat, (ROM for Case 2)
        % Reduction_PODGreedy_Fig12C2.png (Illustration of results)

    case '12c'
        disp('Running Figure 12c ROM test for Case 2...');
        % Input for Figure 10c:
        Pchange = 10;  % Alternatively, Pchange can be 5
        Ntest = 128;
        load('ROM_PODGreedy_m_C2.mat','ROM_PODGreedy_m_C2') % First compute ROM_PODGreedy_m_C2 in Case 12b
        
        Figure_12c_Plot(Pchange,ROM_PODGreedy_m_C2,Ntest);
        % Outputs: Errors_reduction_m_Case2.mat, (Data)
        %          Reduction_Test_PODGreedy_Case2_Fig12c.png (Illustration of results)
    case '10'
        getinfsupcomparison()
        %Output: Inf_sup_comparison.png Illustration
    case 'wrwi'
        Compute_Optimal_Parameters()
        % Output [wr wi]
    otherwise
        error('Undefined CaseID. Please check the input and select a valid case.');
end
