# MOR_For_Seismic_Applications_Script_Files

[![DOI](https://zenodo.org/badge/811891525.svg)](https://zenodo.org/doi/10.5281/zenodo.11520031)
------
README
------

This archive contains the source code to reproduce all numerical results of the publication

"Model order reduction for seismic applications"
by R. Hawkins, M.H. Khalid, K. Smetana, M. Schlottbom.



All scripts are tested with MATLAB 2023b.
ATTENTION: Most of the code is extremely time consuming as, for instance, the reduction in the m-parameters requires very many full order model computations to compare errors.
(The code will probably run for several days depending on the computational architecture.)

## To reproduce Figure 1
%   - Command: run('Figure_1_Plot.m')

%   - Here, 'optimals0_Ricker.mat' has been pre-computed optimal parameters that we use for the numerical Laplace inversion via Weeks’ method

%   - Output: Result_Fig1.png (Illustration of the results)

%   - CaseID: '1' in the main script

## To reproduce Figure 4
Figure 4:
%   - Functionfile: Figure_4_Plot(alpha, Ntrain)

%   - Inputs:

%       alpha: A vector [1.0*pi, 1.5*pi, 2.0*pi] or a number (e.g., 1.0*pi)

%       Ntrain: 1024 (training set size)

%   - Outputs: ROM_POD_s.mat, ROM_SPOD_s.mat, ROM_Greedy_s.mat, 

%              Reduction_POD_and_Greedy_Fig4.png

%   - CaseID: '4'in the main script

## To reproduce Figure 5
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

Additional Description: 

i- ROM_POD_Test_FD.mat (This contains frequency domain errors for the POD over a uniform random test set different from the training set. The .mat file contains 'ErrorX_rel_POD', 'ErrorX_Omega_POD', 'Test_set' ) 
ii -ROM_Greedy_Test_FD.mat (This contains frequency domain errors for the greedy algorithm over a uniform random test set that is different from the training set. The .mat file contains 'ErrorX_rel_Greedy','ErrorX_Omega_Greedy','Test_set') 
iii- ROM_Greedy_Test_TD.mat (This contains time domain errors for the greedy algorithm over 5 uniform receivers. The .mat file contains 'ErrorX_rel_Greedy','ErrorX_Omega_Greedy','Test_set') 
iv- ROM_Greedy_Test_TD.mat (This contains time domain errors for the greedy algorithm over 5 uniform receivers. The .mat file contains 'error_seismo_POD') 
v - ROM_POD_Test_TD.mat (This contains time domain errors for the POD over 5 uniform receivers. The .mat file contains 'error_seismo_POD') 
vi - ROM_SPOD_Test_TD.mat ( contains error_seismo_SPOD )



# To reproduce Figure 8 b and c (POD over a coarse sample)
Use: run('Figure_8_bc_Plot.m')
Output:

POD_Coarse.png (Illustration of results)
- CaseID: '8bc'
  
% Note: The codes for the following cases can take several days to complete. It is recommended to use many cores.
## To reproduce Figure 9

Use: Figure_9abc_12b_Plot(CaseROM,Pchange,Ntrain,Nkmax)

#Inputs For Case 1
CaseROM =  1 (Case as in the paper)
Pchange = 30  (Percentage Change) 
Ntrain = 512 (No of training parameters)
Nkmax = 300 (Maximum number of RB functions)
#Outputs for Case 1
i- 'ROM_PODGreedy_m.mat' (ROM for POD-Greedy algorithm, name file is ROM_PODGreedy_m that is a struct containing maximum values of estimator/true_error/upperbound)
ii- ROM_Greedy_m.mat', (ROM for POD-Greedy algorithm, name file is 'ROM_Greedy_m that is a struct containing maximum values of estimator/true_error/upperbound)
iii - Reduction_PODGreedy_and_Greedy_Fig9.png (Illustration of the results for the effectivity and true error)
- Case ID = 9
  
## To reproduce Figure 11a
Use : Figure_11a_Plot(CaseROM, Pchange, Ntest,ROM_PODGreedy_m,ROM_Greedy_m, Nbasis)
#Required files from Case ID 9 
ROM_PODGreedy_m.mat
ROM_Greedy_m.mat
Inputs:
CaseROM = 1; (Case 1)
Pchange = 30; (Percentage change)
Ntest = 128; (No of random parameters for testing)
Outputs:
Errors_reduction_m.mat (L2 Errors named as 'Error_seismoTDL2_PODGreedy' and 'Error_seismoTDL2_Greedy')
Reduction_Test_PODGreedy_and_Greedy_Fig9a.png (Illustration of the maximum error convergence)

## To reproduce Figure 9b and 9c
Use: Figure_9bc_Plot()
#Required files
ROM_PODGreedy_m.mat
ROM_Greedy_m.mat
Output:
SeismogramHcompTD_PODGreedy_and_Greedy_Fig9bc.png (Illustration of the seismograms and errors)

## To reproduce Figure 10 b (For Case 2)
Use: run('Figure_8_Plot(CaseROM,Pchange,Ntrain,Nkmax)')
Inputs: 
CaseROM = 2;
Pchange = 5; or Pchange = 10;
Ntrain = 512;
Nkmax = 450;
Output: 
ROM_PODGreedy_m_C2.mat (The reduced order models named as 'ROM_PODGreedy_m_C2')
Reduction_PODGreedy_Fig10C2.png (Illustration of the ROM maximum error convergence)

## To reproduce Figure 10 c (Test ROM for Case 2)
Use: run('Figure_10c_Plot(Pchange)')
#Required files
ROM_PODGreedy_m_C2.mat
Input: 
Pchange = 10 (Or Pchange = 5)
Output: 
Errors_reduction_m_Case2.mat,(Error file containing 'basis_ids' and 'Error_seismoTDL2_PODGreedy')
Reduction_Test_PODGreedy_Case2_Fig10c.png (Illustration of mean error convergence)

## To find the optimal parameters use
Use: Compute_Optimal_Parameters()
Output
s0 = [wr wi] (Vector containg optimal parameters)

## To find C_W for the a posteriori error estimator
Use: getCW()
Output:
CW (scalar)

## To compare inf-sup constants

Use: getinfsupcomparison()

Output:
Inf_sup_comparison.png 


For queries: Contact: mhamzakhalid15@gmail.com
