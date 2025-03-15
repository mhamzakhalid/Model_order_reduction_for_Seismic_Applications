function HFM = assemblesystem(HFM,var_para)

bd = [HFM.bdbot;HFM.bdleft;HFM.bdright]; % Dirichlet u=0 boundary condition on bottom left and right
HFM.ibcd = [2*bd - 1; 2*bd]; % indices for both horizontal and vertical displacement 
HFM.DOF_full = 2*HFM.nd; % Total degrees of freedom (DOFs) where mesh_data.nd = 2*numel(HFM.P(:,1))
free_dofs = setdiff(1:HFM.DOF_full,HFM.ibcd); % Internal DOFs excluding boundary conditions
HFM.DOF = numel(free_dofs); % Total number of internal DOFs
HFM.free_dofs = free_dofs;
%%%%% Discrete problem matrices 
currentfolder = pwd;
cd('FOM_DATA/SetupMatrices')
M = getMass(HFM,var_para.rhoi); % Matrix marix of size with rho values from var_param.rhoi
K = getStiffnessLay(HFM,var_para); % Stifness matrix including the values for the Lame parameters
M_rho0 = getMass(HFM,ones(5,1)); % Mass matrix with rho values as 1.
H = getH1norm_Matrix(HFM); % H - norm matrix for Korn's constant
[Ki_lam,Ki_mu] = getStiffnessAffineSplit(HFM); % % Parameter independent Cell array with 5 cells for each layer. Each matrix is mesh_data.DOF x mesh_data.DOF 

HFM.M_rho0 = M_rho0(free_dofs, free_dofs);
HFM.H = H(free_dofs,free_dofs);
HFM.K = K(free_dofs, free_dofs);
HFM.M = M(free_dofs, free_dofs);
 
for lid=1:5
    Ki_lam{lid}  =  Ki_lam{lid}(free_dofs,free_dofs);
    Ki_mu{lid}  =  Ki_mu{lid}(free_dofs,free_dofs);
end
HFM.Ki_lam = Ki_lam;
HFM.Ki_mu = Ki_mu;
HFM.Xh = HFM.M + HFM.K;
cd(currentfolder)
end
