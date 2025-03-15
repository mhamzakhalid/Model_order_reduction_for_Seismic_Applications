function Figure_11a_Plot(CaseROM,Pchange,Ntest,ROM_PODGreedy_m,ROM_Greedy_m,Nbasis)

fprintf('Running Case %d with %2.2f Change\n',CaseROM,Pchange)
CaseROM = string(CaseROM);

Pchange = Pchange/100;
%% Input Parameters 
addpath('FOM_DATA')
HFM = readmeshfiles();
var_para = getvariables(HFM);
% Assembling the discrete problem
HFM = assemblesystem(HFM,var_para);
HFM.Rec = getRecvector(HFM,var_para);
HFM.f = getSourcevector(HFM,var_para);
C_korn = computeKornsBetab(HFM);    
load('optimals0_Ricker.mat')

M_HFM = HFM.M; K_HFM = HFM.K; bi_HFM = HFM.f;Xh = HFM.Xh; DOF = HFM.DOF; 
Ki_HFMlam = HFM.Ki_lam; Ki_HFMmu = HFM.Ki_mu; Rec = HFM.Rec;
Nt = numel(var_para.t);
% Factorization of Xh = M + K 
[L_X,U_X,P_X,Q_X] =lu(Xh); 
% dualnorm of the output function (Receiver location)                        
[Output_Riesz] = Q_X*(U_X\(L_X\(P_X*Rec)));
ldualnorm = sqrt(abs(Rec'*Output_Riesz));
% dual norm of the Gaussian right hand side
source_Riesz=Q_X*(U_X\(L_X\(P_X*bi_HFM)));
fdualnorm = sqrt(abs(bi_HFM'*source_Riesz));

%% Weeks method time approximation
% Step 1: find z(sWeek_j)
% Step 2: Compute a_p = wi/Nz *\sum_{j=-Nz}^{N_z - 1} C_j^p*z(sWeek_j) where
% exp(-1i*p*theta_jhalf)/(1 - exp(1i*theta_jhalf) )
% Step 3: Compute z(t) = sum_{p=0}^{N_z-1} a_p exp( (s0(1) - s0(2))*t)*L_p(2*s0(2)*t)


var_para.width = 1.0*pi;
var_para.t0 = 4*pi/var_para.width;
var_para.s0 =s_optROM_Ricker; % s_optROM_Ricker; %[wr , wi] % Use [0.14 14.3] for longer than 20 final time (or recompute)
var_para.Nz = 608; % Number of frequency samples and Laguerre polynomials
wr = var_para.s0(1);
wi = var_para.s0(2);
jdx = -var_para.Nz:(var_para.Nz-1);
theta_jhalf = (jdx+1/2)*pi/var_para.Nz;
Cj = exp(1i*theta_jhalf);

sWeek = wr-wi*(Cj+1)./(Cj-1); % Same as wr + 1i*wi cot(theta_jhalf/2)
% Due to complex conjugation symmetry property we only need either + or - imaginary parts
sPstF = sWeek(var_para.Nz+1:end); % We take positive values here
smax = 11.7535; % This is obtained from the fixed parameter m case
[var_para.s_LFids] = find(abs(imag(sPstF(:)))<=smax);
s_LF = sPstF(var_para.s_LFids); % We consider only the low-frequencies
Ns_LF = numel(s_LF); % No of low-frequencies.
    
RICKLF = eval(var_para.source(var_para.Amp,var_para.width,var_para.t0,s_LF.')); %Ricker source at low-frequencies (The transpose in s_LF is only for the reasons concerning compatibility matrix/vector operations)

%% Training set for Pchange% change
minLam = abs(var_para.lami-Pchange*(var_para.lami));%
maxLam = abs(var_para.lami+Pchange*(var_para.lami));

minmu= abs(var_para.mui-Pchange*(var_para.mui));
maxmu= abs(var_para.mui+Pchange*(var_para.mui));

var_para.mui_min = minmu;
var_para.lam_min = minLam;
var_para.mui_max = maxmu;
var_para.lam_max = maxLam;

var_para.Ntest = Ntest;
%% Training set
% The solution for the following random parameters have already been computed

switch CaseROM
    case '1'
        rndvar = rand(var_para.Ntest,2); 
        lam_rand = (maxLam'-minLam').*rndvar(:,1) + minLam';
        mu_rand = (maxmu'-minmu').*rndvar(:,2) + minmu';
  case '2'
        % Every parameter goes independent change
         rndvar = rand(var_para.Ntest,5);
         
         lam_rand = ((maxLam'-minLam').*rndvar + minLam');
         mu_rand = ((maxmu'-minmu').*rndvar + minmu');
end
%% Computing seismograms using the full order model
 fprintf('Computing FOM solutions for error comparisons\n')
 addpath('Reduced_basis_methods/')

 Seismo_HFMTD = zeros(var_para.Ntest,Nt,2);
 
  parfor lam_id = 1:var_para.Ntest
      [Seismo_HFMTD(lam_id,:,:)] = getseismosTD_lam(M_HFM,bi_HFM,s_LF,RICKLF,var_para,Rec,[],lam_rand(lam_id,:),mu_rand(lam_id,:),Ki_HFMlam,Ki_HFMmu);
       fprintf('Process FOM solution computation:%2.2f\n',lam_id/Ntest*100)
  end
 fprintf('Done FOM solutions for error comparisons\n')

%% Computing errors

basis_ids = round(linspace(2,Nbasis,15));
Nb = numel(basis_ids);

   
%%%%%%%%%%%%%%%%
Aq = cell(11,1);
Aq{1} = M_HFM;
for qa = 2:6
    Aq{qa} = Ki_HFMlam{qa-1};    
    Aq{qa+5} = Ki_HFMmu{qa-1};      
end

%% POD-Greedy algorithm test
Error_seismoTDL2_PODGreedy = zeros(var_para.Ntest,Nb);


for Nk_id = 1:Nb
        Vk_it = ROM_PODGreedy_m.V(:,1:basis_ids(Nk_id));

       % 
       %  %%%%%%%%%%%%%% Project
         M_ROM = Vk_it'*(Aq{1}*Vk_it);
         bi_ROM =Vk_it'*(bi_HFM);
         K_ROMlam = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(2:6), 'UniformOutput', false);
         K_ROMmu = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(7:11), 'UniformOutput', false);

         %  %%%%%%%%%%%%%% Find RB seismograms
         Seismo_RBMTD_PODGreedy = zeros(var_para.Ntest,Nt,2);
       % Solving the ROM to find seismograms
         parfor lam_id = 1:var_para.Ntest
             [Seismo_RBMTD_PODGreedy(lam_id,:,:)] = getseismosTD_lam(M_ROM,bi_ROM,s_LF,RICKLF,var_para,Rec,Vk_it,lam_rand(lam_id,:),mu_rand(lam_id,:),K_ROMlam,K_ROMmu);
         end
          for lam_id = 1:var_para.Ntest
              % Compute error Laplace domain
              
              % Compute error time domain
              relative_seismoTD = sqrt( (norm(var_para.dt*Seismo_RBMTD_PODGreedy(lam_id,:,1),2))^2 + (norm(var_para.dt*Seismo_RBMTD_PODGreedy(lam_id,:,2),2))^2);
    
              er_FOM_ROMxTD =  Seismo_HFMTD(lam_id,:,1) - Seismo_RBMTD_PODGreedy(lam_id,:,1);
              er_FOM_ROMyTD =  Seismo_HFMTD(lam_id,:,2) - Seismo_RBMTD_PODGreedy(lam_id,:,2);
              Error_seismoTDL2_PODGreedy(lam_id,Nk_id) = (sqrt( (norm(var_para.dt*er_FOM_ROMxTD,2))^2 + (norm(var_para.dt*er_FOM_ROMyTD,2))^2))/...
                                                relative_seismoTD;          
         end 
        fprintf('Process %2.2f Nk: %d max Error %2.2e\n',Nk_id/Nb*100,basis_ids(Nk_id),max(Error_seismoTDL2_PODGreedy(:,Nk_id)))

end

fprintf('ROM Errors for POD-Greedy done\n')
%% Greedy algorithm test
if CaseROM == '1'
    Error_seismoTDL2_Greedy = zeros(var_para.Ntest,Nb);
        
    for Nk_id = 1:Nb
            Vk_it = ROM_Greedy_m.V(:,1:basis_ids(Nk_id));
           %  %%%%%%%%%%%%%% Project
             M_ROM = Vk_it'*(Aq{1}*Vk_it);
             bi_ROM =Vk_it'*(bi_HFM);
             K_ROMlam = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(2:6), 'UniformOutput', false);
             K_ROMmu = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(7:11), 'UniformOutput', false);
    
             %  %%%%%%%%%%%%%% Find RB seismograms
             Seismo_RBMTD_Greedy = zeros(var_para.Ntest,Nt,2);
           % Solving the ROM to find seismograms
             parfor lam_id = 1:var_para.Ntest
                 [Seismo_RBMTD_Greedy(lam_id,:,:)] = getseismosTD_lam(M_ROM,bi_ROM,s_LF,RICKLF,var_para,Rec,Vk_it,lam_rand(lam_id,:),mu_rand(lam_id,:),K_ROMlam,K_ROMmu);
             end
              for lam_id = 1:var_para.Ntest
                  % Compute error Laplace domain
                  
                  % Compute error time domain
                  relative_seismoTD = sqrt( (norm(var_para.dt*Seismo_RBMTD_Greedy(lam_id,:,1),2))^2 + (norm(var_para.dt*Seismo_RBMTD_Greedy(lam_id,:,2),2))^2);
        
                  er_FOM_ROMxTD =  Seismo_HFMTD(lam_id,:,1) - Seismo_RBMTD_Greedy(lam_id,:,1);
                  er_FOM_ROMyTD =  Seismo_HFMTD(lam_id,:,2) - Seismo_RBMTD_Greedy(lam_id,:,2);
                  Error_seismoTDL2_Greedy(lam_id,Nk_id) = (sqrt( (norm(var_para.dt*er_FOM_ROMxTD,2))^2 + (norm(var_para.dt*er_FOM_ROMyTD,2))^2))/...
                                                    relative_seismoTD;          
             end 
            fprintf('Process %2.2f Nk: %d max Error %2.2e\n',Nk_id/Nb*100,basis_ids(Nk_id),max(Error_seismoTDL2_Greedy(:,Nk_id)))
    
    end
    
    fprintf('DONE\n')
    
    save('Errors_reduction_m.mat','basis_ids','Error_seismoTDL2_PODGreedy','Error_seismoTDL2_Greedy')
    
    %%
    figure(1)
    set(gcf,'position',[200,550,550,500])
    
    errorbarlog(basis_ids,mean(Error_seismoTDL2_PODGreedy,1),std(Error_seismoTDL2_PODGreedy,1),':','linewidth',3,'Markersize',8)
    hold on
    errorbarlog(basis_ids,mean(Error_seismoTDL2_Greedy,1),std(Error_seismoTDL2_Greedy,1),':','linewidth',3,'Markersize',8)
    legend('POD$_s$-Greedy$_m$ algorithm','Greedy$_{(s,m)}$ algorithm','Interpreter','latex')
    xlabel('No. of RB functions: $N_k$','Interpreter','latex')
    ylabel('Mean$_{m\in\Xi_{t}}[\|\hat{e}_k(m)\|_{L^2}]_{rel}$','Interpreter','latex')
    hold off
    legend boxoff
    set(gca, 'FontName', 'Arial','FontSize',29)
    
    saveas(gcf,'Reduction_Test_PODGreedy_and_Greedy_Fig11a.png')
else

    figure(1)
    set(gcf,'position',[200,550,550,500])
    errorbarlog(basis_ids,mean(Error_seismoTDL2_PODGreedy,1),std(Error_seismoTDL2_PODGreedy,1),':','linewidth',3,'Markersize',8)
    legend('POD$_s$-Greedy$_m$ algorithm','Interpreter','latex')
    xlabel('No. of RB functions: $N_k$','Interpreter','latex')
    ylabel('Mean$_{m\in\Xi_{t}}[\|\hat{e}_k(m)\|_{L^2}]_{rel}$','Interpreter','latex')
    hold off
    legend boxoff
    set(gca, 'FontName', 'Arial','FontSize',29)
    
    saveas(gcf,'Reduction_Test_PODGreedy_Case2_Fig12.png')


    save('Errors_reduction_mCase2.mat','basis_ids','Error_seismoTDL2_PODGreedy')

end

end

%% Function for plotting
function hh = errorbarlog(varargin)
%ERRORBARLOG Symmetrical error bars for logarithmic Y-axis.
%   ERRORBARLOG(X,Y,E,...) plots the graph of vector X vs. vector Y with
%   a logarithmic Y-axis, using symmetrical bars about the data points, ie:
%   the bars are such that Y is the geometric mean (instead of arithmetic
%   mean) of the lower and upper bars. The total length of the error bar is
%   2E.
%
%   ERRORBARLOG has the same syntax as the original Matlab's ERRORBAR
%   function. The only difference is that while ERRORBAR displays the bars
%   symmetrically for a linear Y-axis (ie: Y is the arithmetic mean of
%   the lower and upper bars), ERRORBARLOG displays them symmetrically
%   for a logarithmic Y-axis.
%   
%   Example:
%      x=logspace(1,3,20);
%      y=5*(1 + 0.5*(rand(1,20)-0.5)).*x.^(-2);
%      errorbarlog(x,y,y/2,'o-');
%
%   F. Moisy
%   Revision: 1.01,  Date: 2006/09/08
%
%   See also ERRORBAR.
% History:
% 2005/05/28: v1.00, first version.
% 2006/09/08: v1.01, help text improved
error(nargchk(3,inf,nargin));
y=varargin{2};
e=varargin{3};
% computes the upper and lower error bars
% ymax and ymin are such that ymax*ymin=y^2 and ymax-ymin=2e.
% u is ymax-y and l is y-ymin.
ymax = e.*(1 + (1+(y./e).^2).^(1/2));
u = ymax - y;
l = 2*e + y - ymax;
h = errorbar(varargin{1:2},l,u,varargin{4:end});
set(gca,'YScale','log'); % set the Y axis in log coordinates.
if nargout>0,
    hh = h;
end
end
