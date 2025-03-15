function Figure_8_Plot(CaseROM,Pchange,Ntrain,Nkmax)    
% Pchange = 0.30
% CaseROM = 'Case 1'
%parpool(8)
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
%% Plotting
blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
cl_colorsdark = {blue, red, black, ...
             green, brown, purple};
blue = [114 147 203]./255;
red = [211 94 96]./255;
black = [128 133 133]./255;
green = [132 186 91]./255;
brown = [171 104 87]./255;
purple = [144 103 167]./255;
cl_colorslight = {blue, red, black, ...
             green, brown, purple};
scalingEstimator = 1;
scalingUpperbound = 1;

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

%% Training set for P% change
minLam = abs(var_para.lami-Pchange*(var_para.lami));%
maxLam = abs(var_para.lami+Pchange*(var_para.lami));

minmu= abs(var_para.mui-Pchange*(var_para.mui));
maxmu= abs(var_para.mui+Pchange*(var_para.mui));

var_para.mui_min = minmu;
var_para.lam_min = minLam;
var_para.mui_max = maxmu;
var_para.lam_max = maxLam;


var_para.Ntrain = Ntrain;
%% Training set
% The solution for the following random parameters have already been computed

switch CaseROM
    case '1'
        rndvar = rand(var_para.Ntrain,2); 
        lam_rand = (maxLam'-minLam').*rndvar(:,1) + minLam';
        mu_rand = (maxmu'-minmu').*rndvar(:,2) + minmu';
  case '2'
        % Every parameter goes independent change
         rndvar = rand(var_para.Ntrain,5);
         
         lam_rand = ((maxLam'-minLam').*rndvar + minLam');
         mu_rand = ((maxmu'-minmu').*rndvar + minmu');
end
%% Computing seismograms using the full order model
 fprintf('Computing FOM solutions for error comparisons\n')
 addpath('Reduced_basis_methods/')

 Seismo_HFMTD = zeros(var_para.Ntrain,Nt,2);
 Seismo_HFMLD = zeros(var_para.Ntrain,Ns_LF,2);
 
  parfor lam_id = 1:var_para.Ntrain
      [Seismo_HFMTD(lam_id,:,:),Seismo_HFMLD(lam_id,:,:)] = getseismosTD_lam(M_HFM,bi_HFM,s_LF,RICKLF,var_para,Rec,[],lam_rand(lam_id,:),mu_rand(lam_id,:),Ki_HFMlam,Ki_HFMmu);
       fprintf('Process FOM solution computation:%2.2f\n',lam_id/Ntrain*100)
  end
 fprintf('Done FOM solutions for error comparisons\n')


%% Running the POD-Greedy algorithm
    %%%%%%%%%%%%%%%% Tools for the a posteriori error estimator
    %%% Check book Page no. 62, Alg 3.4 'Reduced Basis Methods for Partial Differential
    %%% Equations: An Introduction' By Alfio Quarteroni et. al. 
    % source_Riesz = Q_X*(U_X\(L_X\(P_X*Gauss_f)));%Riesz_representer_Xnorm(b_HFM, L, L', perm);%FOM.Xnorm\FOM.Fq{q1};    
    Cqq = source_Riesz'*bi_HFM; 
    THETA_F_Vec = RICKLF.';
    estimator.res_ff = conj(THETA_F_Vec).*THETA_F_Vec*Cqq;
    estimator.THETA_F_Vec = THETA_F_Vec;

    estimator.L_X = L_X ;
    estimator.U_X = U_X;
    estimator.P_X = P_X ;
    estimator.Q_X = Q_X;
    %
    estimator.Cj = abs(1./(1 - exp(1i*theta_jhalf)));
    estimator.CW = 21.57;%(var_para.s0(2))/(var_para.Nz)*sqrt(sum(abs(Lag_intC2),'all'))

    % Lower bound on the inf-sup
    c1 = min([1,1./max(minLam(:)./var_para.lami(:)),1./max(minmu(:)./var_para.mui(:))]);
    c2 = max([1,max(maxLam(:)./var_para.lami(:)),max(maxmu(:)./var_para.mui(:))]);
    c_equivalance = c1/c2;
    estimator.beta_lb = c_equivalance*min([1,wr^2,2*wr/sqrt(4*wr^2+smax^2)]);
    estimator.ldualnorm = ldualnorm;
%%
    
    %%%%%%%%%%%%%%%%
    Aq = cell(11,1);
    Aq{1} = M_HFM;
    for qa = 2:6
        Aq{qa} = Ki_HFMlam{qa-1};    
        Aq{qa+5} = Ki_HFMmu{qa-1};      
    end
var_para.Nkmax = Nkmax;
 fprintf('ROM generation using POD-Greedy algorithm\n')

switch CaseROM
    case '1'
         ROM_PODGreedy_m = getPODGreedyalgorithm(Aq,bi_HFM,Xh,Rec,s_LF,RICKLF,lam_rand,mu_rand,estimator,var_para,Seismo_HFMTD,Seismo_HFMLD)
        
         fprintf('Done generation using POD-Greedy algorithm\n')
        
        save('ROM_PODGreedy_m.mat','ROM_PODGreedy_m','-v7.3')
         
        fprintf('ROM generation using Greedy algorithm\n')
        
        
        ROM_Greedy_m = getGreedyalgorithm_m(Aq,bi_HFM,Xh,Rec,s_LF,RICKLF,lam_rand,mu_rand,estimator,var_para,Seismo_HFMTD,Seismo_HFMLD)
        save('ROM_Greedy_m.mat','ROM_Greedy_m','-v7.3')
        
        fprintf('Done generation using Greedy algorithm\n')

                %%
        % Figure 8 (a) POD-Greedy Algorithm for 30% Change

        figure(1)
        bids = 10:var_para.Nkmax;
        mark_id = round(linspace(1,var_para.Nkmax-10,15));
        subplot(1,3,1)
        semilogy(bids,ROM_PODGreedy_m.maxEstimatorm(bids)*scalingEstimator,'-v','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        hold on
        semilogy(bids,ROM_PODGreedy_m.maxTruem(bids)*scalingUpperbound,'-.>','Color',cl_colorsdark{2},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        semilogy(bids,ROM_PODGreedy_m.maxTrueL2m(bids),':<','Color',cl_colorsdark{3},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        hold off
        set(gca, 'FontName', 'Arial','FontSize',29)
        hold off
        xlabel('No. of RB functions: $N_k$','Interpreter','latex')
        legend('$e_k(m)=\Delta_{k,t}(m)$','$e_k(m)=\widetilde{\Delta}_{k,t}(m)$',...
                '$e_k(m)=\|\hat{e}_k(m)\|_{L^2}$','Interpreter','latex','Location','southwest')
        legend boxoff
        ylabel('$\max_{m\in\Xi_m} [e_k(m)]_{rel}$','Interpreter','latex')
        ylim([1e-4 1e7])
        
        % Figure 8 (a) Effectivity
        subplot(1,3,2)
        
        semilogy(bids,ROM_PODGreedy_m.maxEstimatorm(bids)*scalingEstimator./(ROM_PODGreedy_m.maxTruem(bids)*scalingUpperbound),'-v','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        
        hold on
        semilogy(bids,ROM_PODGreedy_m.maxEstimatorm(bids)*scalingEstimator./ROM_PODGreedy_m.maxTrueL2m(bids),'->','Color',cl_colorsdark{2},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        hold off
        set(gca, 'FontName', 'Arial','FontSize',29)
        ylim([1e-4 1e7])
        hold off
        xlabel('No. of RB functions: $N_k$','Interpreter','latex')
        legend('$\frac{\max_{m\in\Xi_m} \Delta_{k,t}(m)}{\max_{m\in\Xi_m} \widetilde{\Delta}_{k,t}(m)}$',...
                '$\frac{\max_{m\in\Xi_m} \Delta_{k,t}(m)}{\max_{m\in\Xi_m} \|\hat{e}_k(m)\|_{L^2} }$','Interpreter','latex','Location','southwest')
        legend boxoff
        ylabel('Effectivity of $\Delta_{k,t}(m)$','Interpreter','latex')
        
        % Figure 8 (c) Greedy Algorithm for 30% Change
        subplot(1,3,3)
        
        semilogy(bids,ROM_Greedy_m.maxEstimatorm(bids)*scalingEstimator,'-v','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        hold on
        semilogy(bids,ROM_Greedy_m.maxTruem(bids)*scalingUpperbound,'-.>','Color',cl_colorsdark{2},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        semilogy(bids,ROM_Greedy_m.maxTrueL2m(bids),':<','Color',cl_colorsdark{3},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        hold off
        set(gca, 'FontName', 'Arial','FontSize',29)
        hold off
        xlabel('No. of RB functions: $N_k$','Interpreter','latex')
        legend('$e_k(m)=\Delta_{k,t}(m)$','$e_k(m)=\widetilde{\Delta}_{k,t}(m)$',...
                '$e_k(m)=\|\hat{e}_k(m)\|_{L^2}$','Interpreter','latex','Location','southwest')
        legend boxoff
        ylabel('$\max_{m\in\Xi_m} [e_k(m)]_{rel}$','Interpreter','latex')
        ylim([1e-4 1e7])
        
        set(gcf,'position',[200,550,1550,500])
        saveas(gcf,'Reduction_PODGreedy_and_Greedy_Fig9.png')
        close all


    case '2'
        ROM_PODGreedy_m_C2 = getPODGreedyalgorithm(Aq,bi_HFM,Xh,Rec,s_LF,RICKLF,lam_rand,mu_rand,estimator,var_para,Seismo_HFMTD,Seismo_HFMLD)
        save('ROM_PODGreedy_m_C2.mat','ROM_PODGreedy_m_C2','-v7.3')

        fprintf('Done generation using POD-Greedy algorithm\n')
        figure(1)
        bids = 10:var_para.Nkmax;
        mark_id = round(linspace(1,var_para.Nkmax-10,15));
        semilogy(bids,ROM_PODGreedy_m_C2.maxEstimatorm(bids)*scalingEstimator,'-v','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        hold on
        semilogy(bids,ROM_PODGreedy_m_C2.maxTrueL2m(bids),':<','Color',cl_colorsdark{3},'LineWidth',3,'Markersize',12,'Markerindices',mark_id)
        hold off
        set(gca, 'FontName', 'Arial','FontSize',29)
        hold off
        xlabel('No. of RB functions: $N_k$','Interpreter','latex')
        legend('$e_k(m)=\Delta_{k,t}(m)$',...
                '$e_k(m)=\|\hat{e}_k(m)\|_{L^2}$','Interpreter','latex','Location','southwest')
        legend boxoff
        ylabel('$\max_{m\in\Xi_m} [e_k(m)]_{rel}$','Interpreter','latex')
        ylim([1e-4 1e7])
                
        saveas(gcf,'Reduction_PODGreedy_Fig12C2.png')

end
rmpath('FOM_DATA')
rmpath('Reduced_basis_methods/')





end
