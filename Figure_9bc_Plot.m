%%% This file contains the script file to test reduced order models
function Figure_9bc_Plot()

%% Input Parameters 
addpath('FOM_DATA')
HFM = readmeshfiles();
var_para = getvariables(HFM);
% Assembling the discrete problem
HFM = assemblesystem(HFM,var_para);
HFM.Rec = getRecvector(HFM,var_para);
HFM.f = getSourcevector(HFM,var_para);
C_korn = computeKornsBetab(HFM);    

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
load('optimals0_Ricker.mat')

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

%% Training set for 30% change
minLam = abs(var_para.lami-0.30*(var_para.lami));%
maxLam = abs(var_para.lami+0.30*(var_para.lami));

minmu= abs(var_para.mui-0.30*(var_para.mui));
maxmu= abs(var_para.mui+0.30*(var_para.mui));

var_para.mui_min = minmu;
var_para.lam_min = minLam;
var_para.mui_max = maxmu;
var_para.lam_max = maxLam;


%% Training set
 var_para.Ntest = 1;
 rndvar = rand(var_para.Ntest,2);
% 
 lam_rand = (maxLam'-minLam').*rndvar(:,1) + minLam';
 mu_rand = (maxmu'-minmu').*rndvar(:,2) + minmu';
        
rand_lam_mu_id = 1;
%% Updating the stifness matrix for its use in the Newmark-beta method (Uncomment if different rand_lam_mu_id is selected)
K_splitHFM = cellfun(@(lam, mu, Kl, Ku) lam * Kl + mu * Ku, num2cell(lam_rand(rand_lam_mu_id,:)), num2cell(mu_rand(rand_lam_mu_id,:)), Ki_HFMlam.', Ki_HFMmu.', 'UniformOutput', false);
 HFM.K = K_splitHFM{1} + K_splitHFM{2} + K_splitHFM{3} + K_splitHFM{4} + K_splitHFM{5};
% Computing time domain seismograms using Newmark-beta method
[seismo_NBMx,seismo_NBMy] = solveTimeNewmark(HFM,var_para);

Nb = 300; % Maximum number of RB functions to be used 
   
%%%%%%%%%%%%%%%%
Aq = cell(11,1);
Aq{1} = M_HFM;
for qa = 2:6
    Aq{qa} = Ki_HFMlam{qa-1};    
    Aq{qa+5} = Ki_HFMmu{qa-1};      
end

%% POD-Greedy algorithm test
% loading the ROM
load('ROM_PODGreedy_m.mat','ROM_PODGreedy_m')
%%
 addpath('Reduced_basis_methods/')

Vk_it = ROM_PODGreedy_m.V(:,1:Nb);

% 
%  %%%%%%%%%%%%%% Project
 M_ROM = Vk_it'*(Aq{1}*Vk_it);
 bi_ROM =Vk_it'*(bi_HFM);
 K_ROMlam = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(2:6), 'UniformOutput', false);
 K_ROMmu = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(7:11), 'UniformOutput', false);

 %  %%%%%%%%%%%%%% Find RB seismograms

% Solving the ROM to find seismograms
[Seismo_PGTD] = getseismosTD_lam(M_ROM,bi_ROM,s_LF,RICKLF,var_para,Rec,Vk_it,lam_rand(rand_lam_mu_id,:),mu_rand(rand_lam_mu_id,:),K_ROMlam,K_ROMmu);
% Compute error Laplace domain

% Compute the L2 error time domain seismogram
relative_seismoTD = sqrt( (norm(var_para.dt*seismo_NBMx,2))^2 + (norm(var_para.dt*seismo_NBMy,2))^2);


er_NBM_PGxTD =  seismo_NBMx - Seismo_PGTD(1,:,1);
er_NBM_PGyTD =  seismo_NBMy - Seismo_PGTD(1,:,2);
Error_seismoTDL2_PODGreedy_NBM = (sqrt( (norm(var_para.dt*er_NBM_PGxTD,2))^2 + (norm(var_para.dt*er_NBM_PGyTD,2))^2))/...
                                relative_seismoTD;  
%
fprintf('Error POD-Greedy vsvs NBM: %2.2e\n',Error_seismoTDL2_PODGreedy_NBM)



% Greedy algorithm test

load('ROM_Greedy_m.mat','ROM_Greedy_m')

 %


Vk_it = ROM_Greedy_m.V(:,1:Nb);

       
% 
%  %%%%%%%%%%%%%% Project
 M_ROM = Vk_it'*(Aq{1}*Vk_it);
 bi_ROM =Vk_it'*(bi_HFM);
 K_ROMlam = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(2:6), 'UniformOutput', false);
 K_ROMmu = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(7:11), 'UniformOutput', false);

 %  %%%%%%%%%%%%%% Find RB seismograms

% Solving the ROM to find seismograms
[Seismo_GTD] = getseismosTD_lam(M_ROM,bi_ROM,s_LF,RICKLF,var_para,Rec,Vk_it,lam_rand(rand_lam_mu_id,:),mu_rand(rand_lam_mu_id,:),K_ROMlam,K_ROMmu);
% Compute error Laplace domain

% Compute the L2 error time domain seismogram
relative_seismoTD = sqrt( (norm(var_para.dt*seismo_NBMx,2))^2 + (norm(var_para.dt*seismo_NBMy,2))^2);
   

er_NBM_GredxTD =  seismo_NBMx - Seismo_GTD(1,:,1);
er_NBM_GredyTD =  seismo_NBMy - Seismo_GTD(1,:,2);
Error_seismoTDL2_Greedy_NBM = (sqrt( (norm(var_para.dt*er_NBM_GredxTD,2))^2 + (norm(var_para.dt*er_NBM_GredyTD,2))^2))/...
                                relative_seismoTD;  

fprintf('Error Greedy vs FOM vs NBM: %2.2e\n',Error_seismoTDL2_Greedy_NBM)


fprintf('DONE\n')


% Plotting the results

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

figure(2)
set(gcf,'position',[1000,550,550,500])
subplot(3,2,[1 2 3 4])
plot(var_para.t,(seismo_NBMx)/norm(seismo_NBMx,inf),'-v','Color','k','LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
hold on
plot(var_para.t,(Seismo_PGTD(1,:,1))/norm(seismo_NBMx,inf),':s','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
plot(var_para.t,(Seismo_GTD(1,:,1))/norm(seismo_NBMx,inf),'-.o','Color',cl_colorsdark{2},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))

hold off
set(gca, 'FontName', 'Arial','FontSize',16)
xlabel('Time (sec)','Interpreter','latex')
ylabel('Horizontal Component','Interpreter','latex')
legend('Newmark-beta method','POD-Greedy approx.','Greedy approx.','Interpreter','latex','numcolumns',2,'Location','southwest')

legend boxoff
ylim([-1 1])

subplot(3,2,[5 6])

plot(var_para.t,(er_NBM_PGxTD)/norm(seismo_NBMx,inf),':s','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
hold on
plot(var_para.t,(er_NBM_GredxTD)/norm(seismo_NBMx,inf),'-.o','Color',cl_colorsdark{2},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
hold off
ylim([-2e-3 2e-3])
set(gca, 'FontName', 'Arial','FontSize',16)
xlabel('Time (sec)','Interpreter','latex')
ylabel('Error','Interpreter','latex')
saveas(gcf,'SeismogramHcompTD_PODGreedy_and_Greedy_Fig9bc.png')

%%
figure(3)
set(gcf,'position',[1000,550,550,500])
subplot(3,2,[1 2 3 4])
plot(var_para.t,(seismo_NBMy)/norm(seismo_NBMy,inf),'-v','Color','k','LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
hold on
plot(var_para.t,(Seismo_PGTD(1,:,2))/norm(seismo_NBMy,inf),':s','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
plot(var_para.t,(Seismo_GTD(1,:,2))/norm(seismo_NBMy,inf),'-.o','Color',cl_colorsdark{2},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))

hold off
set(gca, 'FontName', 'Arial','FontSize',16)
xlabel('Time (sec)','Interpreter','latex')
ylabel('Vertical Component','Interpreter','latex')
legend('Newmark-beta method','POD-Greedy approx.','Greedy approx.','Interpreter','latex','numcolumns',2,'Location','southwest')
ylim([-1 1])

legend boxoff
subplot(3,2,[5 6])

plot(var_para.t,(er_NBM_PGyTD)/norm(seismo_NBMy,inf),':s','Color',cl_colorsdark{1},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
hold on
plot(var_para.t,(er_NBM_GredyTD)/norm(seismo_NBMy,inf),'-.o','Color',cl_colorsdark{2},'LineWidth',3,'Markersize',12,'Markerindices',round(linspace(1,numel(var_para.t),10)))
hold off
ylim([-2e-3 2e-3])
set(gca, 'FontName', 'Arial','FontSize',16)
xlabel('Time (sec)','Interpreter','latex')
ylabel('Error','Interpreter','latex')

rmpath('FOM_DATA')
saveas(gcf,'SeismogramVcompTD_PODGreedy_and_Greedy_Fig9bc.png')
close all

end