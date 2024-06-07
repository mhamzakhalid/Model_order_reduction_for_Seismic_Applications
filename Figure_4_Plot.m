function Figure_4_Plot(width_all,Ntest)

%parpool(8)
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
Rec = HFM.Rec;
% Factorization of Xh = M + K 
[L_X,U_X,P_X,Q_X] =lu(Xh); 
% dualnorm of the output function (Receiver location)                        
[Output_Riesz] = Q_X*(U_X\(L_X\(P_X*Rec)));
ldualnorm = sqrt(abs(Rec'*Output_Riesz));
% dual norm of the Gaussian right hand side
source_Riesz=Q_X*(U_X\(L_X\(P_X*bi_HFM)));
fdualnorm = sqrt(abs(bi_HFM'*source_Riesz));

% width_all = [1*pi 1.5*pi];%[1.0 1.5 2.0]*pi;

%%
% Ntest = 512;

Nb = 15; %maximum iterations to reach the total number of RB functions
load('ROM_POD_s.mat','ROM_POD')
load('ROM_Greedy_s.mat','ROM_Greedy')
load('ROM_SPOD_s.mat','ROM_SPOD')

Test_set = cell(numel(width_all),1);
ErrorX_Omega_POD = zeros(Nb,Ntest,numel(width_all));
ErrorX_rel_POD = zeros(Nb,Ntest,numel(width_all));
ErrorX_Omega_Greedy = zeros(Nb,Ntest,numel(width_all));
ErrorX_rel_Greedy = zeros(Nb,Ntest,numel(width_all));
nb_val = [max(find(ROM_Greedy{1}.maxDelta>1e-1)) max(find(ROM_Greedy{2}.maxDelta>1e-1)) max(find(ROM_Greedy{3}.maxDelta>1e-1))];
basis_id_wid = zeros(3,15);
fprintf('Training set size: %d\n',Ntest)
for wid = 1:numel(width_all)
    fprintf('Running for alpha=%2.2fpi\n',width_all(wid))
    width = width_all(wid);
    t0 = 4*pi/width;
    tolerance = 1e-4;
    basis_id = round(linspace(2,nb_val(wid),Nb));
    basis_id_wid(wid,:) = basis_id;
    [smax] = getsmax(width_all(wid),tolerance,fdualnorm,ldualnorm,K_HFM,M_HFM,var_para.Amp,C_korn);
    fprintf('smax=%2.2f \n',smax)

    % Test test for the reduced basis
    s_test_real = (0.62-0.12).*rand(Ntest,1) + 0.12;
    s_test_imag = (smax - 0.01).*rand(Ntest,1) + 0.01;
    s_test = s_test_real + 1i*s_test_imag;
    Test_set{wid} = s_test;
    RICK = eval(var_para.source(var_para.Amp,width_all(wid),t0,s_test));
    % Computing solutions using the full order model for testing
    Ntest = numel(s_test);
    uh_test = zeros(DOF,Ntest);
    parfor sid = 1:Ntest
        uh_test(:,sid) = (s_test(sid)^2*M_HFM + K_HFM)\(RICK(sid)*bi_HFM);
        fprintf('Process FOM solution computation: %2.2f\n',sid/Ntest*100)
    end
    seismo_HFMLD = zeros(2,Ntest);
    seismo_HFMLD(1,:) = Rec(1:2:end,1)'*uh_test(1:2:end,:);
    seismo_HFMLD(2,:) = Rec(2:2:end,1)'*uh_test(2:2:end,:);
    %% Frequency domain error for the POD method
    for Nk = 1:numel(basis_id)
        Vk = ROM_POD{wid}.V(:,1:basis_id(Nk));
        M_ROM_POD = Vk'*(M_HFM*Vk);
        K_ROM_POD = Vk'*(K_HFM*Vk);
        bi_ROM_POD = Vk'*bi_HFM;
        uk_test_POD = zeros(basis_id(Nk),Ntest);
        parfor sid = 1:Ntest
            uk_test_POD(:,sid) = (s_test(sid)^2*M_ROM_POD + K_ROM_POD)\(RICK(sid)*bi_ROM_POD);
            fprintf('Process: %2.2f\n',sid/Ntest*100)
        end
        uk_h_POD  = Vk*uk_test_POD;
        Seismo_RBMLD_POD = zeros(2,Ntest);
        Seismo_RBMLD_POD(1,:) = Rec(1:2:end,1)'*uk_h_POD(1:2:end,:);
        Seismo_RBMLD_POD(2,:) = Rec(2:2:end,1)'*uk_h_POD(2:2:end,:);
        %%%% Compute error
           
        % Full domain
        
        normalizor_Omega = sqrt(abs(sum((uh_test'*(Xh*uh_test)), 1)));
        er_POD = uh_test - uk_h_POD;
        ErrorX_Omega_POD(Nk,:,wid) = sqrt(abs(sum((er_POD'*(Xh*er_POD)), 1)))/norm(normalizor_Omega,inf);
        
        % Seismograms
        normalizor = norm(vecnorm([Seismo_RBMLD_POD(1,:) ; Seismo_RBMLD_POD(2,:)],2,1),inf);
                  
        ErrorX_rel_POD(Nk,:,wid) =  vecnorm([seismo_HFMLD(1,:) - Seismo_RBMLD_POD(1,:) ;seismo_HFMLD(2,:) - Seismo_RBMLD_POD(2,:)],2,1)/normalizor;
  
    end
    %% Frequency domain error for the Greedy algorithm

    for Nk = 1:numel(basis_id)
        Vk = ROM_Greedy{wid}.V(:,1:basis_id(Nk));
        M_ROM_Greedy = Vk'*(M_HFM*Vk);
        K_ROM_Greedy = Vk'*(K_HFM*Vk);
        bi_ROM_Greedy = Vk'*bi_HFM;
        uk_test_Greedy = zeros(basis_id(Nk),Ntest);
        parfor sid = 1:Ntest
            uk_test_Greedy(:,sid) = (s_test(sid)^2*M_ROM_Greedy + K_ROM_Greedy)\(RICK(sid)*bi_ROM_Greedy);
            fprintf('Process: %2.2f\n',sid/Ntest*100)
        end
        uk_h_Greedy  = Vk*uk_test_Greedy;
        Seismo_RBMLD_Greedy = zeros(2,Ntest);
        Seismo_RBMLD_Greedy(1,:) = Rec(1:2:end,1)'*uk_h_Greedy(1:2:end,:);
        Seismo_RBMLD_Greedy(2,:) = Rec(2:2:end,1)'*uk_h_Greedy(2:2:end,:);
        %%%% Compute error           
        % Full domain        
        normalizor_Omega = sqrt(abs(sum((uh_test'*(Xh*uh_test)), 1)));
        er = uh_test - uk_h_Greedy;
        ErrorX_Omega_Greedy(Nk,:,wid) = sqrt(abs(sum((er'*(Xh*er)), 1)))/norm(normalizor_Omega,inf);     
        % Seismograms
        normalizor = norm(vecnorm([Seismo_RBMLD_Greedy(1,:) ; Seismo_RBMLD_Greedy(2,:)],2,1),inf);                  
        ErrorX_rel_Greedy(Nk,:,wid) =  vecnorm([seismo_HFMLD(1,:) - Seismo_RBMLD_Greedy(1,:) ;seismo_HFMLD(2,:) - Seismo_RBMLD_Greedy(2,:)],2,1)/normalizor;
    end
   
end
save('ROM_POD_Test_FD.mat','ErrorX_rel_POD','ErrorX_Omega_POD','Test_set','-v7.3')
save('ROM_Greedy_Test_FD.mat','ErrorX_rel_Greedy','ErrorX_Omega_Greedy','Test_set','-v7.3')



%% Computing time domain errors
fprintf('Running ROM Time domain Greedy\n')
[error_seismo_Greedy] = gettimedomainerrors(HFM,ROM_Greedy,var_para,width_all,basis_id_wid);
fprintf('Done ROM Time domain Greedy\n')
fprintf('Running ROM Time domain POD\n')
[error_seismo_POD] = gettimedomainerrors(HFM,ROM_POD,var_para,width_all,basis_id_wid);
fprintf('Done ROM Time domain Greedy\n')
fprintf('Running ROM Time domain S-POD\n')
[error_seismo_SPOD] = gettimedomainerrors(HFM,ROM_SPOD,var_para,width_all,basis_id_wid);
fprintf('Done ROM Time domain S-POD\n')
save('ROM_Greedy_Test_TD.mat','error_seismo_Greedy','-v7.3')
save('ROM_POD_Test_TD.mat','error_seismo_POD','-v7.3')
save('ROM_SPOD_Test_TD.mat','error_seismo_SPOD','-v7.3')

%% Plotting
blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
cl_colorsdark = {blue, red,green, brown, purple, black};

blue = [114 147 203]./290;
red = [211 94 96]./255;
black = [128 133 133]./255;
green = [132 186 91]./250;
brown = [171 104 87]./255;
purple = [144 103 167]./255;
cl_colorslight = {blue, red, green,purple, black, ...
              brown};
marker_style = {"<","o","s"};    
mark_id = round(linspace(2,200,15));

figure(1)
for wid = 1:numel(width_all)
    % Time domain errors
    subplot(1,3,3)
    errorbarlog(basis_id_wid(wid,:),mean(error_seismo_POD(:,:,wid),2 ),std(error_seismo_POD(:,:,wid),0,2),'-','Color',cl_colorsdark{wid},'Linewidth',3.5 )
    hold on
    errorbarlog(basis_id_wid(wid,:),mean(error_seismo_Greedy(:,:,wid),2 ),std(error_seismo_Greedy(:,:,wid),0,2),':','Color',cl_colorsdark{wid},'Linewidth',3.5 )

    % Frequency domain Errors
    
   
    % % For the Greedy algorithm
    ErrorX_rel_Greedy_meanbar = ErrorX_Omega_Greedy; %
    subplot(1,3,2)
    E1 = ErrorX_rel_Greedy_meanbar(:,:,wid);
    errorbarlog(basis_id_wid(wid,:),mean(E1,2 ),std(E1,1,2 ),'-.<','Color',cl_colorsdark{wid},'Linewidth',3 )
    hold on
    E1_s = ErrorX_rel_Greedy(:,:,wid);
    semilogy(basis_id_wid(wid,:),mean(E1_s,2),'-.<','Color',cl_colorslight{wid},'Linewidth',3 ,'Markersize',10)
    % % For the POD method
    ErrorX_rel_POD_meanbar = ErrorX_Omega_POD; 
    subplot(1,3,1)
    E2 = ErrorX_rel_POD_meanbar(:,:,wid);
    errorbarlog(basis_id_wid(wid,:),mean(E2,2 ),std(E2,1,2 ),'-.<','Color',cl_colorsdark{wid},'Linewidth',3 )
    hold on
    E2_s = ErrorX_rel_POD(:,:,wid);
    semilogy(basis_id_wid(wid,:),mean(E2_s,2),'-.<','Color',cl_colorslight{wid},'Linewidth',3 ,'Markersize',12)
end

    subplot(1,3,1)
    hold off
    ylim([1e-6 1e0])
    
    set(gca, 'FontName', 'Arial','FontSize',25)
    xlabel('No of RB functions: $N_k$','Interpreter','latex')
    ylabel('Mean$_{s\in \tilde{I}_t} [e_k(s)]_{0}$','Interpreter','latex')
    legend('POD (full domain) $\alpha=1.0\pi$','POD (Seismo) $\alpha=1.0\pi$',...
            'POD (full domain)$\alpha=1.5\pi$','POD (Seismo) $\alpha=1.5\pi$',...
            'POD (full domain) $\alpha=2.0\pi$','POD (Seismo) $\alpha=2.0\pi$',... 
            'interpreter','latex','Location','northeast','FontSize',20); 
    legend boxoff 

    subplot(1,3,2)
    hold off
    ylim([1e-6 1e0])
    
    set(gca, 'FontName', 'Arial','FontSize',25)
    xlabel('No of RB functions: $N_k$','Interpreter','latex')
    ylabel('Mean$_{s\in \tilde{I}_t} [e_k(s)]_{0}$','Interpreter','latex')
    legend('Greedy (full domain) $\alpha=1.0\pi$','Greedy (Seismo) $\alpha=1.0\pi$',...
            'Greedy (full domain)$\alpha=1.5\pi$','Greedy (Seismo) $\alpha=1.5\pi$',...
            'Greedy (full domain) $\alpha=2.0\pi$','Greedy (Seismo) $\alpha=2.0\pi$',... 
            'interpreter','latex','Location','northeast','FontSize',20); 
    legend boxoff 

    subplot(1,3,3)       
    xlabel('No. of RB functions: $N_k$','Interpreter','latex')
    ylabel('Mean$_i [\|e_{k,i}\|_{L^2}]_{rel}$','Interpreter','latex')
    set(gca, 'FontName', 'Arial','FontSize',29)
    legend('POD $\alpha=1.0\pi$','Greedy $\alpha=1.0\pi$',...
            'POD $\alpha=1.5\pi$','Greedy $\alpha=1.5\pi$',...
            'POD $\alpha=2.0\pi$','Greedy $\alpha=2.0\pi$',... 
            'interpreter','latex','Location','northeast','FontSize',20); 
    legend boxoff 
    hold off


set(gcf,'position',[200,550,1550,500])
saveas(gcf,'Reduction_POD_and_Greedy_Fig4.png')

close all

figure(2)
for wid=1:numel(width_all)
    subplot(1,2,2)
    errorbarlog(basis_id_wid(wid,:),mean(error_seismo_POD(:,:,wid),2 ),std(error_seismo_POD(:,:,wid),0,2),'-','Color',cl_colorsdark{wid},'Linewidth',3.5 )
    hold on
    errorbarlog(basis_id_wid(wid,:),mean(error_seismo_SPOD(:,:,wid),2 ),std(error_seismo_SPOD(:,:,wid),0,2),':','Color',cl_colorslight{wid},'Linewidth',3.5 )
    subplot(1,2,1)
    semilogy(ROM_POD{wid}.Sigma,'-','Color',cl_colorsdark{wid},'Linewidth',3.5 )
    hold on
    semilogy(ROM_SPOD{wid}.Sigma,':','Color',cl_colorslight{wid},'Linewidth',3.5 )
   
end
    subplot(1,2,2)
    hold off
    legend('POD $\alpha=1.0\pi$','SPOD $\alpha=1.0\pi$',...
            'POD $\alpha=1.5\pi$','SPOD $\alpha=1.5\pi$',...
            'POD $\alpha=2.0\pi$','SPOD $\alpha=2.0\pi$',... 
            'interpreter','latex','Location','northeast','FontSize',20); 
    legend boxoff 
    hold off
    xlabel('No. of RB functions: $N_k$','Interpreter','latex')
    ylabel('Mean$_i [\|e_{k,i}\|_{L^2}]_{rel}$','Interpreter','latex')
    set(gca, 'FontName', 'Arial','FontSize',29)

    subplot(1,2,1)
    legend('POD $\alpha=1.0\pi$','SPOD $\alpha=1.0\pi$',...
            'POD $\alpha=1.5\pi$','SPOD $\alpha=1.5\pi$',...
            'POD $\alpha=2.0\pi$','SPOD $\alpha=2.0\pi$',... 
            'interpreter','latex','Location','northeast','FontSize',20); 
    legend boxoff 
    hold off
    xlabel('No. of RB functions: $N_k$','Interpreter','latex')
    ylabel('Singular values $\gamma_k$','Interpreter','latex')
    set(gca, 'FontName', 'Arial','FontSize',29)
    hold off

    set(gcf,'position',[200,550,1550,500])
    saveas(gcf,'Reduction_POD_and_SPOD_Supplementary.png')

end
%% Function files
function [error_seismo,seismo_TD_Weekx,seismo_TD_Weeky] = gettimedomainerrors(HFM,ROM,var_para,width_all,basis_id_wid)
% We consider the following 5 receivers within our region of interest
Multiple_Rec_x = [0.736161035200000 ;...
                  0.842559309800000;...
                  0.947519769900000;...
                  1.052480230000000;...
                  1.157440690200000]*1.0e4;
Multiple_Rec_y = 2.0e4*ones(5,1);
HFM.Rec = zeros(HFM.DOF,5);
for rec_id = 1:5
   var_para.x0r = Multiple_Rec_x(rec_id); var_para.y0r = Multiple_Rec_y(rec_id);
   HFM.Rec(:,rec_id) = getRecvector(HFM,var_para);
end

C_korn = computeKornsBetab(HFM);    

M_HFM = HFM.M; K_HFM = HFM.K; bi_HFM = HFM.f;Xh = HFM.Xh; DOF = HFM.DOF; 
Rec = HFM.Rec;
% Factorization of Xh = M + K 
[L_X,U_X,P_X,Q_X] =lu(Xh); 
% dualnorm of the output function (Receiver location)                        
[Output_Riesz] = Q_X*(U_X\(L_X\(P_X*Rec(:,3))));
ldualnorm = sqrt(abs(Rec(:,3)'*Output_Riesz));
% dual norm of the Gaussian right hand side
source_Riesz=Q_X*(U_X\(L_X\(P_X*bi_HFM)));
fdualnorm = sqrt(abs(bi_HFM'*source_Riesz));

%%% Newmark-beta method
seismo_NBMx = zeros(size(Rec,2),numel(var_para.t),numel(width_all));
seismo_NBMy = zeros(size(Rec,2),numel(var_para.t),numel(width_all));
addpath('DATA_for_s') 
load('seismo_NBM_5_rec_3_width.mat') % These are pre-computed solutions. Otherwise uncomment the following
% t_NBM = zeros(3,1);
% for wid = 1:numel(width_all)
%     var_para.width = width_all(wid);
%     var_para.t0 = 4*pi/var_para.width;
%     tic
%     [seismo_NBMx(:,:,wid),seismo_NBMy(:,:,wid)] = solveTimeNewmark(HFM,var_para);
%     t_NBM(wid) = toc
% 
% end
% save('seismo_NBM_5_rec_3_width.mat','seismo_NBMx','seismo_NBMy')

%%% Weeks method time approximation
% Step 1: find z(sWeek_j)
% Step 2: Compute a_p = wi/Nz *\sum_{j=-Nz}^{N_z - 1} C_j^p*z(sWeek_j) where
% exp(-1i*p*theta_jhalf)/(1 - exp(1i*theta_jhalf) )
% Step 3: Compute z(t) = sum_{p=0}^{N_z-1} a_p exp( (s0(1) - s0(2))*t)*L_p(2*s0(2)*t)

% Loading optimal parameters for Weeks method
load('optimals0_Ricker.mat')
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

%%%
Nr = size(Rec,2);
error_seismo = zeros(size(basis_id_wid,2),Nr,numel(width_all));
for wid = 1:numel(width_all)
    width = width_all(wid);
    t0 = 4*pi/width;
    tolerance = 1e-4;
    [smax] = getsmax(width_all(wid),tolerance,fdualnorm,ldualnorm,K_HFM,M_HFM,var_para.Amp,C_korn);
    % We truncate sPstF to only our frequency values of interest
    [var_para.s_LFids] = find(abs(imag(sPstF(:)))<=smax);
    s_LF = sPstF(var_para.s_LFids);
    Ns_LF = length(s_LF);
    % Consider the number of basis functions over which we compute
    % seismograms
    basis_id = basis_id_wid(wid,:);                                                                                  
    RICK = eval(var_para.source(var_para.Amp,width_all(wid),t0,s_LF));
    % Frequency domain error
    for Nk = 1:size(basis_id,2)
        Vk = ROM{wid}.V(:,1:basis_id(Nk));
        M_ROM = Vk'*(M_HFM*Vk);
        K_ROM = Vk'*(K_HFM*Vk);
        bi_ROM = Vk'*bi_HFM;
        uk_test = zeros(basis_id(Nk),Ns_LF);
        parfor sid = 1:Ns_LF
            uk_test(:,sid) = (s_LF(sid)^2*M_ROM + K_ROM)\(RICK(sid)*bi_ROM);
        end
        uk_h  = Vk*uk_test;
        
        for rec_id = 1:size(Rec,2)
            Seismo_LD_sLFx = Rec(1:2:end,rec_id)'*uk_h(1:2:end,:);
            Seismo_LD_sLFy = Rec(2:2:end,rec_id)'*uk_h(2:2:end,:);
            % Compute error
            Seismo_LD_sPstx = zeros(1,var_para.Nz);
            Seismo_LD_sPsty = zeros(1,var_para.Nz);
            %
            Seismo_LD_sPstx(1,var_para.s_LFids) = Seismo_LD_sLFx;
            Seismo_LD_sPsty(1,var_para.s_LFids) = Seismo_LD_sLFy;
            % Complex conjugation to obtain the +/- counterparts of z(s)
            Seismo_LD_sNegx = fliplr(real(Seismo_LD_sPstx) - 1i*imag(Seismo_LD_sPstx));
            seismo_LD_sWeekx = [Seismo_LD_sNegx,Seismo_LD_sPstx];% Horizontal component
            Seismo_LD_sNegy = fliplr(real(Seismo_LD_sPsty) - 1i*imag(Seismo_LD_sPsty));
            seismo_LD_sWeeky = [Seismo_LD_sNegy,Seismo_LD_sPsty];% Veritcal component
            %
            % Step 2: Find expansion coefficients a_p = wi/Nz* \sum_{j=-N_z}^{N_z-1} C_j^p z(sWeek_j)
            % We do it using FFT but can be done using a midpoint evalulation as well.
            p = 0:var_para.Nz-1; 
            FFTSamplesx = zeros(1,2*var_para.Nz);            % Twice the samples as the number of coefficients
            FFTSamplesx(1,jdx+var_para.Nz+1) = (2*wi./(1-Cj)).*seismo_LD_sWeekx;
            %Note the order: FFTSamples(1,1:2*N)
            TempCoefx = fftshift(fft(fftshift(FFTSamplesx)))/(2*var_para.Nz); 
            %Use only part of the TempCoef
            Lag_Coefx = real(TempCoefx(var_para.Nz+1:2*var_para.Nz).*exp(-1i*p*pi/(2*var_para.Nz)))';
            
            FFTSamplesy = zeros(1,2*var_para.Nz);            % Twice the samples as the number of coefficients
            FFTSamplesy(1,jdx+var_para.Nz+1) = (2*wi./(1-Cj)).*seismo_LD_sWeeky;
            %Note the order: FFTSamples(1,1:2*N)
            TempCoefy = fftshift(fft(fftshift(FFTSamplesy)))/(2*var_para.Nz); 
            %Use only part of the TempCoef
            Lag_Coefy = real(TempCoefy(var_para.Nz+1:2*var_para.Nz).*exp(-1i*p*pi/(2*var_para.Nz)))';
            
            
            % Step 3: Clenshaw algorithm to compute z(t)
            seismo_TD_Weekx = doClenshaw(var_para,Lag_Coefx,wr,wi);
            seismo_TD_Weeky = doClenshaw(var_para,Lag_Coefy,wr,wi);
            
            % Computing errors
            relative_seismoTD = sqrt( (norm(var_para.dt*seismo_NBMx(rec_id,:,wid),2))^2 + (norm(var_para.dt*seismo_NBMy(rec_id,:,wid),2))^2);

            er_FOM_ROMxTD =  seismo_TD_Weekx - seismo_NBMx(rec_id,:,wid);
            er_FOM_ROMyTD =  seismo_TD_Weeky - seismo_NBMy(rec_id,:,wid);
      
            error_seismo(Nk,rec_id,wid) = (sqrt( (norm(var_para.dt*er_FOM_ROMxTD,2))^2 + (norm(var_para.dt*er_FOM_ROMyTD,2))^2))/...
                                            relative_seismoTD;
        end
            fprintf('Process: %2.2f\n',Nk/size(basis_id,2)*100)

    end
end
rmpath('DATA_for_s')
end

function seismo_TD_Week = doClenshaw(var_para,Lag_Coef,wr,wi)

t = var_para.t(:);
Nt = numel(t);
Cpresent = zeros(Nt,1);
Cptwo = zeros(Nt,1);
Cpone = Lag_Coef(var_para.Nz)*ones(Nt,1);
for kidx=(var_para.Nz-1):-1:1
 Cpresent = ((2*kidx-1-(2*wi*t))/kidx).*Cpone - (kidx/(kidx+1)).*Cptwo + Lag_Coef(kidx)*ones(Nt,1);
 Cptwo = Cpone; 
 Cpone = Cpresent;  
end

seismo_TD_Week = (exp((wr-wi)*t).*Cpresent)';

end
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