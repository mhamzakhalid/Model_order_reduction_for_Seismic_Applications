function Figure_4_Plot(width_all,Ntrain)
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

M_HFM = HFM.M; K_HFM = HFM.K; bi_HFM = HFM.f;Xh = HFM.Xh; DOF = HFM.DOF;  Rec = HFM.Rec;
% Factorization of Xh = M + K 
[L_X,U_X,P_X,Q_X] =lu(Xh); 
% dualnorm of the output function (Receiver location)                        
[Output_Riesz] = Q_X*(U_X\(L_X\(P_X*Rec)));
ldualnorm = sqrt(abs(Rec'*Output_Riesz));
% dual norm of the Gaussian right hand side
source_Riesz=Q_X*(U_X\(L_X\(P_X*bi_HFM)));
fdualnorm = sqrt(abs(bi_HFM'*source_Riesz));


%%
Aq = cell(2,1);
Aq{1} = M_HFM; Aq{2} = K_HFM;
ROM_POD = cell(numel(width_all),1);
ROM_SPOD = cell(numel(width_all),1);
ROM_Greedy = cell(numel(width_all),1);
    
fprintf('Training set size: %d\n',Ntrain)
addpath('Reduced_basis_methods/')
for wid = 1:numel(width_all)
    fprintf('Running for alpha=%2.2fpi\n',width_all(wid)/pi)

    tolerance = 1e-4;
    [smax] = getsmax(width_all(wid),tolerance,fdualnorm,ldualnorm,K_HFM,M_HFM,var_para.Amp,C_korn);
    fprintf('smax=%2.2f \n',smax)

    %s_train = Training_set(:,wid);
    s_train_real = (0.62-0.12).*rand(Ntrain,1) + 0.12;
    s_train_imag = (smax - 0.01).*rand(Ntrain,1) + 0.01;
    s_train = s_train_real + 1i*s_train_imag;
    % Source term Ricker in frequency domain
    t0 = 4*pi/width_all(wid);
    RICK = eval(var_para.source(var_para.Amp,width_all(wid),t0,s_train));
    % Computing solutions using the full order model
    Ntrain = numel(s_train);
    uh_train = zeros(DOF,Ntrain);
    parfor sid = 1:Ntrain
        uh_train(:,sid) = (s_train(sid)^2*M_HFM + K_HFM)\(RICK(sid)*bi_HFM);
        fprintf('Process: %2.2f\n',sid/Ntrain*100)
    end
    fprintf('Solution computation done\n')
    ROM_POD{wid} = POD(HFM.Xh,uh_train,tolerance,250);
    fprintf('POD Done done\n')
    
    ROM_SPOD{wid} = POD(HFM.Xh,[real(uh_train),imag(uh_train)],tolerance,250);
    fprintf('SPOD Done done\n')
    
    seismo_HFMLD = zeros(2,Ntrain);
    seismo_HFMLD(1,:) = Rec(1:2:end,:)'*uh_train(1:2:end,:);
    seismo_HFMLD(2,:) = Rec(2:2:end,:)'*uh_train(2:2:end,:);
    
    ROM_Greedy{wid} = WeakGreedy_vecnorm(Aq,Xh,bi_HFM,Rec,var_para,s_train,uh_train,seismo_HFMLD,RICK,L_X,U_X,P_X,Q_X,C_korn);
    fprintf('Greedy Done done\n')
    
end

rmpath('FOM_DATA')
rmpath('Reduced_basis_methods/')

save('ROM_POD_s.mat','ROM_POD','-v7.3')
save('ROM_SPOD_s.mat','ROM_SPOD','-v7.3')
save('ROM_Greedy_s.mat','ROM_Greedy','-v7.3')

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

for wid = 1:numel(width_all)
    basis_id_wid = round(linspace(2,max(find(ROM_Greedy{wid}.maxDelta>1e-1)),15));
    % Plot illustrating the decay in the Singular values for choices of width = [1pi 1.5pi 2pi]
    subplot(1,3,1)
    semilogy(mark_id,ROM_POD{wid}.Sigma(mark_id),'-','Marker',marker_style{wid},'Color',cl_colorsdark{wid},'Linewidth',3.5,'Markersize',15)
    hold on
    xlim([0 200])
    xlabel('Singular values index: $k$','interpreter','latex')
    ylabel('Singular values: $\gamma_k$','interpreter','latex')
    legend('$\alpha=1.0\pi$','$\alpha=1.5\pi$','$\alpha=2.0\pi$','interpreter','latex','Location','northeast')
    set(gca, 'FontName', 'Arial','FontSize',29)
    legend boxoff
        
    
    % Ploting the Estimator for the Greedy algorithm for choices of width = [1pi 1.5pi 2pi]
    subplot(1,3,2)
    semilogy(basis_id_wid,ROM_Greedy{wid}.maxDelta(basis_id_wid),'-','Marker',marker_style{wid},'Color',cl_colorsdark{wid},'Linewidth',3.5,'Markersize',15)
    hold on
    xlim([0 200])
    yticks([1e-3 1e-1  1e2 1e4])
    ylim([1e-4 1e4])
    
    xlabel('No. of RB functions: $N_k$','interpreter','latex')
    ylabel('$\max_{s \in \tilde{I}_s} [\Delta_{k,f}(s)]_{0}$','interpreter','latex')
    legend('$\alpha=1.0\pi$','$\alpha=1.5\pi$','$\alpha=2.0\pi$','tol=$10^{-1}$','interpreter','latex','Location','northeast')
    set(gca, 'FontName', 'Arial','FontSize',29)
    legend boxoff
    % Ploting the True error for the Greedy algorithm for choices of width = [1pi 1.5pi 2pi]
    subplot(1,3,3)   
    semilogy(basis_id_wid,ROM_Greedy{wid}.maxTrue(basis_id_wid),':','Marker',marker_style{wid},'Color',cl_colorsdark{wid},'Linewidth',3.5,'Markersize',15)
    hold on
    xlim([0 200])
    xlabel('No. of RB functions: $N_k$','interpreter','latex')
    ylabel('$\max_{s \in \tilde{I}_s} [\|z_h(s) - z_k(s)\|_2]_{0}$','interpreter','latex')
    legend('$\alpha=1.0\pi$','$\alpha=1.5\pi$','$\alpha=2.0\pi$','interpreter','latex','Location','northeast')
    set(gca, 'FontName', 'Arial','FontSize',29)
    legend boxoff
    ylim([1e-4 1e4])

end
set(gcf,'position',[200,550,1550,500])
saveas(gcf,'Reduction_POD_and_Greedy_Fig4.png')


close all

end
