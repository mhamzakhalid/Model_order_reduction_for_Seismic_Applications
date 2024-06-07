function getinfsupcomparison()
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
Ki_HFMlam = HFM.Ki_lam; Ki_HFMmu = HFM.Ki_mu; Rec = HFM.Rec;
% Factorization of Xh = M + K 
[L_X,U_X,P_X,Q_X] =lu(Xh); 
% dualnorm of the output function (Receiver location)                        
[Output_Riesz] = Q_X*(U_X\(L_X\(P_X*Rec)));
ldualnorm = sqrt(abs(Rec'*Output_Riesz));
% dual norm of the Gaussian right hand side
source_Riesz=Q_X*(U_X\(L_X\(P_X*bi_HFM)));
fdualnorm = sqrt(abs(bi_HFM'*source_Riesz));

addpath('Reduced_basis_methods/')
s_Reals = round(linspace(0.1,2,5),2);
s_imag = round(linspace(0.1,20,5),2);
beta_h = zeros(numel(s_Reals),numel(s_imag));
beta_lb = beta_h;
beta_star = beta_h;
for sr_id = 1:numel(s_Reals)
         s_val = s_Reals(sr_id) + 1i*s_imag;
         beta_lb(sr_id,:) = getBetaLB(Ki_HFMlam,Ki_HFMmu,M_HFM,var_para.lami,var_para.mui,s_val,'3',Xh,var_para,C_korn,var_para.lami,var_para.lami,var_para.mui,var_para.mui);
                 
         beta_h(sr_id,:) = getBetaLB(Ki_HFMlam,Ki_HFMmu,M_HFM,var_para.lami,var_para.mui,s_val,'2',Xh,var_para,C_korn,var_para.lami,var_para.lami,var_para.mui,var_para.mui);

         beta_star(sr_id,:) = getBetaLB(Ki_HFMlam,Ki_HFMmu,M_HFM,var_para.lami,var_para.mui,s_val,'4',Xh,var_para,C_korn,var_para.lami,var_para.lami,var_para.mui,var_para.mui);
end
rmpath('Reduced_basis_methods/')

save('Inf_sup.mat','beta_h','beta_star','beta_lb')
%% Plotting of the figure

figure(1)
subplot(1,3,1)
heatmap(round(s_imag,2),round(s_Reals,2),round(beta_h,2))
xlabel('s_R')
ylabel('s_I')
title('beta_h')
set(gca, 'FontName', 'Arial','FontSize',29)
subplot(1,3,2)
heatmap(round(s_imag,2),round(s_Reals,2),round(beta_lb,2))
xlabel('s_R')
ylabel('s_I')
title('beta_lb')
set(gca, 'FontName', 'Arial','FontSize',29)
subplot(1,3,3)
heatmap(round(s_imag,2),round(s_Reals,2),round(beta_star,2))
xlabel('$s_R$')
ylabel('$s_I$')
set(gca, 'FontName', 'Arial','FontSize',29)
xlabel('s_R')
ylabel('s_I')
title('beta_s (Uniform)')
set(gcf,'position',[200,550,1550,500])
saveas(gcf,'Inf_sup_comparison.png')


close all




end
