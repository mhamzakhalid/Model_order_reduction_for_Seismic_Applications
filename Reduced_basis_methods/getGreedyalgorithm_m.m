function ROM_Greedy_m = getGreedyalgorithm_m(Aq,bi_HFM,Xh,Rec,s_LF,RICKLF,lam_rand,mu_rand,estimator,var_para,Seismo_HFMTD,Seismo_HFMLD)

    res_ff = estimator.res_ff;
    Cj = estimator.Cj;
    CW = estimator.CW;
    beta_lb = estimator.beta_lb;
    L_X = estimator.L_X;
    U_X = estimator.U_X;
    P_X = estimator.P_X;
    Q_X = estimator.Q_X;
    THETA_F_Vec = estimator.THETA_F_Vec;
    ldualnorm = estimator.ldualnorm;
    clear estimator 
    Nk = 0;
    Nkmax = var_para.Nkmax;
    epsilontol = 1e-1;
    delta_k = 1.;
    DOF = size(Aq{1},2);
    Vk = zeros(DOF,var_para.Nkmax);
    
    maxTruekm = zeros(1,Nkmax);
    maxTruekL2m= zeros(1,Nkmax);
    maxTrueks= zeros(1,Nkmax);
    maxEstkm= zeros(1,Nkmax);
    maxEstks= zeros(1,Nkmax);
    maxEstk_idm = 1;
    maxEstk_ids = 1;
    parameters = zeros(Nkmax,2);
    parametersEst = zeros(Nkmax,2);
        
    Ns_LF= numel(s_LF);
    Nt = numel(var_para.t);
    relative_seismo = zeros(var_para.Ntrain,Ns_LF,Nkmax);
    Error_seismo_L2 = zeros(var_para.Ntrain,Ns_LF,Nkmax);
    Error_seismoTD_UB= zeros(var_para.Ntrain,1,Nkmax);
     
    Deltakt = zeros(var_para.Ntrain,var_para.Nkmax);
    Deltakf = zeros(var_para.Ntrain,Ns_LF,var_para.Nkmax);
    Deltakf_rel = zeros(var_para.Ntrain,Ns_LF,var_para.Nkmax);
    Deltakt_rel= zeros(var_para.Ntrain,var_para.Nkmax);
    relative_TD_term = zeros(var_para.Ntrain,var_para.Nkmax);
    

    t_onlineE = zeros(Nkmax,1);
    t_offlineE= zeros(Nkmax,1);
    t_seismorb= zeros(Nkmax,1);
    
    
    % Running the Greedy algorithm
    %
    while Nk<Nkmax && delta_k>epsilontol
        %%%%%%%%%%%%%% Update RB 
        Nk = Nk + 1;
           
            K_splitHFM = cellfun(@(lam, mu, Kl, Ku) lam * Kl + mu * Ku, num2cell(lam_rand(maxEstk_idm,:)), num2cell(mu_rand(maxEstk_idm,:)), Aq(2:6).', Aq(7:11).', 'UniformOutput', false);
            K_HFM = K_splitHFM{1} + K_splitHFM{2} + K_splitHFM{3} + K_splitHFM{4} + K_splitHFM{5};
            uh = (s_LF(maxEstk_ids)^2*Aq{1} + K_HFM)\(RICKLF(maxEstk_ids)*bi_HFM);
            if Nk==1
                zeta = Gram_Schmidt([], uh, Xh);        
            else
                zeta = Gram_Schmidt(Vk(:,1:Nk-1), uh, Xh);
            end
            Vk(:,Nk) = zeta;
            Vk_it = Vk(:,1:Nk);
            

       %  clear K_HFM uh zeta POD_basis
       % 
       %  %%%%%%%%%%%%%% Project
         M_ROM = Vk_it'*(Aq{1}*Vk_it);
         bi_ROM =Vk_it'*(bi_HFM);
         K_ROMlam = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(2:6), 'UniformOutput', false);
         K_ROMmu = cellfun(@(Ki) Vk_it'*(Ki*Vk_it), Aq(7:11), 'UniformOutput', false);
       %  %[M_ROM,K_ROMlam,K_ROMmu] = projectRBGPU(Aq{1},Aq(2:6),Aq(7:11),Vk_it);
       %  %%%%%%%%%%%%%% Find RB seismograms
         tic 
         Seismo_RBMTD = zeros(var_para.Ntrain,Nt,2);
         Seismo_RBMLD = zeros(var_para.Ntrain,Ns_LF,2);
         uk = zeros(Nk,Ns_LF,var_para.Ntrain);
       % 
         parfor lam_id = 1:var_para.Ntrain
             [Seismo_RBMTD(lam_id,:,:),Seismo_RBMLD(lam_id,:,:),uk(:,:,lam_id)] = getseismosTD_lam(M_ROM,bi_ROM,s_LF,RICKLF,var_para,Rec,Vk_it,lam_rand(lam_id,:),mu_rand(lam_id,:),K_ROMlam,K_ROMmu);
         end
         t_seismorb(Nk) = toc;
       % 
       %  %%%%%%%%%%%%%% Finding update parameter
     
         tic
        % Offline a Posteriori error estimator
        [dqq, Eqq]  = offline_residual_comp(Aq,bi_HFM,L_X,U_X,P_X,Q_X,Vk_it);

        t_offlineE(Nk) = toc;

        tic
        deltak_LF = zeros(var_para.Ntrain,Ns_LF);
         parfor lam_id =1:var_para.Ntrain
            deltak_LF(lam_id,:) = evalEstimateovers(s_LF,lam_rand(lam_id,:),mu_rand(lam_id,:),uk(:,:,lam_id),dqq,Eqq,res_ff,THETA_F_Vec,beta_lb);
         end
        t_onlineE(Nk) = toc;
          deltakFull = zeros(var_para.Ntrain,var_para.Nz);

          deltakFull(:,var_para.s_LFids) = deltak_LF;
          deltakFull = Cj.*[fliplr(deltakFull),deltakFull];
          sumdeltaFull = CW*sum(deltakFull,2);
          Deltakt(:,Nk) =  ldualnorm*sumdeltaFull;
          Deltakf(:,:,Nk) = ldualnorm*deltak_LF;%deltak_LF;
          
          %%%%%%%%%% True errors
    
    
        relative_seismoTD = zeros(var_para.Ntrain,1);
        er_FOM_ROMxTD = zeros(var_para.Ntrain,Nt);
        er_FOM_ROMyTD = zeros(var_para.Ntrain,Nt);
        Error_seismoTDL2= zeros(var_para.Ntrain,1);
           %
          for lam_id = 1:var_para.Ntrain
              % Compute error Laplace domain
              relative_seismo(lam_id,:,Nk) = abs(Seismo_RBMLD(lam_id,:,1) + Seismo_RBMLD(lam_id,:,2));
              Error_seismo_L2(lam_id,:,Nk) =  (abs(Seismo_HFMLD(lam_id,:,1) - Seismo_RBMLD(lam_id,:,1) +...
                                               (Seismo_HFMLD(lam_id,:,2) - Seismo_RBMLD(lam_id,:,2))) )/max(relative_seismo(lam_id,:,Nk));
                
              Deltakf_rel(lam_id,:,Nk) = Deltakf(lam_id,:,Nk)/max(relative_seismo(lam_id,:,Nk));
              % Compute error time domain
              relative_seismoTD(lam_id,1) = sqrt( (norm(var_para.dt*Seismo_RBMTD(lam_id,:,1),2))^2 + (norm(var_para.dt*Seismo_RBMTD(lam_id,:,2),2))^2);
    
              er_FOM_ROMxTD(lam_id,:) =  Seismo_HFMTD(lam_id,:,1) - Seismo_RBMTD(lam_id,:,1);
              er_FOM_ROMyTD(lam_id,:) =  Seismo_HFMTD(lam_id,:,2) - Seismo_RBMTD(lam_id,:,2);
              Error_seismoTDL2(lam_id,1) = (sqrt( (norm(var_para.dt*er_FOM_ROMxTD(lam_id,:),2))^2 + (norm(var_para.dt*er_FOM_ROMyTD(lam_id,:),2))^2))/...
                                                relative_seismoTD(lam_id,1);          
              % Relative part
              seismo_P_rel = zeros(1,var_para.Nz); 
              seismo_P_rel(1,var_para.s_LFids) = relative_seismo(lam_id,:,Nk) ;
              seismo_full_rel = Cj.*[fliplr(seismo_P_rel),seismo_P_rel];
              relative_TD_term(lam_id,Nk) = CW*sum(seismo_full_rel);
              
              seismo_P = zeros(1,var_para.Nz); 
              seismo_P(1,var_para.s_LFids) = Error_seismo_L2(lam_id,:,Nk)*max(relative_seismo(lam_id,:,Nk));
              seismo_full = Cj.*[fliplr(seismo_P),seismo_P];
              Error_seismoTD_UB(lam_id,1,Nk) = CW*sum(seismo_full)/relative_seismoTD(lam_id,1);%/relative_TD_term;
              
             Deltakt_rel(lam_id,Nk) = Deltakt(lam_id,Nk)/relative_seismoTD(lam_id,1);
 
          end 
           %%% For true error m parameter
            [maxTruekm(Nk),maxTruek_idm] = max(Error_seismoTD_UB(:,1,Nk));
            [maxTruekL2m(Nk),maxTruek_idL2m] = max(Error_seismoTDL2);
            [maxEstkm(Nk),maxEstk_idm] = max(Deltakt_rel(:,Nk));
            delta_k = maxEstkm(Nk);
            %%% For true error s parameter
            [maxTrueks(Nk),maxTruek_ids] = max(Error_seismo_L2(maxTruek_idm,:,Nk));
            %Deltakf = zeros(var_para.Ntrain,Ns_LF);
            [maxEstks(Nk),maxEstk_ids] = max(Deltakf(maxEstk_idm,:,Nk));
    
    
            parameters(Nk,1:2) = [maxTruek_idm,maxTruek_ids];
            parametersEst(Nk,1:2) = [maxEstk_idm,maxEstk_ids];
            fprintf('NRB: %d Max True:%2.2e Index: %d \n',Nk,  maxTruekm(Nk),maxTruek_idm)
            fprintf('NRB: %d Est True:%2.2e Index: %d \n',Nk,  maxEstkm(Nk),maxEstk_idm)
            fprintf('NRB: %d L2 True:%2.2e Index: %d \n',Nk,  maxTruekL2m(Nk),maxTruek_idL2m)
        
            fprintf(' Time RB: %2.2f \n Time offline_Est %2.2f \n Time Online_Est: %2.2f \n', t_seismorb(Nk),t_offlineE(Nk),t_onlineE(Nk))
             
        
    end
    
    ROM_Greedy_m.maxTrueL2m = maxTruekL2m;
    ROM_Greedy_m.maxTruem = maxTruekm;
    ROM_Greedy_m.maxEstimatorm = maxEstkm;
    ROM_Greedy_m.V = Vk;

end