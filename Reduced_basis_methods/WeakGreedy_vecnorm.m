function ROM_Greedy = WeakGreedy_vecnorm(Aq,Xh,bi_HFM,Rec,var,s_train,uh,seismo_HFMLD,RICKLF,L_X,U_X,P_X,Q_X,C_korn)

Ns_train = length(s_train);

% Dual norm of the output       
[Output_Riesz] = Q_X*(U_X\(L_X\(P_X*Rec)));
dual_norm_output = sqrt(abs(Rec'*Output_Riesz));

source_Riesz=Q_X*(U_X\(L_X\(P_X*bi_HFM)));
Cqq = source_Riesz'*bi_HFM; 
THETA_F_Vec = RICKLF.';
res_ff = conj(THETA_F_Vec).*THETA_F_Vec*Cqq;

beta_lb= getBetaLB([],[],[],[],[],s_train,'4',[],var,[],var.lami,var.lami,var.mui,var.mui);

Nk = 0;% Start with index
Nkmax = 180;%maximum number of Basis
epsilontol = 1e-1;% Tolerance for the estimator
delta_k = 1.; % Intializing the while loop for 
Vk = [];
skid = 1;

maxDeltak =  zeros(1,Nkmax);
maxTruek =  zeros(1,Nkmax);
maxTruek_Omega=  zeros(1,Nkmax);
maxDeltak_Omega=  zeros(1,Nkmax);

while Nk<Nkmax && delta_k>epsilontol
    Nk = Nk + 1;   
    tic
    zeta = Gram_Schmidt(Vk, uh(:,skid), Xh);
    Vk = full([Vk zeta]);
     % Offline a posteriori error estimator
    [dqq, Eqq]  = offline_residual_comp(Aq,bi_HFM,L_X,U_X,P_X,Q_X,Vk);
    % Project
    M_ROM = Vk'*(Aq{1}*Vk);
    K_ROM = Vk'*(Aq{2}*Vk);
    bi_ROM = Vk'*(bi_HFM);
    uk = zeros(size(Vk,2),Ns_train);
    parfor sid = 1:Ns_train
        uk(:,sid) = (s_train(sid)^2*M_ROM + K_ROM)\(RICKLF(sid)*bi_ROM);
    end
    uk_h  = Vk*uk;
    Seismo_RBMLD = zeros(2,Ns_train);
    Seismo_RBMLD(1,:) = Rec(1:2:end,1)'*uk_h(1:2:end,:);
    Seismo_RBMLD(2,:) = Rec(2:2:end,1)'*uk_h(2:2:end,:);
    % Compute error
      % Seismogram
      normalizor = norm(vecnorm([Seismo_RBMLD(1,:) ; Seismo_RBMLD(2,:)],2,1),inf);                  
      ErrorX_rel =  vecnorm([seismo_HFMLD(1,:) - Seismo_RBMLD(1,:) ;...
                                  seismo_HFMLD(2,:) - Seismo_RBMLD(2,:)],2,1)/normalizor;
      % Estimator online seismogram       
      residual_dualnorm = evalEstimateovers_fixedm(s_train,uk,dqq,Eqq,res_ff,THETA_F_Vec,beta_lb);
      Deltak = dual_norm_output*residual_dualnorm/normalizor; 
      % Full domain
      normalizor_Omega = norm(sqrt(abs(sum((uh'*(Xh*uh)), 1))),inf); 
      er = uh - uk_h;
      ErrorX_Omega = sqrt(abs(sum((er'*(Xh*er)), 1)));
 
      Deltak_Omega = residual_dualnorm/normalizor_Omega; % Estimator full domain scaled
      ErrorX_Omega = ErrorX_Omega/normalizor_Omega; % Error full domain scaled
   
    %%
       % Finding the update parameter
       [maxDeltak(Nk),maxDeltak_id] = max(Deltak);
       %%% For estimator
       delta_k = maxDeltak(Nk); skid = maxDeltak_id;
       %%% For true error        
       [maxTruek(Nk)] = max(ErrorX_rel);
      
       [maxTruek_Omega(Nk)] = max(ErrorX_Omega);
              
       [maxDeltak_Omega(Nk)] = max(Deltak_Omega);

       fprintf('NRB: %d Max Estimate:%2.2e\n',Nk,  maxDeltak(Nk))
       fprintf('NRB: %d Max True Seismo:%2.2e\n',Nk,  maxTruek(Nk))
       fprintf('NRB: %d Max True Full domain:%2.2e\n',Nk,  maxTruek_Omega(Nk))

end
%
ROM_Greedy.maxDelta = maxDeltak;
ROM_Greedy.maxTrue = maxTruek;
ROM_Greedy.V = Vk;
ROM_Greedy.maxDelta_Omega = maxDeltak_Omega;
ROM_Greedy.maxTrue_Omega = maxTruek_Omega;

end


