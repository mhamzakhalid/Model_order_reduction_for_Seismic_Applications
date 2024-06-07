function Compute_Optimal_Parameters()
%% Input Parameters 
addpath('FOM_DATA')
addpath('Optimal_parameters')
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

%width_all = 1*pi;%[1.0 1.5 2.0]*pi;
width = 2*pi;
%%
Ntrain = 1024;
fprintf('Training set size: %d\n',Ntrain)
addpath('Reduced_basis_methods/')
    fprintf('Running for alpha=%2.2fpi\n',width/pi)

    tolerance = 1e-4;
    [smax] = getsmax(width,tolerance,fdualnorm,ldualnorm,K_HFM,M_HFM,var_para.Amp,C_korn);
    fprintf('smax=%2.2f \n',smax)

    %s_train = Training_set(:,wid);
    s_train_real = (0.62-0.12).*rand(Ntrain,1) + 0.12;
    s_train_imag = (smax - 0.01).*rand(Ntrain,1) + 0.01;
    s_train = s_train_real + 1i*s_train_imag;
    % Source term Ricker in frequency domain
    t0 = 4*pi/width;
    RICK = eval(var_para.source(var_para.Amp,width,t0,s_train));
    % Computing solutions using the full order model
    Ntrain = numel(s_train);
    uh_train = zeros(DOF,Ntrain);
    parfor sid = 1:Ntrain
        uh_train(:,sid) = (s_train(sid)^2*M_HFM + K_HFM)\(RICK(sid)*bi_HFM);
        fprintf('Process: %2.2f\n',sid/Ntrain*100)
    end
    fprintf('Solution computation done\n')
    ROM_POD = POD(HFM.Xh,uh_train,tolerance,250);
    fprintf('POD Done done\n')
    

    Vk = ROM_POD.V(:,1:176);
    var_para.Nz = 608;
    var_para.width = 2.0*pi;
    var_para.t0 = 4*pi/var_para.width;
    fprintf('Optimal Parameter computation\n')
    [wr,wi] =  getoptimalparameterWeeks(M_HFM,K_HFM,Vk,bi_HFM,smax,Rec,var_para);
    fprintf('Optimal Parameter computation done\n')
    
    s_optROM_Ricker = [wr wi]
    save('optimal_s0.mat','s_optROM_Ricker')
rmpath('Reduced_basis_methods/')

end

