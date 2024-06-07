%%% This function file is used in the POD-Greedy and the Greedy algorithm 
% to compute seismograms and solution for changes in the Lame' parameters. 

function [Seismo_RBMTD,Seismo_RBMLD,uh] = getseismosTD_lam(M_HFM,bi_HFM,s_LF,RICKLF,var_para,Rec,V,lam_rand,mu_rand,Ki_lam,Ki_mu)

K_splitHFM = cellfun(@(lam, mu, Kl, Ku) lam * Kl + mu * Ku, num2cell(lam_rand), num2cell(mu_rand), Ki_lam.', Ki_mu.', 'UniformOutput', false);
K_HFM = K_splitHFM{1} + K_splitHFM{2} + K_splitHFM{3} + K_splitHFM{4} + K_splitHFM{5};
%toc
if ~(isempty(V))
    DOF = size(V,2);
else
    DOF = size(K_HFM,2);
end
course_idS = var_para.s_LFids;
Ns_LF = length(s_LF);
uh = zeros(DOF,Ns_LF);
for sid = 1:Ns_LF
    uh(:,sid) =(s_LF(sid)^2*M_HFM + K_HFM)\(RICKLF(sid)*bi_HFM);
end
if ~(isempty(V))
Snapshot_WeekRBPst_LF = V*uh;
else
Snapshot_WeekRBPst_LF = uh;
end
Seismo_RBMLD = zeros(1,Ns_LF,2);
Seismo_RBMTD = zeros(1,length(var_para.t),2);

% Seismo reconstrucion
    Snapshot_WeekRBPst_LFx = Rec(1:2:end,1)'*Snapshot_WeekRBPst_LF(1:2:end,:);
    Snapshot_WeekRBPst_LFy = Rec(2:2:end,1)'*Snapshot_WeekRBPst_LF(2:2:end,:);

    Seismo_RBMLD(1,:,1) = Snapshot_WeekRBPst_LFx;
    Seismo_RBMLD(1,:,2) = Snapshot_WeekRBPst_LFy;


    Snapshot_WeekRBPstx = zeros(1,var_para.Nz);
    Snapshot_WeekRBPstx(:,course_idS) = Snapshot_WeekRBPst_LFx;

    Snapshot_WeeRBkNegx = fliplr(real(Snapshot_WeekRBPstx(:,:)) - 1i*imag(Snapshot_WeekRBPstx(:,:)));
    seismoLD_RBMx = [Snapshot_WeeRBkNegx,Snapshot_WeekRBPstx];%[Snapshot_WeekPstF,Snapshot_WeekNeg];%

    Snapshot_WeekRBPsty = zeros(1,var_para.Nz);
    Snapshot_WeekRBPsty(:,course_idS) = Snapshot_WeekRBPst_LFy;

    Snapshot_WeeRBkNegy = fliplr(real(Snapshot_WeekRBPsty(:,:)) - 1i*imag(Snapshot_WeekRBPsty(:,:)));   
    seismoLD_RBMy = [Snapshot_WeeRBkNegy,Snapshot_WeekRBPsty];

    Seismo_RBMTD(1,:,1) = (Getweeksmethod(seismoLD_RBMx, var_para, var_para.s0(1), var_para.s0(2)))';                 
    Seismo_RBMTD(1,:,2)= (Getweeksmethod(seismoLD_RBMy, var_para, var_para.s0(1), var_para.s0(2)))';


end


function [ut] = Getweeksmethod(F, var_para, wr, wi)
    
    t   = var_para.t(:);                          % Make sure t is a column vector.    
    Nz = var_para.Nz;
    ap = GetLagCoef(F, Nz, wi);
    ut = DoClenshaw(Nz,t,wr,wi,real(ap));

end
%% Useful function files

function Lag_Coef = GetLagCoef(seismo_LD, Nz, wi)

FFTSamples = zeros(1,2*Nz,'double');            % Twice the samples as the number of coefficients
jdx = -Nz:(Nz-1);
theta_jhalf = (jdx+1/2)*pi/Nz;
Wtemp = exp(1i*theta_jhalf);

FFTSamples(jdx+Nz+1) = (2*wi./(1-Wtemp)).*seismo_LD;
TempCoef = fftshift(fft(fftshift(FFTSamples)))/(2*Nz); 

%Use only part of the TempCoef
Lag_Coef = TempCoef(Nz+1:2*Nz).*exp(-1i*(0:Nz-1)*pi/(2*Nz));

end 



function ftimevec = DoClenshaw(Nz,t,wr,wi,Lcf)

Ntimes = length(t);

Cpresent = zeros(Ntimes,1);
Cptwo = zeros(Ntimes,1);
Cpone = Lcf(Nz)*ones(Ntimes,1);

for kidx=(Nz-1):-1:1
 Cpresent = ((2*kidx-1-(2*wi*t))/kidx).*Cpone - (kidx/(kidx+1)).*Cptwo + Lcf(kidx)*ones(Ntimes,1);
 Cptwo = Cpone; 
 Cpone = Cpresent;  
end

ftimevec = exp((wr-wi)*t).*Cpresent;

end

