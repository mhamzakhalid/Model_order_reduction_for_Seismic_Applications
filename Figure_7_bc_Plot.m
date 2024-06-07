clear
clc
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

Nt = length(var_para.t);
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
%% Weeks method time approximation
% Step 1: find z(sWeek_j)
% Step 2: Compute a_p = wi/Nz *\sum_{j=-Nz}^{N_z - 1} C_j^p*z(sWeek_j) where
% exp(-1i*p*theta_jhalf)/(1 - exp(1i*theta_jhalf) )
% Step 3: Compute z(t) = sum_{p=0}^{N_z-1} a_p exp( (s0(1) - s0(2))*t)*L_p(2*s0(2)*t)

% Loading optimal parameters for Weeks method
addpath('DATA_for_s/')
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

var_para.width = 1.0*pi;
var_para.t0 = 4*pi/var_para.width;
tolerance = 1e-4;
[var_para.smax] = getsmax(var_para.width,tolerance,fdualnorm,ldualnorm,K_HFM,M_HFM,var_para.Amp,C_korn);
[var_para.s_LFids] = find(abs(imag(sPstF(:)))<=var_para.smax);

s_LF = sPstF(var_para.s_LFids); % Low frequencies
Ns_LF = length(s_LF);
RICKLF = eval(var_para.source(var_para.Amp,var_para.width,var_para.t0,s_LF));

%% Training set
s_train = s_LF(1:2:end);
N_train = numel(s_train);
% Step 1: Finding z(sWeekj)
%Ricker Laplace domain for s_LF (low-frequencies)  
tic
uh_train = zeros(DOF,N_train);
RICK_strain = eval(var_para.source(var_para.Amp,var_para.width,var_para.t0,s_train));
% Solve full-order model for low-frequencies
parfor sid = 1:N_train
    uh_train(:,sid) = (s_train(sid)^2*M_HFM + K_HFM)\(bi_HFM*RICK_strain(sid));
end
addpath('Reduced_basis_methods/')
ROM_POD = POD(HFM.Xh,uh_train,tolerance,80);

basis_id_wid = round(linspace(2,80,10));
%% Computing errors 
[error_seismo_POD] = gettimedomainerrors(HFM,ROM_POD,var_para,var_para.width,basis_id_wid);
%%
rmpath('Reduced_basis_methods/')
figure(1)
subplot(1,2,1)
semilogy(ROM_POD.Sigma,'-k','Linewidth',3)
xlabel('Singular values index: $k$','interpreter','latex')
ylabel('Singular values: $\gamma_k$','interpreter','latex')
set(gca, 'FontName', 'Arial','FontSize',27)
xlim([0 80])
subplot(1,2,2)
errorbarlog(basis_id_wid,mean(error_seismo_POD(:,:,1),2 ),std(error_seismo_POD(:,:,1),0,2),'-k','Linewidth',3.5 )
set(gca, 'FontName', 'Arial','FontSize',27)
xlabel('No. of RB functions: $N_k$','Interpreter','latex')
    ylabel('$[\|\hat{e}_k\|_{L^2} ]_{rel}$','Interpreter','latex')
 %   legend boxoff
set(gcf,'position',[200,550,1550,500])
saveas(gcf,'POD_Coarse.png')
rmpath('FOM_DATA')

close all

%%
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
        Vk = ROM.V(:,1:basis_id(Nk));
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
