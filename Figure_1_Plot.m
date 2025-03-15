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
rmpath('FOM_DATA')

Nt = length(var_para.t);

M_HFM = HFM.M; K_HFM = HFM.K; bi_HFM = HFM.f;Xh = HFM.Xh; DOF = HFM.DOF; 
Ki_HFMlam = HFM.Ki_lam; Ki_HFMmu = HFM.Ki_mu; Rec = HFM.Rec;
%% Weeks method time approximation
% Step 1: find z(sWeek_j)
% Step 2: Compute a_p = wi/Nz *\sum_{j=-Nz}^{N_z - 1} C_j^p*z(sWeek_j) where
% exp(-1i*p*theta_jhalf)/(1 - exp(1i*theta_jhalf) )
% Step 3: Compute z(t) = sum_{p=0}^{N_z-1} a_p exp( (s0(1) - s0(2))*t)*L_p(2*s0(2)*t)

% Loading optimal parameters for Weeks method
addpath('DATA_for_s/')
load('optimals0_Ricker.mat')
var_para.s0 =s_optROM_Ricker; % s_optROM_Ricker; %[wr , wi] 
var_para.Nz = 608; % Number of frequency samples and Laguerre polynomials
wr = var_para.s0(1);
wi = var_para.s0(2);
jdx = -var_para.Nz:(var_para.Nz-1);
theta_jhalf = (jdx+1/2)*pi/var_para.Nz;
Cj = exp(1i*theta_jhalf);

sWeek = wr-wi*(Cj+1)./(Cj-1); % Same as wr + 1i*wi cot(theta_jhalf/2)
% Due to complex conjugation symmetry property we only need either + or - imaginary parts
sPstF = sWeek(var_para.Nz+1:end); % We take positive values here

%%
widths = [1*pi 1.5*pi 2*pi];
uh_widths = zeros(DOF,var_para.Nz,numel(widths));
energy= zeros(numel(widths),var_para.Nz);
var_para.smax = 25;
for w_id = 1:numel(widths)
    % Choice of Ricker wavelet parameters
    % Ricker wavelet parameters
    var_para.width = widths(w_id); % Central width
    var_para.t0 = 4*pi/var_para.width; % Time delay
    % Get maximum cut-off frequency for tolerance
    tolerance = 1e-4;
    %[var_para.smax] = getsmax(var_para.width,tolerance,fdualnorm,ldualnorm,K_HFM,M_HFM,var_para.Amp,C_korn);
    %fprintf('Cut-off omega_max:%2.2f for tolerance:%2.1e\n',var_para.smax,tolerance)    
    % We truncate sPstF to only our frequency values of interest
    [var_para.s_LFids] = find(abs(imag(sPstF(:)))<=var_para.smax);
    s_LF = sPstF(var_para.s_LFids);
    Ns_LF = length(s_LF);

    % Step 1: Finding z(sWeekj)
    %Ricker Laplace domain for s_LF (low-frequencies)  
    RICKLF = eval(var_para.source(var_para.Amp,var_para.width,var_para.t0,s_LF));
    tic
    uh_train = zeros(DOF,Ns_LF);

    % Solve full-order model for low-frequencies
    parfor sid = 1:Ns_LF
        uh_train(:,sid) = (s_LF(sid)^2*M_HFM + K_HFM)\(bi_HFM*RICKLF(sid));
    end
    uh_widths(:,1:Ns_LF,w_id) = uh_train;
    t_FOM = toc;
    for sid = 1:Ns_LF
        energy(w_id,sid) = sqrt(abs(uh_train(:,sid)'*(Xh*uh_train(:,sid)))); 
    end
end

%%
seismo_NBMx = zeros(numel(widths),numel(var_para.t));
seismo_NBMy = zeros(numel(widths),numel(var_para.t));
for w_id = 1:numel(widths)
    var_para.width = widths(w_id);
    var_para.t0 = 4*pi/var_para.width;
    [seismo_NBMx(w_id,:),seismo_NBMy(w_id,:)] = solveTimeNewmark(HFM,var_para);
end


%% maximum frequencies considered
smax_v = linspace(1,var_para.smax,25);
error_L2 = zeros(numel(widths),numel(smax_v));
for wid= 1:numel(widths)
    
    uh_train = uh_widths(:,1:Ns_LF,wid);
    for smax_id = 1:numel(smax_v)

        Seismo_TD_HFM = zeros(2,numel(var_para.t)); 
        [Seismo_TD_HFM(1,:),Seismo_TD_HFM(2,:)]= getTDseismo(var_para,uh_train,Rec,smax_v(smax_id),sPstF);
        relative_seismoTD = sqrt( (norm(var_para.dt*seismo_NBMx(w_id,:),2))^2 + (norm(var_para.dt*seismo_NBMy(w_id,:),2))^2);
        
        er_FOM_ROMxTD =  Seismo_TD_HFM(1,:) - seismo_NBMx(wid,:);%seismo_NBMx(w_id,:) - Seismo_TD_HFM(1,:);
        er_FOM_ROMyTD =  Seismo_TD_HFM(2,:) - seismo_NBMy(wid,:);
        error_L2(wid,smax_id) = (sqrt( (norm(var_para.dt*er_FOM_ROMxTD,2))^2 + (norm(var_para.dt*er_FOM_ROMyTD,2))^2))/...
                                        relative_seismoTD;  
    end
end

%% Plot
%%
cl_colorsdark = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4667    0.6745    0.1882]};

blue_lt  = [0    0.1961    0.3882];
red_lt = [0.6314    0.2471    0.0824];
green_lt = [0.3333    0.4902    0.1333];
cl_colorslight = {blue_lt,red_lt,green_lt};
figure(1)



for wid = 1:numel(widths)
    var_para.width = 1.0*pi; 
    var_para.t0 = 4*pi/var_para.width;
    % Time domain source
    subplot(1,3,1)
    ft = var_para.ft(var_para,var_para.t,var_para.t0);
    plot(var_para.t,ft,'-d','Color',cl_colorsdark{wid},'LineWidth',3)
    hold on
    subplot(1,3,2)

    source1 = eval(var_para.source(var_para.Amp,var_para.width,var_para.t0,s_LF));
    u1_normX = energy(1,1:Ns_LF);
    
    semilogy(imag(s_LF),abs(source1)/norm(source1,inf),':d','Color',cl_colorsdark{wid},'LineWidth',3, ...
           'Markersize',15,'MarkerIndices',50:25:numel(s_LF)-150)
    hold on
    
    semilogy(imag(s_LF),u1_normX/norm(u1_normX,inf),'-d','Color',cl_colorslight{wid}, ...
        'LineWidth',3,'Markersize',15,'MarkerIndices',50:35:numel(s_LF)-70)   
    set(gca, 'FontName', 'Arial','FontSize',25)
    ylim([1e-4 2e0])
    xlim([0 22])

end
hold off
xlabel('$s_I$','Interpreter','latex')
ylabel('${[|q(s;\alpha)| ; \|u(s;\alpha)\|_X}]_{0}$','Interpreter','latex')
legend('${|q(s;1.0\pi)|}$',...
       '${\|u(s;1.0\pi)\|_X} $',...
       '${|q(s;1.5\pi)|}$',...
       '${\|u(s;1.5\pi)\|_X} $',...
       '${|q(s;2.0\pi)|}$',...
       '${\|u(s;2.0\pi)\|_X} $',...
       'interpreter','latex','Location','northeast')
    legend boxoff

subplot(1,3,1)
xlabel('Time (sec)')
ylabel('Ricker wavelet')
hold off
set(gcf,'units','normalized','outerposition',[0 0 1 1])

subplot(1,3,3)
p_id = 1:2:25;
for wid=1:numel(widths)
    semilogy(smax_v(2:2:12),error_L2(wid,2:2:12),'--d','Linewidth',3,'Color',cl_colorsdark{wid},'Markersize',15,'Markerindices',1:1:50)
    xlabel('$s_{max}$','Interpreter','latex')
    ylabel('$[\|\hat{e}(s_{max})\|_{L^2} ]_{rel}$','Interpreter','latex')
    legend boxoff
    set(gca, 'FontName', 'Arial','FontSize',25)
    ylim([6e-4 3e0])
    xlim([0 22])
end
    legend('$\alpha=1.0\pi$','$\alpha=1.5\pi$','$\alpha=2.0\pi$','Interpreter','latex')
saveas(gcf,'Result_Fig1.png')



rmpath('DATA_for_s/')

%%
function [seismo_TD_Weekx,seismo_TD_Weeky]= getTDseismo(var_para,uh_train,Rec,s_max,sPstF)
% Calculate output on low-frequencies
Seismo_LD_sLFx = Rec(1:2:end,1)'*uh_train(1:2:end,:);
Seismo_LD_sLFy = Rec(2:2:end,1)'*uh_train(2:2:end,:);
% Post-processing to fill entries with smax>tolerance with value 0
Seismo_LD_sPstx = zeros(1,var_para.Nz);
Seismo_LD_sPsty = zeros(1,var_para.Nz);
%
Seismo_LD_sPstx(1,var_para.s_LFids) = Seismo_LD_sLFx;
Seismo_LD_sPsty(1,var_para.s_LFids) = Seismo_LD_sLFy;
Seismo_LD_sPstx = Seismo_LD_sPstx.*(imag(sPstF)<=s_max);
Seismo_LD_sPsty = Seismo_LD_sPsty.*(imag(sPstF)<=s_max);
% Complex conjugation to obtain the +/- counterparts of z(s)
Seismo_LD_sNegx = fliplr(real(Seismo_LD_sPstx) - 1i*imag(Seismo_LD_sPstx));
seismo_LD_sWeekx = [Seismo_LD_sNegx,Seismo_LD_sPstx];% Horizontal component
Seismo_LD_sNegy = fliplr(real(Seismo_LD_sPsty) - 1i*imag(Seismo_LD_sPsty));
seismo_LD_sWeeky = [Seismo_LD_sNegy,Seismo_LD_sPsty];% Veritcal component

wr = var_para.s0(1);
wi = var_para.s0(2);
jdx = -var_para.Nz:(var_para.Nz-1);
theta_jhalf = (jdx+1/2)*pi/var_para.Nz;
Cj = exp(1i*theta_jhalf);
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


function [seismo_x,seismo_y] = solveTimeNewmark(HFM,var_para)


ft = var_para.ft(var_para,var_para.t,var_para.t0);
n = HFM.DOF; 
dt = var_para.dt; % time-step
K = HFM.K; % Stiffness matrix
M = HFM.M; % Mass matrix
b = HFM.f; % Source location vector
u_old = zeros(n,1);
v_old = zeros(n,1);
Ntime = size(var_para.t,2);
gamma = 1/2; % Newmark-Beta constants
beta = 1/4; % Newmark-Beta constants (1/4 for implicit)

seismo_x = zeros(size(HFM.Rec,2),length(var_para.t));
seismo_y = zeros(size(HFM.Rec,2),length(var_para.t));
MK = M + beta*dt^2*K;        
[L,U,P,Q] = lu(MK); % Sparse LUPQ decomposition
a_old0 = L\(P*(b*ft(1)));
a_old1 = U\a_old0;
a_old = Q*a_old1;
for i=1:Ntime-1
   F = sparse(b*ft(i+1)); 
   utilda = u_old + v_old*dt + a_old *(0.5 - beta)*dt^2;
   vtilda = v_old + a_old*(1-gamma)*dt;
   a_new0 = (L)\(P*(F - K*utilda));
   a_new1 = U\a_new0;
   a_old = Q*a_new1;
   u_old  = utilda + a_old*beta*dt^2;
   v_old = vtilda + a_old*gamma*dt;
   seismo_x(:,i+1) = HFM.Rec(1:2:end,:)'*u_old(1:2:end,1);
   seismo_y(:,i+1) = HFM.Rec(2:2:end,:)'*u_old(2:2:end,1);
   if mod(i, 100) == 0
        fprintf('Process: %2.4f percent \n',i/(Ntime-1)*100)
   end
end
end
