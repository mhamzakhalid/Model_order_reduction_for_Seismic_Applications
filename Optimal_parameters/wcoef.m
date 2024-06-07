

function ap = wcoef(F,Nz,sr,si,Mk,Kk,bk,Vk,Rec,var_para)


s0 = [sr,si];
s = getScontour(s0,Nz);
sPst = s(Nz+1:end);
% Finding the set of low-frequencies
[course_idS] = find(abs(imag(sPst(:)))<=var_para.si_max);%Indices for low-freq.
s_LF = sPst(course_idS);% set of low-freq.
% Get the right hand side
RickFDq = eval(var_para.source(var_para.Amp,var_para.width,var_para.t0,s_LF)).';
% Solve the ROMs 
ur_h = zeros(size(Vk,1)/1,length(s_LF));
parfor sid = 1:length(s_LF)
    ur = (s_LF(sid)^2*Mk + Kk)\(RickFDq(sid,1)*bk);    
    ur_h(:,sid) = Vk*ur;
end
% Seismogram
Snapshot_WeekRBPst_LF = Rec(:,1)'*ur_h;
% Seismo reconstrucion
Snapshot_WeekRBPst = zeros(1,Nz);
Snapshot_WeekRBPst(:,course_idS) = Snapshot_WeekRBPst_LF; %Filtered version of the solutions
Snapshot_WeeRBkNeg = fliplr(real(Snapshot_WeekRBPst(:,:)) - 1i*imag(Snapshot_WeekRBPst(:,:))); % Complex conugate symmetry  
seismoLD_RBM = [Snapshot_WeeRBkNeg,Snapshot_WeekRBPst]; %         

ap = Getap(seismoLD_RBM, Nz, si);
ap = ap/norm(abs(ap),inf);
end

function Lag_Coef = Getap(F, Nz, si)

FFTSamples = zeros(1,2*Nz,'double');            % Twice the samples as the number of coefficients

jdx = -Nz:(Nz-1);
theta_jhalf = (jdx+1/2)*pi/Nz;
Wtemp = exp(1i*theta_jhalf);
FFTSamples(jdx+Nz+1) = (2*si./(1-Wtemp)).*F;
%Note the order: FFTSamples(1,1:2*N)
TempCoef = fftshift(fft(fftshift(FFTSamples)))/(2*Nz); 
%Use only part of the TempCoef
Lag_Coef = TempCoef.*exp(-1i*(-Nz:Nz-1)*pi/(2*Nz));

end


function s = getScontour(s0,Nz)
jdx = -Nz:(Nz-1);
theta_jhalf = (jdx+1/2)*pi/Nz;
cstExp = exp(1i*theta_jhalf);
s = s0(1)-s0(2)*(cstExp+1)./(cstExp-1); 
end
