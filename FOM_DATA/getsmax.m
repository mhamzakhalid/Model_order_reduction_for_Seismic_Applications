function [smax,RICK_smax,s_v] = getsmax(width,tolerance,fdualnorm,ldualnorm,K_HFM,M_HFM,Amp,C_korn)
OPTS.tol    = 1e-6;
OPTS.disp   = 0;

sr = 0.35;
si = linspace(width,50,50);
s_v = sr + 1i*si;
a = width;
%Amp = 1;
t0 = 4*pi/width;
RICK = zeros(length(s_v),1);
for sid = 1:length(s_v)
    s = s_v(sid);
    %RICKLF(sid,1) = getRHSnoSplit(var,s);
    RICK(sid,1) = Amp*exp(-a^2*t0^2/4)/a^3*(  a^3*t0 +2*s*a - 2*s.^2*sqrt(pi).*erfc(vpa((-a^2*t0 + 2*s)/(2*a))).*exp(vpa((-a^2*t0 + 2*s).^2/(2*a).^2)) );
end
% Computing stability constant
beta_h = zeros(length(s_v),1);

%lambda0 = eigs(0.5*( K_HFM + K_HFM'),Xh, 1, 'sm', OPTS); 
lambda0 = C_korn/(1+C_korn);

for sid = 1:length(s_v)
    s = s_v(sid);
    beta_h(sid,1)  = GetbetaLB(s,lambda0);
end
                
smax_id = min(find(abs(fdualnorm*ldualnorm*RICK.*1./(beta_h))<=tolerance));
smax = abs(min(s_v(smax_id)));
RICK_smax = abs(RICK(smax_id)/norm(RICK,inf));
end

function betaLB =GetbetaLB(s,lam0)
tau = @(sr,si) 1 - (1 - (sr^2 - si^2))/( (sr^2-si^2 - 1)^2 + 4*sr^2*si^2 );
sr = real(s); si = imag(s);
if lam0<=tau(sr,si) && tau(sr,si)<=1
    betaLB = abs(s^2*(1-tau(sr,si)) + tau(sr,si));
elseif lam0>tau(sr,si)
     betaLB = abs(s^2*(1-lam0) + lam0);
elseif tau(sr,si)>1
     betaLB = 1;
else
   fprint('Case not defined') 
end

end