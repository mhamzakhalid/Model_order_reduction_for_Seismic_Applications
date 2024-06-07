function var_para = getvariables(HFM)
var_para.sigx = 0.008*10000;  % Gaussian width  x
var_para.sigy = 0.008*10000;  % Gaussian width  y
var_para.sigmx = 1e1; % For scaling of Gaussian 
var_para.sig0 = 1e1; % For scaling of Gaussian 

var_para.sigm = 1e1;
var_para.sigmy = 1e1; 
% Parameters Sinc time signal    
var_para.Amp = 1e5; % Amplitude (Scaling factor)

var_para.gx = @(x,y,x0,y0,sig0,sigx,sigy) 1/(2*pi*sig0^2)*...
                                           exp(-1*( 1/(2*sigx^2)*(x - x0).^2  + 1/(2*sigy^2)*(y - y0).^2 ));




% Reading material properties
% - id layer
% - 2 rho/density (kg/m3)
% - 3 vs (km/s) Shear-wave velocity
% - 4 vp (km/s) P-Wave velocity
% - 5 mu/G/shear modulus
% - 6 kappa/bulk modulus
% - 7 lambda(lame parameter)
var_para.kappa = cell2mat(HFM.materProp(:,6));
var_para.vp = cell2mat(HFM.materProp(:,4));
var_para.vs = cell2mat(HFM.materProp(:,3));

var_para.mui  = cell2mat(HFM.materProp(:,5));
var_para.lami = cell2mat(HFM.materProp(:,7));
var_para.rhoi = cell2mat(HFM.materProp(:,2));

var_para.ft = @(var,t,t0) var.Amp*(1 - 1/2*var.width^2*(t - t0).^2).*exp(-var.width^2 *(t - t0).^2/4);
var_para.source =@(Amp,a,t0,s) Amp*exp(-a^2*t0^2/4)/a^3*(  a^3*t0 +2*s*a - 2*s.^2*sqrt(pi).*erfc(vpa((-a^2*t0 + 2*s)/(2*a))).*exp(vpa((-a^2*t0 + 2*s).^2/(2*a).^2)) );

%% Time stepping 
var_para.dt = 1e-03; % delta t
var_para.t = 0:var_para.dt:20;% Time steps 

%% Receivers
var_para.x0r = 9.475197699000000e+03;%9.4752e+03;
var_para.y0r = 20000;
fprintf('Receiver location (x_r,y_r)=(%2.1e,%2.1e)\n',var_para.x0r,var_para.y0r)
%% Source locations
var_para.x0s = 9.915810261499999e+03;
var_para.y0s = 1.641749665800000e+04;
fprintf('Source location (x_0,y_0)=(%2.1e,%2.1e)\n',var_para.x0s,var_para.y0s)

end