

function beta_i = getBetaLB(Ki_HFM_lam,Ki_HFM_mu,M_HFM,lam,mu,sPst,Caseid,Xh,var,C_korn,maxLam,minLam,minmu,maxmu)

OPTS.tol    = 1e-2;
OPTS.disp   = 0;
OPTS.maxit  = 100;

%% Lower bound
beta_i = zeros(1,length(sPst));
c1 = min([1,1./max(minLam(:)./var.lami(:)),1./max(minmu(:)./var.mui(:))]);
c2 = max([1,max(maxLam(:)./var.lami(:)),max(maxmu(:)./var.mui(:))]);
c_equivalance = c1/c2;


switch Caseid
            case '1'
                K_split = cell(5,1);
                for lid = 1:5
                K_split{lid} = lam(lid)*Ki_HFM_lam{lid} + mu(lid)*Ki_HFM_mu{lid};
                end
                K_full = K_split{1} + K_split{2} + K_split{3} + K_split{4} + K_split{5};               
                lambdaa = eigs(0.5*( K_full + K_full'),Xh, 1, 'sm', OPTS);  
                for sid = 1:length(sPst)
                    s = sPst(sid);
                    beta_i(1,sid)  = c_equivalance*GetbetaLB(s,lambdaa);
                end
                fprintf('Running Beta_LB(s,m)\n')
            case '2'         
                K_split = cell(5,1);
                for lid = 1:5
                K_split{lid} = lam(lid)*Ki_HFM_lam{lid} + mu(lid)*Ki_HFM_mu{lid};
                end
                K_full = K_split{1} + K_split{2} + K_split{3} + K_split{4} + K_split{5};
                Xh = K_full + M_HFM;
                % Computing stability constant
                parfor (sid = 1:length(sPst))
                    s = sPst(sid);
                    Ah = K_full + s^2*M_HFM;
                    beta_i(1,sid)  = get_betah(Xh,Ah);
                end
          case '3'
            lambda0 = C_korn/(1+C_korn);

            for sid = 1:length(sPst)
                s = sPst(sid);
                beta_i(1,sid)  = c_equivalance*GetbetaLB(s,lambda0);
    
            end      
            case '4'         

        for sid = 1:length(sPst)
            sr = real(sPst(sid)); si = imag(sPst(sid));
            beta_i(1,sid)  = c_equivalance*min([1,sr^2,2*sr/sqrt(4*sr^2+si^2)]);

        end      
end

end
function betaLB =GetbetaLB(s,lam)
tau = @(sr,si) 1 - (1 - (sr^2 - si^2))/( (sr^2-si^2 - 1)^2 + 4*sr^2*si^2 );
sr = real(s); si = imag(s);
if lam<=tau(sr,si) && tau(sr,si)<=1
    betaLB = abs(s^2*(1-tau(sr,si)) + tau(sr,si));
elseif lam>tau(sr,si)
     betaLB = abs(s^2*(1-lam) + lam);
elseif tau(sr,si)>1
     betaLB = 1;
else
   fprint('Case not defined') 
end

end

function [beta] = get_betah(Xh,Ah)

% Options for s^2 M+K   
OPTS.tol    = 6e-3;
OPTS.issym  = 0;
OPTS.disp   = 1;
OPTS.maxit  = 200;
  
[D]     = eigs(@(x)Linsolve(x,Ah,Xh), size(Ah,1), Xh, 1, 'sm', OPTS);
beta    = sqrt(abs(D));


end

function y = Linsolve(x, Ah, Xh)
% y = Ah \ ( Xh * (Ah' \ x) )
y1 = Ah'\x;
y2 = Xh*y1;
y   = Ah\y2;

end

