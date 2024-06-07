%%
function getCW()
    addpath('FOM_DATA')
    load('optimals0_Ricker.mat')
    var.s0 = s_optROM_Ricker;
    wr = var.s0(1);
    wi = var.s0(2);
    Nz = 608;
    T = 20;
    Lap2 = zeros(Nz,Nz);
    parfor p_id = 0:Nz-1  
        [Lap2(p_id+1,:)] = SolveintegralCw(p_id,Nz,wr,wi,T); 
    end
    Lag_intC2 = (Lap2' - diag(diag(Lap2))) + Lap2;
    CW = wi/Nz*sqrt(sum(abs(Lag_intC2),'all'));    
    save('CW.mat','CW')   
    rmpath('FOM_DATA')
    end    


function [Lap_p,resi] = SolveintegralCw(p_id,Nz,sr,si,T)
Lap_p = zeros(1,Nz);
resi = zeros(1,Nz);
 for q_id = p_id+1:Nz  
       integrand =@(t)  exp(2*(sr-si)*t).*Laguerre_p(2*si*t,p_id).*Laguerre_p(2*si*t,q_id-1);
       [Lap_p(1,q_id),resi(1,q_id)] = quadgk(integrand,0,T,'AbsTol',1e-8,'RelTol',1e-10);

 end   

end
function Lag_p = Laguerre_p(x, p)
    Nx = numel(x);  
    Lag = zeros(p+1,Nx);
    Lag(1,:) = ones(1,Nx);
    if  1 <= p 
      Lag(2,:) = ones(1,Nx) - x;
    end
    for k_id = 2 : p
        numerator =  ( 2 * k_id - 1 - x ) .* Lag(k_id,:) - ( k_id - 1 ) .* Lag(k_id-1,:);
        Lag(k_id+1,:) = numerator./k_id;
    end
    Lag_p = Lag(p+1,:);
end


