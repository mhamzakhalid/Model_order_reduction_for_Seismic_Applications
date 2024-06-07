function [seismo_x,seismo_y,u] = solveTimeNewmark(HFM,var)


ft = var.ft(var,var.t,var.t0);
n = HFM.DOF; 
dt = var.dt; % time-step
K = HFM.K;
M = HFM.M; 
b = HFM.f; % matrices and point source vector 

u_old = zeros(n,1);
v_old = zeros(n,1);
Ntime = size(var.t,2);
gamma = 1/2; % Newmark-Beta constants
beta = 1/4; % Newmark-Beta constants (1/4 for implicit)

seismo_x = zeros(size(HFM.Rec,2),length(var.t));
seismo_y = zeros(size(HFM.Rec,2),length(var.t));
MK = M + beta*dt^2*K;        
[L,U,P,Q] = lu(MK);
a_old0 = L\(P*(b*ft(1)));
a_old1 = U\a_old0;
a_old = Q*a_old1;
u = zeros(n,length(var.t));
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
  u(:,i+1) = u_old;
  fprintf('Process: %2.4f percent \n',i/(Ntime-1)*100)

end
