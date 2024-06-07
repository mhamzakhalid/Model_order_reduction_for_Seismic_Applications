% We tailor the code to our problem. 
% Copyright (c) 2011, Patrick Kano, Moysey Brio
% All rights reserved.
%  The function [so,bo] = wpar2(F,t,N,sig0,sigmax,bmax,Ntol) estimates 
%  the optimal parameters s (= sigma) and b for the Weeks method.  
%  The function requires no information on the singularities of the transform.
%
%  Input:
%  F      -- Transform; symbolic function of s. (We replace this with the
%  reduced order solutions)
%  t      -- Ordinate where f(t) are to be computed (scalar).
%  N      -- Number of terms in the Weeks expansion.
%  sig0   -- Laplace convergence asbscissa.
%  sigmax -- Maximum possible value of sigma.
%  bmax   -- Maximum possible value of b.
%  Ntol   -- Determines the tolerance in the optimization routine fmin.
%            Large Ntol => small tolerance.  Recommended: Ntol = 20 to 50.
%  Output:
%  (so, bo) -- Estimated parameters to be used in weeks.m or weekse.m.


function [wr,wi] = wpar2(F, t, N, sig0, sigmax, bmax, Ntol,Mk,Kk,bk,V,Rec,var)

tols = (sigmax-sig0)/Ntol; 
tolb = bmax/Ntol;
options = optimset('fminbnd');
sooptnew = optimset(options,'TolX',tols,'Display','iter');
booptnew = optimset(options,'TolX',tolb,'Display','iter');
wr = fminbnd('werr2e',sig0,sigmax,sooptnew,F,t,N,sig0,sigmax,bmax,tolb,Mk,Kk,bk,V,Rec,var);
wi = fminbnd('werr2t',0,bmax,booptnew,F,N,wr,Mk,Kk,bk,V,Rec,var); 

end