% We tailor the code to our problem. 
% Copyright (c) 2011, Patrick Kano, Moysey Brio
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function E = werr2e(sig,F,t,N,sig0,sigmax,bmax,tolb,Mk,Kk,bk,Vk,Rec,var)

%   The function E = errorb(sig,F,t,N,sig0,sigmax,cc,tolb) computes
%   the error bound E = truncation plus conditioning error
%   on the optimal curve b = b(sig).

options = optimset('fminbnd');
optnew = optimset(options,'TolX',tolb);
b  = fminbnd('werr2t',sig0,bmax,optnew,F,N,sig,Mk,Kk,bk,Vk,Rec,var);  % Estimate optimal b=bopt(sig)

M  = 2*N;
ap  = wcoef(F,M,sig,b,Mk,Kk,bk,Vk,Rec,var);
a1 = ap(2*N+1:3*N); sa1 = sum(abs(a1));
a2 = ap(3*N+1:4*N); sa2 = sum(abs(a2));
E  = exp(sig*t)*(sa2+eps*sa1);
E  = log(E);



end