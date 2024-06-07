%% This function files uses POD method on a given snapshot matrix to find ROM
%%% Inputs
% Xh: Scalar product matrix 
% tol: Maximum tolerance to select singular values/vectors for the POD
% basis
% Nkmax: Maximum number of Reduced basis
%%% Outout
%ROM.V : POD basis functions
%ROM.Sigma : Singular values 

function ROM = POD(Xh,Snapshot,tol,Nkmax)
Ns = size(Snapshot,2);
% Solving the Eigenvalue problem by using SVD-based approach
C = Snapshot'*(Xh*Snapshot)*1/Ns;
[phi, Sigma]  = svd(full(C),'econ');
Sigma         = diag(Sigma);                                               
% Extract relevant basis     
Sigma_N = cumsum(Sigma);
if ~isempty(Nkmax)
     N = Nkmax;
else
     N  = find(tol^2>=1-Sigma_N/sum(Sigma),1,'first');
end 
phi   =  phi(:,1:N);
% Assemble POD basis
phi   =  Snapshot*phi;    
for cc = 1 : N
    phi(:,cc) = 1/(sqrt(Ns*Sigma(cc))) * phi(:,cc);
end 
% Output
ROM.V = phi;
ROM.Sigma = sqrt(Sigma);

end