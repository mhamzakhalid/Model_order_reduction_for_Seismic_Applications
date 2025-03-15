function deltak_LF = evalEstimateovers_fixedm(s_LF,uk,dqq,Eqq,res_ff,THETA_F_Vec,beta_lb)
Ns_LF = length(s_LF);
 
 % size 2 x Nz
THETA_A_MAT = zeros(2,Ns_LF);
THETA_A_MAT(1,:) = s_LF.^2;
THETA_A_MAT(2,:) = ones(1,Ns_LF);
% 2 real(sum_i=1^2  theta_a*conj(theta_f)*uN'*d_i)     
res_af = real(sum((THETA_A_MAT'.*THETA_F_Vec.') .*(uk(:,:,1)'*dqq)  ,2)).';             
% 2 real(sum_i=1^2  theta_a*theta_a'*uN'*e_ij*uN)             
uk_4d = repmat(uk(:,:,1),[1,1,2,2]);
THETA_AtimesTheta_A = bsxfun(@times, conj(permute(THETA_A_MAT, [1 3 2])), permute(THETA_A_MAT, [3 1 2]));
UmE_matU = permute(pagemtimes(pagectranspose(uk_4d),pagemtimes(Eqq,uk_4d)), [3 4 1 2]);
%res_aa2 = getT2(THETA_AtimesTheta_A,UmE_matU,Ns_LF);
res_aa = zeros(1,Ns_LF);
for sid = 1:Ns_LF        
    res_aa(1,sid) = sum(THETA_AtimesTheta_A(:,:,sid).*UmE_matU(:,:,sid,sid) ,'all');        
end
res = res_aa - 2*real(res_af) + res_ff;

deltak_LF = sqrt(abs(res))./min(beta_lb);

end
