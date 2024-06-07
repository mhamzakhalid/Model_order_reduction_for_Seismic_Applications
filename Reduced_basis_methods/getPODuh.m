function uh = getPODuh(Aq,bi_HFM,s_LF,RICKLF,lam_rand,mu_rand)
K_splitHFM = cellfun(@(lam, mu, Kl, Ku) lam * Kl + mu * Ku, num2cell(lam_rand(1,:)), num2cell(mu_rand(1,:)), Aq(2:6).', Aq(7:11).', 'UniformOutput', false);
K_HFM = K_splitHFM{1} + K_splitHFM{2} + K_splitHFM{3} + K_splitHFM{4} + K_splitHFM{5};
Ns_LF = length(s_LF);
uh = zeros(size(bi_HFM,1),Ns_LF);
M_HFM = Aq{1};
parfor sid=1:Ns_LF
    uh(:,sid) = (s_LF(sid)^2*M_HFM + K_HFM)\(RICKLF(sid)*bi_HFM);
end

end