
function  [ dqq, Eqq]  = offline_residual_comp(Aq,b_HFM,L_X,U_X,P_X,Q_X,V)
Qa = numel(Aq);


Nk = size(V,2);
Eqq = zeros(Nk,Nk,Qa,Qa);

C = cellfun(@(x) x*V, Aq,'UniformOutput',0);

Z = cellfun(@(x) Q_X*(U_X\(L_X\(P_X*x))), C,'UniformOutput',0);
dqq = cellfun(@(x) x'*b_HFM, Z,'UniformOutput',0);
dqq = horzcat(dqq{:});
Zt = cellfun(@ctranspose,Z,'UniformOutput',0);
for q1 = 1 : Qa
        T = cellfun(@(x) Zt{q1}*x, C,'UniformOutput',0);
        Tnew = cat(4, T{:});
        Eqq(:,:,q1,:) = permute(Tnew, [1, 2, 4, 3]);
end

end