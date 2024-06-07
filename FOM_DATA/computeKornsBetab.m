function C_korn = computeKornsBetab(HFM)

H = HFM.M + HFM.H;
OPTS.tol    = 1e-6;
OPTS.disp   = 0;
OPTS.maxit  = 500;

C_korn = eigs(0.5*( HFM.K + HFM.K'),0.5*( H + H'), 1, 'sm', OPTS);   

end