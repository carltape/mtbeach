function hdur = m02hdur(M0)
%M02HDUR input seismic moment in N-m and output half-duration in seconds

M0 = M0(:);
%K = 2.4e-6;                 % Ekstrom, 1996: Dahlen and Tromp (1999), p. 178
K = 1.05e-8 * 10^(7/3);     % Ekstrom et al. (2012), Eq. 1, converted for N-m

hdur = K * M0.^(1/3);

%--------------------------------------------------------------------------

if 0==1
    Mw = linspace(2,5,100);
    M0 = mw2m0(1,Mw);
    hdur = m02hdur(M0);
    
    figure; nr=2; nc=1;
    
    subplot(nr,nc,1); plot(Mw,hdur,'b'); grid on;
    xlabel('Moment magnitude, Mw');
    ylabel('Empirical half-duration, hdur, s');
    
    % SPECFEM2D: f0 = 1/hdur
    subplot(nr,nc,2); plot(Mw,1./hdur,'b'); grid on;
    xlabel('Moment magnitude, Mw');
    ylabel('Empirical central frequency, f0, Hz');
end

%==========================================================================