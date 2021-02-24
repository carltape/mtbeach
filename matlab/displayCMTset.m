function [M0,Mw,M0tot,Mwtot] = displayCMTset(M,tshift,hdur,zele_m)
%DISPLAYCMTSET display information about a set of moment tensors
%
% It was written with finite-source CMTSOLUTION files in mind.
%
% calls CMT2m0.m, m02mw.m
% called by read_CMTSOLUTION_finite.m
%

if nargin==1, tshift = []; hdur = []; zele_m = []; end

% compute M0 and Mw for each subsource
M0 = CMT2m0(1,M);
Mw = m02mw(1,M0);
M0tot = sum(M0);
Mwtot = m02mw(1,M0tot);
    
% number of subsources
n = length(Mw);

disp('----------------------------------');
disp('summary from displayCMTset.m');
%disp(sprintf('file: %s',fname));
disp(sprintf('number of subsources: %i',n));
disp(sprintf('min/median/max Mw of %i subsources: %.2f / %.2f / %.2f',...
    n,min(Mw),median(Mw),max(Mw)));
disp(sprintf('total M0: %.3e N-m',M0tot));
disp(sprintf('total Mw: %.3f',Mwtot));
disp(sprintf('min/median/max elevation of %i subsources:',n));
disp(sprintf('    %.3f km / %.3f km / %.3f km',...
    min(zele_m)/1000,median(zele_m)/1000,max(zele_m)/1000));
disp(sprintf('min/median/max tshift of %i subsources:',n));
disp(sprintf('         %.1f s / %.1f s / %.1f s',...
    min(tshift),median(tshift),max(tshift)));
disp(sprintf('      %.2f min / %.2f min / %.2f min',...
    min(tshift)/60,median(tshift)/60,max(tshift)/60));
disp(sprintf('duration of rupture of %i subsources:',n));
disp(sprintf('        %5.2f s',max(tshift)-min(tshift)));
disp(sprintf('        %5.2f min',(max(tshift)-min(tshift))/60 ));
disp(sprintf('min/median/max hdur of %i subsources:',n));
disp(sprintf('         %.1f s / %.1f s / %.1f s',min(hdur),median(hdur),max(hdur)));
disp('----------------------------------');

%==========================================================================
