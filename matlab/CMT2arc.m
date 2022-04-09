function [thetadiff] = CMT2arc(M1,M2,bfigure)
%CMT2ARC compute the arc distance between two moment tensor source types
%
% This angle will be equal to the matrix angle between two moment tensors
% having the same frame (basis).
%
% INPUT
%   M1,M2     6 x n matrices of moment tensors
%
% calls CMTdecom.m (plot_histo.m, lam2lune.m for plotting(
%
% Carl Tape, 2022-04-08
%

if nargin==2, bfigure = false; end

% check sizes
[a,b] = size(M1); if a~=6, M1 = M1'; n1 = a; else n1 = b; end
[c,d] = size(M2); if c~=6, M2 = M2'; n2 = c; else n2 = d; end
if n1~=n2, error('M1 and M2 must have same number of columns'); else n = n1; end

lam1 = CMTdecom(M1);
lam2 = CMTdecom(M2);

norm1 = zeros(1,n);
norm2 = zeros(1,n);
for ii=1:n
    norm1(ii) = norm(lam1(:,ii));
    norm2(ii) = norm(lam2(:,ii));
end

cosx = dot(lam1,lam2) ./ (norm1 .* norm2);
% correction for possible numerical errors
% this correction is needed for comparing U's that are very close to each other
ipos = cosx > 1;
ineg = cosx < -1;
disp(sprintf('%i/%i dot products > 1',sum(ipos),n));
disp(sprintf('%i/%i dot products < 1',sum(ineg),n));
cosx(ipos) = 1;
cosx(ineg) = -1;
% angle
thetadiff = 180/pi * acos(cosx);

if bfigure
    figure; plot_histo(thetadiff,[0:5:180]);

    [gamma1,delta1] = lam2lune(lam1);
    [gamma2,delta2] = lam2lune(lam2);
    if n < 20
        figure; hold on;
        for ii=1:n
            plot(gamma1(ii),delta1(ii),'ro');
            plot(gamma2(ii),delta2(ii),'bo');
            plot([gamma1(ii) gamma2(ii)],[delta1(ii) delta2(ii)],'k');
            text(gamma2(ii),delta2(ii),sprintf('%.0f',thetadiff(ii)));
        end
        axis equal, axis([-30 30 -90 90]); grid on;
    end
end

%==========================================================================

if 0==1
    clear, close all, clc
    % uniform moment tensors (fixed magnitude)
    n = 10000;
    [M1,v,w,kappa,sigma,h,lam] = uniformMT(n);
    [M2,v,w,kappa,sigma,h,lam] = uniformMT(n);
    
    thetadiff = CMT2arc(M1,M2,true);
end

%==========================================================================