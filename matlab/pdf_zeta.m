function [zeta,fzeta] = pdf_zeta(bfixedphi,n)
%PDF_ZETA probability density function for zeta angle (radians) for uniform moment tensors
%
% INPUT
%   fixedphi    =false for all moment tensors; =true for fixed phi
% optional
%   n           number of points
%
% OUTPUT
%   zeta        crack fraction within CDC model, in RADIANS [0,pi/2]
%   fzeta       f(zeta)
%
% See TapeTape2013 "The classical model for moment tensors"
% called by plotMT_phizeta.m
%
% Carl Tape, 2022-09-09
%

if nargin==1, n = 500; end

zeta = linspace(0,pi/2,n);
cosz = cos(zeta);
sinz = sin(zeta);

if bfixedphi
    % UniformParamForFobeChangesTOSS.pdf
    fzeta = 3/2 * cosz.^3;
else
    % PDFsForZetaAndTheta.pdf, AlphaZetaThetaForCarloTOSS.pdf
    fzeta = 4 * cosz.^3 .* sinz;
end
