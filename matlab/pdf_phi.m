function [phi,fphi] = pdf_phi(n)
%PDF_PHI probability density function for phi angle (radians) for uniform moment tensors
%
% INPUT
% optional
%   n       number of points
%
% OUTPUT
%   phi     azimuthal angle on the lune, RADIANS [-pi,pi]
%   fphi    f(phi)
%
% See TapeTape2013 "The classical model for moment tensors"
% called by plotMT_phizeta.m
%
% Carl Tape, 2022-09-09
%

if nargin==0, n = 501; end

phi = linspace(-pi,pi,n);
% AlphaZetaThetaForCarloTOSS.pdf
fphi = 2 ./ (5*pi - 3*pi*cos(2*phi));
