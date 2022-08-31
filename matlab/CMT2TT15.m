function [rho,v,w,kappa,sigma,h,K,N,S,thetadc,lam,U] = CMT2TT15(M,bdisplay)
%CMT2TT converts a moment tensor to six parameters of TapeTape2015
%
% INPUT
%   M           6 x n moment tensors in CMT convention (UP-SOUTH-EAST)
%               M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%   bdisplay    OPTIONAL (if present, will display details)
%               
% OUTPUT
%   rho         norm of the moment tensor (note: rho = sqrt(2)*M0)
%   v           coordinate proportional to lune longitude gamma
%   w           coordinate proportional to lune latitude delta
%   kappa       strike angle, degrees: [0,360]
%   sigma       slip (or rake) angle, degrees
%   h           cos(dip): [0, 1]
% optional:
%   K           strike vector (SOUTH-EAST-UP)
%   N           normal vector (SOUTH-EAST-UP)
%   S           slip vector (SOUTH-EAST-UP)
%   thetadc     angle to DC
%   lam         eigenvalues
%   U           basis (SOUTH-EAST-UP)
%
% Reverse program for TT152CMT.m
% See TapeTape2015 "A geometric setting for moment tensors"
%
% The one modification is that here we return w = 3*pi/8 - u, rather than
% the u coordinate in TT2015. w is analagous to lune latitude delta,
% whereas u is analagous to lune colatitude beta
%
% Carl Tape, 2015-08-07
%

if nargin==1, bdisplay=false; end

[gamma,delta,M0,kappa,theta,sigma,K,N,S,thetadc,lam,U] = CMT2TT(M,bdisplay);

rho = sqrt(2)*M0;
[v,w] = lune2rect(gamma,delta);
h = cos(theta*pi/180);
