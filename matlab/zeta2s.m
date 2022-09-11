function s = zeta2s(zeta)
%ZETA2 convert zeta angle to s
%
% INPUT
%   zeta    zeta angle, in RADIANS [0,pi/2]
%
% OUTPUT
%   s       parameter [0,1]
%
% This is the reverse function of s2zeta.m
% See examples in s2zeta.m
%
% Carl Tape, 2022-09-09
%

s = 1/8 * (9 * sin(zeta) + sin(3*zeta));
