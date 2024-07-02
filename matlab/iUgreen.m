%
% function igreen = iUgreen(v)
% Carl Tape, 11-March-2011
%
% See Latex notes on beachballs.
% 
% calls xxx
% called by Ueigvec.m
%

function igreen = iUgreen(v)

vx = v(1);
vy = v(2);
vz = v(3);
vmag = sqrt(vx^2 + vy^2 + vz^2);

EPSVAL = 1e-12;

igreen = NaN;
if all(abs([vx vy vz]) <= EPSVAL)
    % center of sphere
    igreen = 0;
    
elseif and( and(abs(vy) <= EPSVAL, abs(vz) <= EPSVAL), and(vx > EPSVAL, vx < 1-EPSVAL) )
    % radial segment from center to (1,0,0)
    igreen = 1;
    
elseif and( and(abs(vx) <= EPSVAL, abs(vz) <= EPSVAL), and(vy > EPSVAL, vy < 1-EPSVAL) )
    % radial segment from center to (0,1,0)
    igreen = 2;
    
elseif and( and(abs(vx) <= EPSVAL, abs(vy) <= EPSVAL), and(vz > EPSVAL, vz < 1-EPSVAL) )
    % radial segment from center, to (0,0,1)
    igreen = 3;
    
elseif and( all([vx vy vz] > EPSVAL), abs(vmag - 1) <= EPSVAL )
    % 1/8 sphere patch on surface of sphere
    igreen = 4;
    
elseif and( and(all([vx vz] > EPSVAL), vy < -EPSVAL), abs(vmag - 1) <= EPSVAL )
    % 1/8 sphere patch on surface of sphere
    igreen = 5;
    
elseif and( and(all([vy vz] > EPSVAL), abs(vx) <= EPSVAL), vmag < 1-EPSVAL )
    % 1/4 circle 
    igreen = 6;
    
elseif and( all([vx vz] > EPSVAL), vmag < 1-EPSVAL )
    % interior of G
    igreen = 7;
end

%==========================================================================
