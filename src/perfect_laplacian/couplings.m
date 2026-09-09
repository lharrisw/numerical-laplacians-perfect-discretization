function c = couplings(r1,r2,varargin)
%COUPLINGS Returns couplings of the perfect Laplacian
%   C = COUPLINGS(R1,R2,N) returns the couplings/coefficients of the 
%   translationally-invariant perfect Laplacian. This is accomplished from
%   performing the renormalization group transformation for the discrete
%   action for a free scalar field. The arguments R1 and R2 refer to the
%   distance from the point (R1,R2) in the horizontal and vertical 
%   directions, respectively. The argument N determines the index of
%   summation for a double sum. The couplings are non normalized.
%
%   C = COUPLINGS(R1,R2,N) returns the coupling constant for position
%   (R1,R2) with 4*N + 2 summation terms.
%
%   C = COUPLINGS(R1,R2,N1,N) returns the coupling constant for the perfect
%   Laplacian near a wall, where N1 is the distance between (R1,R2) and the
%   wall with zero boundary conditions.
%
%   C = COUPLINGS(R1,R2,N1,N2,N) returns the coupling constant for the
%   perfect Laplacian near a corner, where N1 is the number of units 
%   between (R1,R2) and the positive x-axis and N2 is the number of units
%   between (R1,R2) and the positive y-axis with N specifying the indices
%   of summation.
%
%   Example: Compute the coupling for (r1,r2) = (-2,1) of the
%            translationally-invariant perfect Laplacian that is one unit
%            away from the wall with N = 10.
%
%       c = couplings(-2,1,10);

% Copyright 2020 Latimer De'Shone Harris-Ward.

% translationally invariant perfect laplacian
rho = @(r1,r2,n) integral2(@(q1,q2)exp(1i*(r1*q1+r2*q2))...
       ./const(q1,q2,n),-pi,pi,-pi,pi)/(2*pi).^2;

% perfect laplacian near a wall; Hauswirth Equation (2.25)
rhow = @(r1,r2,n1,n) rho(r1,r2,n)-rho(r1+2*n1+1,r2,n);
   
% perfect laplacian near a corner; Hauswirth Equation (2.26)
rhoc = @(r1,r2,n1,n2,n) rho(r1,r2,n)-rho(r1+2*n1+1,r2,n)- ...
                        rho(r1,r2+2*n2+1,n)+rho(r1+2*n1+1,r2+2*n2+1,n);

if nargin < 3
    error('Not enough input arguments.')
end

if nargin == 3
    % COUPLINGS(R1,R2,N)
    N = varargin{1};
    if N <= 0
        error('The index of summation must be a positive integer.')
    end
    c = real(rho(r1,r2,N));
elseif nargin == 4
    % COUPLINGS(R1,R2,N1,N)
    n1 = varargin{1};
    N = varargin{2};
    if n1 < 0
        error('The distance must be greater than or equal to zero.')
    end
    if N <= 0
        error('The index of summation must be a positive integer.')
    end
    c = real(rhow(r1,r2,n1,N));
elseif nargin == 5
    % COUPLINGS(R1,R2,N1,N2,N)
    n1 = varargin{1};
    n2 = varargin{2}; 
    N = varargin{3};
    if n1 < 0 || n2 < 0 || (n1 < 0 && n2 < 0)
        error('The distances must be greater than zero.')
    end
    if N <= 0
        error('The index of summation must be a positive integer.')
    end
    c = real(rhoc(r1,r2,n1,n2,N));
end
   
% function for 1/rho(q)
function rhoinv = const(q1,q2,n)

k = 2; % blocking kernel number
rhoinv = 0; % initialization

% represents equation 2.21 in Hauswirth, 17 in Hasenfratz, and the system
% 37 in Hasenfratz 2
for s = -n:n
    for t = -n:n
        rhoinv = rhoinv + 1./(((q1+2*pi*s).^2+(q2+2*pi*t).^2)...
                          .*(q1+2*pi*s).^2.*(q2+2*pi*t).^2);
    end
end
rhoinv = 1/(3*k)+16*sin(q1/2).^2.*sin(q2/2).^2.*rhoinv;
end
end
