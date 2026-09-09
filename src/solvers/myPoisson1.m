function uperf = myPoisson1(f,h,varargin)
%MYPOISSON Solves the Laplace or Poisson equation on a square.
%   UPERF = MYPOISSON1(F,H,...) is an algorithm that solves the Laplace 
%   or Poisson equation on the square using the five-point stencil, the 
%   modified nine-point stencil, or the perfect modified 9-point stencil.
%   The perfect modified nine-point stencil utilizes coefficients obtained
%   from the normalized translationally-invariant perfect Laplace operator.
%   The algorithm used is the vectorized Jacobi's algorithm.
%
%   The matrix F is a 2-D matrix of size (N + 2) x (N + 2), where N >= 3,
%   and is the right-hand side of the Poisson equation. The value H is the
%   step size given by L/(N + 1), where L is the side length of the square.
%   The matrix U is the same size as F and contains the boundary
%   information.
%
%   UPERF = MYPOISSON1(F,H) solves the Poisson equation with the default
%   setting of the 5-point Laplacian using 100 iterations for homogeneous
%   Dirichlet boundary conditions.
%
%   UPERF = MYPOISSON1(F,H,OPTIONS) solves the Poisson equation with the
%   OPTIONS specified. The OPTIONS that UPERF accepts are: MaxIter and
%   Stencil. The boundary conditions are homogeneous Dirichlet conditions.
% 
%   UPERF = MYPOISSON1(F,H,U) solves the Laplace or the Poisson equation
%   with specified boundary conditions in U, the initialization matrix. The
%   default setting is the 5-point Laplacian using 100 iterations.
%
%   UPERF = MYPOISSON1(F,H,U,OPTIONS) solves the Laplace or Poisson
%   equation with the specified boundary conditions and OPTIONS:
%   
%       Option value:       Output:
%               5               5-point stencil, 100 iterations
%              '9m'             modified 9-point stencil, 100 iterations
%              '9pm'            perfect modified 9-point stencil, 
%                               100 iterations
%
%   Example: Solve Poisson's equation on the square with homogeneous
%            boundary conditions with a right-hand side of 
%            f(x,y) = 2*(x*(x-a)+y*(y-a)), 0 <= x,y <= a, for a = 1
%            using the modified 9-point stencil and 150 iterations.
%
%       N = 25;
%       maxiter = 150;
%       a = 1;
%       h = a/(N + 1);
%       [X,Y] = meshgrid(0:h:a,0:h:a);
%       f = 2*(X.*(X-a)+Y.*(Y-a));
%       u = myPoisson1(f,h,'Maxiter',150,'Stencil','9m');
%       mesh(X,Y,f)
%
%   See also MYPOISSON2, MYPOISSON3, PERFPOISSON.

%  Copyright 2020 Latimer De'Shone Harris-Ward.        

% coefficients for the 9-point perfect Laplacian
p = 0.618835977336115;
q = 0.190582011331943;
c = 3.237671954672230;

% defining parameters
N = length(f(1,:)) - 2;
I = 2:N+1;
J = 2:N+1;
uperf = zeros(N+2,N+2);

% checking the size of u and f
[m,n] = size(f);
if m ~= n || m < n
    error('The matrix f must be square.')
end

% the default algorithm
if nargin == 2
    % UPERF(F,H) 
    % solution using 5 point standard/normalized perfect laplacian
    for iter = 1:100
        uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                           uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
    end
elseif nargin == 3
    % UPERF(F,H,U)
    uperf = varargin{1};
    if ~strcmpi('double',class(uperf))
        error("Input argument must be of type 'double.'")
    end
    [m,n] = size(uperf);
    if m ~= n || (m == 1 && n == 1)
        error('The input u must be a square matrix.')
    end
    for iter = 1:100
        uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                           uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
    end    
elseif nargin == 4
    % U(F,H,STRING,VALUE)
    % number of maximum iterations
    if strcmpi('MaxIter',varargin{1})
        maxiter = varargin{2};
        for iter = 1:maxiter
            uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                               uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
        end
    end
    % determining the stencil method
    if strcmpi('Stencil',varargin{1})
        if varargin{2} == 5
            for iter = 1:100
                uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                                   uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
            end
        elseif strcmpi('9m',varargin{2})
            for iter=1:100
                uperf(I,J) = 0.05*(4*uperf(I-1,J)+4*uperf(I+1,J)+ ...
                                   4*uperf(I,J-1)+4*uperf(I,J+1)+ ...
                                     uperf(I-1,J-1)+uperf(I-1,J+1)+ ... 
                                     uperf(I+1,J-1)+uperf(I+1,J+1)- ...
                                     6*h^2*(2*f(I,J)/3+(f(I+1,J)+...
                                     f(I-1,J)+f(I,J+1)+f(I,J-1))/12));
            end
        elseif strcmpi('9pm',varargin{2})
            for iter=1:100
                uperf(I,J) = (p*uperf(I-1,J)+p*uperf(I+1,J)+ ...
                              p*uperf(I,J-1)+p*uperf(I,J+1)+ ...
                              q*uperf(I-1,J-1)+q*uperf(I-1,J+1)+ ...
                              q*uperf(I+1,J-1)+q*uperf(I+1,J+1)- ...
                              h^2*(2*f(I,J)/3+(f(I+1,J)+f(I-1,J)+ ...
                              f(I,J+1)+f(I,J-1))/12))/c;
            end 
        else
            error("The only options are 5, '9m', or '9pm.'")
        end
    end
elseif nargin == 5
    % U(F,H,U,STRING,VALUE)
    % number of maximum iterations
    uperf = varargin{1};
    if ~strcmpi('double',class(uperf))
        error("Input argument must be of type 'double.'")
    end
    [m,n] = size(uperf);
    if m ~= n || (m == 1 && n == 1)
        error('The input u must be a square matrix.')
    end
    if strcmpi('MaxIter',varargin{2})
        maxiter = varargin{3};
        for iter = 1:maxiter
            uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                               uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
        end
    end
    % determining the stencil method
    if strcmpi('Stencil',varargin{2})
        if varargin{3} == 5
            for iter = 1:100
                uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                                   uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
            end
        elseif strcmpi('9m',varargin{3})
            for iter=1:100
                uperf(I,J) = 0.05*(4*uperf(I-1,J)+4*uperf(I+1,J)+ ...
                                   4*uperf(I,J-1)+4*uperf(I,J+1)+ ...
                                     uperf(I-1,J-1)+uperf(I-1,J+1)+ ... 
                                     uperf(I+1,J-1)+uperf(I+1,J+1)- ...
                                     6*h^2*(2*f(I,J)/3+(f(I+1,J)+...
                                     f(I-1,J)+f(I,J+1)+f(I,J-1))/12));
            end
        elseif strcmpi('9pm',varargin{3})
            for iter=1:100
                uperf(I,J) = (p*uperf(I-1,J)+p*uperf(I+1,J)+ ...
                              p*uperf(I,J-1)+p*uperf(I,J+1)+ ...
                              q*uperf(I-1,J-1)+q*uperf(I-1,J+1)+ ...
                              q*uperf(I+1,J-1)+q*uperf(I+1,J+1)- ...
                              h^2*(2*f(I,J)/3+(f(I+1,J)+f(I-1,J)+ ...
                              f(I,J+1)+f(I,J-1))/12))/c;
            end 
        else
            error("The only options are 5, '9m', or '9pm.'")
        end
    end
elseif nargin == 6
    % UPERF(F,H,STRING,VALUE,STRING,VALUE)
    if strcmpi('Maxiter',varargin{1})
        maxiter = varargin{2};
        if ~strcmpi('Stencil',varargin{3})
            error("The only second option is the string 'Stencil'.");
        end
        if varargin{4} == 5
            for iter = 1:maxiter
                uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                                   uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
            end
        elseif strcmpi('9m',varargin{4})
            for iter=1:maxiter
                uperf(I,J) = 0.05*(4*uperf(I-1,J)+4*uperf(I+1,J)+ ...
                                   4*uperf(I,J-1)+4*uperf(I,J+1)+ ...
                                     uperf(I-1,J-1)+uperf(I-1,J+1)+ ... 
                                     uperf(I+1,J-1)+uperf(I+1,J+1)- ...
                                     6*h^2*(2*f(I,J)/3+(f(I+1,J)+...
                                     f(I-1,J)+f(I,J+1)+f(I,J-1))/12));
            end            
        elseif strcmpi('9pm',varargin{4})
            for iter=1:maxiter
                uperf(I,J) = (p*uperf(I-1,J)+p*uperf(I+1,J)+ ...
                              p*uperf(I,J-1)+p*uperf(I,J+1)+ ...
                              q*uperf(I-1,J-1)+q*uperf(I-1,J+1)+ ...
                              q*uperf(I+1,J-1)+q*uperf(I+1,J+1)- ...
                              h^2*(2*f(I,J)/3+(f(I+1,J)+f(I-1,J)+ ...
                              f(I,J+1)+f(I,J-1))/12))/c;
            end 
        else
            error("The only options are 5, '9m', or '9pm.'")
        end
    elseif strcmpi('Stencil',varargin{1})
        maxiter = varargin{4};
        if ~strcmpi('MaxIter',varargin{3})
            error("The only second option is the string 'MaxIter'.");
        end
        if varargin{2} == 5
            for iter = 1:maxiter
                uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                                   uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
            end
        elseif strcmpi('9m',varargin{2})
            for iter=1:maxiter
                uperf(I,J) = 0.05*(4*uperf(I-1,J)+4*uperf(I+1,J)+ ...
                                   4*uperf(I,J-1)+4*uperf(I,J+1)+ ...
                                     uperf(I-1,J-1)+uperf(I-1,J+1)+ ... 
                                     uperf(I+1,J-1)+uperf(I+1,J+1)- ...
                                     6*h^2*(2*f(I,J)/3+(f(I+1,J)+...
                                     f(I-1,J)+f(I,J+1)+f(I,J-1))/12));
            end
        elseif strcmpi('9pm',varargin{2})
            for iter=1:maxiter
                uperf(I,J) = (p*uperf(I-1,J)+p*uperf(I+1,J)+ ...
                              p*uperf(I,J-1)+p*uperf(I,J+1)+ ...
                              q*uperf(I-1,J-1)+q*uperf(I-1,J+1)+ ...
                              q*uperf(I+1,J-1)+q*uperf(I+1,J+1)- ...
                              h^2*(2*f(I,J)/3+(f(I+1,J)+f(I-1,J)+ ...
                              f(I,J+1)+f(I,J-1))/12))/c;
            end
        else
            error("The only options are 5, '9m', or '9pm.'")
        end
    end
elseif nargin == 7
   % UPERF(F,H,U,STRING,VALUE,STRING,VALUE)
    uperf = varargin{1};
    if ~strcmpi('double',class(uperf))
        error("Input argument must be of type 'double.'")
    end
    [m,n] = size(uperf);
    if m ~= n || (m == 1 && n == 1)
        error('The input u must be a square matrix.')
    end
    if strcmpi('Maxiter',varargin{2})
        maxiter = varargin{3};
        if ~strcmpi('Stencil',varargin{4})
            error("The only option is the string 'Stencil'.");
        end
        if varargin{5} == 5
            for iter = 1:maxiter
                uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                                   uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
            end
        elseif strcmpi('9m',varargin{5})
            for iter=1:maxiter
                uperf(I,J) = 0.05*(4*uperf(I-1,J)+4*uperf(I+1,J)+ ...
                                   4*uperf(I,J-1)+4*uperf(I,J+1)+ ...
                                     uperf(I-1,J-1)+uperf(I-1,J+1)+ ... 
                                     uperf(I+1,J-1)+uperf(I+1,J+1)- ...
                                     6*h^2*(2*f(I,J)/3+(f(I+1,J)+...
                                     f(I-1,J)+f(I,J+1)+f(I,J-1))/12));
            end            
        elseif strcmpi('9pm',varargin{5})
            for iter=1:maxiter
                uperf(I,J) = (p*uperf(I-1,J)+p*uperf(I+1,J)+ ...
                              p*uperf(I,J-1)+p*uperf(I,J+1)+ ...
                              q*uperf(I-1,J-1)+q*uperf(I-1,J+1)+ ...
                              q*uperf(I+1,J-1)+q*uperf(I+1,J+1)- ...
                              h^2*(2*f(I,J)/3+(f(I+1,J)+f(I-1,J)+ ...
                              f(I,J+1)+f(I,J-1))/12))/c;
            end 
        else
            error("The only options are 5, '9m', or '9pm.'")
        end
    elseif strcmpi('Stencil',varargin{2})
        maxiter = varargin{5};
        if ~strcmpi('MaxIter',varargin{4})
            error("The only second option is the string 'MaxIter'.");
        end
        if varargin{3} == 5
            for iter = 1:maxiter
                uperf(I,J) = 0.25*(uperf(I-1,J)+uperf(I+1,J)+ ...
                                   uperf(I,J-1)+uperf(I,J+1)-h^2*f(I,J));
            end
        elseif strcmpi('9m',varargin{3})
            for iter=1:maxiter
                uperf(I,J) = 0.05*(4*uperf(I-1,J)+4*uperf(I+1,J)+ ...
                                   4*uperf(I,J-1)+4*uperf(I,J+1)+ ...
                                     uperf(I-1,J-1)+uperf(I-1,J+1)+ ... 
                                     uperf(I+1,J-1)+uperf(I+1,J+1)- ...
                                     6*h^2*(2*f(I,J)/3+(f(I+1,J)+...
                                     f(I-1,J)+f(I,J+1)+f(I,J-1))/12));
            end
        elseif strcmpi('9pm',varargin{3})
            for iter=1:maxiter
                uperf(I,J) = (p*uperf(I-1,J)+p*uperf(I+1,J)+ ...
                              p*uperf(I,J-1)+p*uperf(I,J+1)+ ...
                              q*uperf(I-1,J-1)+q*uperf(I-1,J+1)+ ...
                              q*uperf(I+1,J-1)+q*uperf(I+1,J+1)- ...
                              h^2*(2*f(I,J)/3+(f(I+1,J)+f(I-1,J)+ ...
                              f(I,J+1)+f(I,J-1))/12))/c;
            end
        else
            error("The only options are 5, '9m', or '9pm.'")
        end
    end 
else
    error('Incorrect number of input arguments.')
end
end
