function uperf = myPoisson3(f,h,maxiter)
%MYPOISSON3 Solves the Poisson equation on a square.
%   UPERF = MYPOISSON3(F,H,MAXITER) solves the Poisson equation on the
%   square with homogeneous Dirichlet boundary conditions. The argument F
%   is the right-hand side of the Poisson equation. The argument H is
%   the step size of the x-interval and MAXITER is the maximum number of
%   iterations. The algorithm uses Jacobi's method and utilizes the
%   single normalized 25-point perfect Laplacian.
%
%   Example: Solve Poisson's equation on the square with homogeneous
%            boundary conditions with a right-hand side of 
%            f(x,y) = exp(-((x-a)^2+(x-a)^2)), 0 <= x,y <= 2*a,
%            for a = 1 using the 9-point stencil with 150 iterations.
%
%       N = 25;
%       maxiter = 150;
%       a = 1;
%       h = a/(N + 1);
%       [X,Y] = meshgrid(0:h:2*a,0:h:2*a);
%       f = exp(-((x-a).^2+(x-a).^2));
%       u = myPoisson3(f,h,maxiter);
%       mesh(X,Y,f)
%
%   See also MYPOISSON1, MYPOISSON2, PERFPOISSON.

%  Copyright 2020 Latimer De'Shone Harris-Ward.

N = length(f(1,:)) - 2;
uperf = zeros(N+2,N+2);

% rho1
rho1 = -[.00162 -.00068 -.00200 -.00068 .00162;
        -.00068 -.19024 -.61773 -.19024 -.00068;
        -.00200 -.61773 1 -.61773 -.00200;
        -.00068 -.19024 -.61773 -.19024 -.00068;
         .00162 -.00068 -.00200 -.00068 .00162];

% perfect Laplacian
for iter = 1:maxiter
    r1 = 3.23881;
    % rho1 laplacian
    for ii=4:N-1
        for jj=4:N-1
            m1 = [uperf(ii-2,jj-2) uperf(ii-2,jj-1) uperf(ii-2,jj) uperf(ii-2,jj+1) uperf(ii-2,jj+2);
                  uperf(ii-1,jj-2) uperf(ii-1,jj-1) uperf(ii-1,jj) uperf(ii-1,jj+1) uperf(ii-1,jj+2);
                  uperf(ii,jj-2)   uperf(ii,jj-1)   h^2*f(ii,jj)   uperf(ii,jj+1)   uperf(ii,jj+2);
                  uperf(ii+1,jj-2) uperf(ii+1,jj-1) uperf(ii+1,jj) uperf(ii+1,jj+1) uperf(ii+1,jj+2);
                  uperf(ii+2,jj-2) uperf(ii+2,jj-1) uperf(ii+2,jj) uperf(ii+2,jj+1) uperf(ii+2,jj+2)];
            uperf(ii,jj) = rho1(:)'*m1(:)/r1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho2 laplacian
    for jj=3:N
        % top row
        m2 = [uperf(1,jj-2) uperf(1,jj-1) uperf(1,jj) uperf(1,jj+1) uperf(1,jj+2);
              uperf(2,jj-2) uperf(2,jj-1) uperf(2,jj) uperf(2,jj+1) uperf(2,jj+2);
              uperf(3,jj-2) uperf(3,jj-1) h^2*f(3,jj) uperf(3,jj+1) uperf(3,jj+2);
              uperf(4,jj-2) uperf(4,jj-1) uperf(4,jj) uperf(4,jj+1) uperf(4,jj+2);
              uperf(5,jj-2) uperf(5,jj-1) uperf(5,jj) uperf(5,jj+1) uperf(5,jj+2)];
        uperf(3,jj) = rho1(:)'*m2(:)/r1;
        % left column
        m3 = [uperf(jj-2,1) uperf(jj-2,2) uperf(jj-2,3) uperf(jj-2,4) uperf(jj-2,5);
              uperf(jj-1,1) uperf(jj-1,2) uperf(jj-1,3) uperf(jj-1,4) uperf(jj-1,5);
              uperf(jj,1)   uperf(jj,2)   h^2*f(jj,3)   uperf(jj,4)   uperf(jj,5);
              uperf(jj+1,1) uperf(jj+1,2) uperf(jj+1,3) uperf(jj+1,4) uperf(jj+1,5);
              uperf(jj+2,1) uperf(jj+2,2) uperf(jj+2,3) uperf(jj+2,4) uperf(jj+2,5)];
        uperf(jj,3) = rho1(:)'*m3(:)/r1;
        % bottom row
        m4 = [uperf(N-2,jj-2) uperf(N-2,jj-1) uperf(N-2,jj) uperf(N-2,jj+1) uperf(N-2,jj+2);
              uperf(N-1,jj-2) uperf(N-1,jj-1) uperf(N-1,jj) uperf(N-1,jj+1) uperf(N-1,jj+2);
              uperf(N,jj-2)   uperf(N,jj-1)   h^2*f(N,jj)   uperf(N,jj+1)   uperf(N,jj+2);
              uperf(N+1,jj-2) uperf(N+1,jj-1) uperf(N+1,jj) uperf(N+1,jj+1) uperf(N+1,jj+2);
              uperf(N+2,jj-2) uperf(N+2,jj-1) uperf(N+2,jj) uperf(N+2,jj+1) uperf(N+2,jj+2)];
        uperf(N,jj) = rho1(:)'*m4(:)/r1;
        % right column
        m5 = [uperf(jj-2,N-2) uperf(jj-2,N-1) uperf(jj-2,N) uperf(jj-2,N+1) uperf(jj-2,N+2);
              uperf(jj-1,N-2) uperf(jj-1,N-1) uperf(jj-1,N) uperf(jj-1,N+1) uperf(jj-1,N+2);
              uperf(jj,N-2)   uperf(jj,N-1)   h^2*f(jj,N)   uperf(jj,N+1)   uperf(jj,N+2);
              uperf(jj+1,N-2) uperf(jj+1,N-1) uperf(jj+1,N) uperf(jj+1,N+1) uperf(jj+1,N+2);
              uperf(jj+2,N-2) uperf(jj+2,N-1) uperf(jj+2,N) uperf(jj+2,N+1) uperf(jj+2,N+2)];
        uperf(jj,N) = rho1(:)'*m5(:)/r1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho3 laplacian
    for jj = 3:N
        % top row
        m6 = [0 0 0 0 0;
              0 0 0 0 0;
              uperf(2,jj-2) uperf(2,jj-1) h^2*f(2,jj) uperf(2,jj+1) uperf(2,jj+2);
              uperf(3,jj-2) uperf(3,jj-1) uperf(3,jj) uperf(3,jj+1) uperf(3,jj+2);
              uperf(4,jj-2) uperf(4,jj-1) uperf(4,jj) uperf(4,jj+1) uperf(4,jj+2)];
        uperf(2,jj) = rho1(:)'*m6(:)/r1;
        % left column
        m7 = [0 0 uperf(jj-2,2) uperf(jj-2,3) uperf(jj-2,4);
              0 0 uperf(jj-1,2) uperf(jj-1,3) uperf(jj-1,4);
              0 0 h^2*f(jj,2)   uperf(jj,3)   uperf(jj,4);
              0 0 uperf(jj+1,2) uperf(jj+1,3) uperf(jj+1,4);
              0 0 uperf(jj+2,2) uperf(jj+2,3) uperf(jj+2,4)];
        uperf(jj,2) = rho1(:)'*m7(:)/r1;
        % bottom row 
        m8 = [uperf(N-1,jj-2) uperf(N-1,jj-1) uperf(N-1,jj) uperf(N-1,jj+1) uperf(N-1,jj+2);
              uperf(N,jj-2)   uperf(N,jj-1)   uperf(N,jj)   uperf(N,jj+1)   uperf(N,jj+2);
              uperf(N+1,jj-2) uperf(N+1,jj-1) h^2*f(N+1,jj) uperf(N+1,jj+1) uperf(N+1,jj+2);
              0 0 0 0 0;
              0 0 0 0 0];
        uperf(N+1,jj) = rho1(:)'*m8(:)/r1;
        % right column
        m9 = [uperf(jj-2,N-1) uperf(jj-2,N) uperf(jj-2,N+1) 0 0;
              uperf(jj-1,N-1) uperf(jj-1,N) uperf(jj-1,N+1) 0 0;
              uperf(jj,N-1)   uperf(jj,N)   h^2*f(jj,N+1)   0 0;
              uperf(jj+1,N-1) uperf(jj+1,N) uperf(jj+1,N+1) 0 0;
              uperf(jj+2,N-1) uperf(jj+2,N) uperf(jj+2,N+1) 0 0];
        uperf(jj,N+1) = rho1(:)'*m9(:)/r1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho4 laplacian;
    % top left sub-corner
    m10 = [uperf(1,1) uperf(1,2) uperf(1,3) uperf(1,4) uperf(1,5);
           uperf(2,1) uperf(2,2) uperf(2,3) uperf(2,4) uperf(2,5);
           uperf(3,1) uperf(3,2) h^2*f(3,3) uperf(3,4) uperf(3,5);
           uperf(4,1) uperf(4,2) uperf(4,3) uperf(4,4) uperf(4,5);
           uperf(5,1) uperf(5,2) uperf(5,3) uperf(5,4) uperf(5,5)];
    uperf(3,3) = rho1(:)'*m10(:)/r1;
    % bottom left sub-corner
    m11 = [uperf(N-2,1) uperf(N-2,2) uperf(N-2,3) uperf(N-2,4) uperf(N-2,5);
           uperf(N-1,1) uperf(N-1,2) uperf(N-1,3) uperf(N-1,4) uperf(N-1,5);
           uperf(N,1)   uperf(N,2)   h^2*f(N,3)   uperf(N,4)   uperf(N,5);
           uperf(N+1,1) uperf(N+1,2) uperf(N+1,3) uperf(N+1,4) uperf(N+1,5);
           uperf(N+2,1) uperf(N+2,2) uperf(N+2,3) uperf(N+2,4) uperf(N+2,5)];
    uperf(N,3) = rho1(:)'*m11(:)/r1;
    % top right sub-corner
    m12 = [uperf(1,N-2) uperf(1,N-1) uperf(1,N) uperf(1,N+1) uperf(1,N+2);
           uperf(2,N-2) uperf(2,N-1) uperf(2,N) uperf(2,N+1) uperf(2,N+2);
           uperf(3,N-2) uperf(3,N-1) h^2*f(3,N) uperf(3,N+1) uperf(3,N+2);
           uperf(4,N-2) uperf(4,N-1) uperf(4,N) uperf(4,N+1) uperf(4,N+2);
           uperf(5,N-2) uperf(5,N-1) uperf(5,N) uperf(5,N+1) uperf(5,N+2)];
    uperf(3,N) = rho1(:)'*m12(:)/r1;
    % bottom right sub-corner
    m13 = [uperf(N-2,N-2) uperf(N-2,N-1) uperf(N-2,N) uperf(N-2,N+1) uperf(N-2,N+2);
           uperf(N-1,N-2) uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) uperf(N-1,N+2);
           uperf(N,N-2)   uperf(N,N-1)   h^2*f(N,N)   uperf(N,N+1)   uperf(N,N+2);
           uperf(N+1,N-2) uperf(N+1,N-1) uperf(N+1,N) uperf(N+1,N+1) uperf(N+1,N+2);
           uperf(N+2,N-2) uperf(N+2,N-1) uperf(N+2,N) uperf(N+2,N+1) uperf(N+2,N+2)];
    uperf(N,N) = rho1(:)'*m13(:)/r1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho5 laplacian
    % top row %
    % top right row
    m14 = [0 0 0 0 0;
           uperf(1,N-2) uperf(1,N-1) uperf(1,N) uperf(1,N+1) 0;
           uperf(2,N-2) uperf(2,N-1) h^2*f(2,N) uperf(2,N+1) 0;
           uperf(3,N-2) uperf(3,N-1) uperf(3,N) uperf(3,N+1) 0;
           uperf(4,N-2) uperf(4,N-1) uperf(4,N) uperf(4,N+1) 0];
    uperf(2,N) = rho1(:)'*m14(:)/r1;
    % top left row
    m15 = [0 0 0 0 0;
           0 uperf(1,2) uperf(1,3) uperf(1,4) uperf(1,5);
           0 uperf(2,2) h^2*f(2,3) uperf(2,4) uperf(2,5);
           0 uperf(3,2) uperf(3,3) uperf(3,4) uperf(3,5);
           0 uperf(4,2) uperf(4,3) uperf(4,4) uperf(4,5)];
    uperf(2,3) = rho1(:)'*m15(:)/r1;
    % bottom row %
    % bottom right row
    m16 = [uperf(N-1,N-2) uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) 0;
           uperf(N,N-2)   uperf(N,N-1)   uperf(N,N)   uperf(N,N+1) 0;
           uperf(N+1,N-2) uperf(N+1,N-1) h^2*f(N+1,N) uperf(N+1,N+1) 0;
           uperf(N+2,N-2) uperf(N+2,N-1) uperf(N+2,N) uperf(N+2,N+1) 0;
           0 0 0 0 0];
    uperf(N+1,N) = rho1(:)'*m16(:)/r1;
    % bottom left row
    m17 = [0 uperf(N-1,2) uperf(N-1,3) uperf(N-1,4) uperf(N-1,5);
           0 uperf(N,2)   uperf(N,3)   uperf(N,4)   uperf(N,5);
           0 uperf(N+1,2) h^2*f(N+1,3) uperf(N+1,4) uperf(N+1,5);
           0 uperf(N+2,2) uperf(N+2,3) uperf(N+2,4) uperf(N+2,5);
           0 0 0 0 0];
    uperf(N+1,3) = rho1(:)'*m17(:)/r1;
    % left column %
    % top left column
    m18 = [0 0 0 0 0;
           0 uperf(2,1) uperf(2,2) uperf(2,3) uperf(2,4);
           0 uperf(3,1) h^2*f(3,2) uperf(3,3) uperf(3,4);
           0 uperf(4,1) uperf(4,2) uperf(4,3) uperf(4,4);
           0 uperf(5,1) uperf(5,2) uperf(5,3) uperf(5,4)];
    uperf(3,2) = rho1(:)'*m18(:)/r1;
    % bottom left column
    m19 = [0 uperf(N-2,1) uperf(N-2,2) uperf(N-2,3) uperf(N-2,4);
           0 uperf(N-1,1) uperf(N-1,2) uperf(N-1,3) uperf(N-1,4);
           0 uperf(N,1)   h^2*f(N,2)   uperf(N,3)   uperf(N,4);
           0 uperf(N+1,1) uperf(N+1,2) uperf(N+1,3) uperf(N+1,4);
           0 0 0 0 0];
    uperf(N,2) = rho1(:)'*m19(:)/r1;
    % right column %
    % top right column
    m20 = [0 0 0 0 0;
           uperf(2,N-1) uperf(2,N) uperf(2,N+1) uperf(2,N+2) 0;
           uperf(3,N-1) uperf(3,N) h^2*f(3,N+1) uperf(3,N+2) 0;
           uperf(4,N-1) uperf(4,N) uperf(4,N+1) uperf(4,N+2) 0;
           uperf(5,N-1) uperf(5,N) uperf(5,N+1) uperf(5,N+2) 0];
    uperf(3,N+1) = rho1(:)'*m20(:)/r1;
    % bottom right column
    m21 = [uperf(N-2,N-1) uperf(N-2,N) uperf(N-2,N+1) uperf(N-2,N+2) 0;
           uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) uperf(N-1,N+2) 0;
           uperf(N,N-1)   uperf(N,N)   h^2*f(N,N+1)   uperf(N,N+2)   0;
           uperf(N+1,N-1) uperf(N+1,N) uperf(N+1,N+1) uperf(N+1,N+2) 0;
           0 0 0 0 0];
    uperf(N,N+1) = rho1(:)'*m21(:)/r1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho 6 laplacian
    % top right corner
    m22 = [0 0 0 0 0;
           uperf(1,N-1) uperf(1,N) uperf(1,N+1) uperf(1,N+2) 0;
           uperf(2,N-1) uperf(2,N) h^2*f(2,N+1) uperf(2,N+2)  0;
           uperf(3,N-1) uperf(3,N) uperf(3,N+1) uperf(3,N+2) 0;
           uperf(4,N-1) uperf(4,N) uperf(4,N+1) uperf(4,N+2) 0];
    uperf(2,N+1) = rho1(:)'*m22(:)/r1;
    % top left corner
    m23 = [0 0 0 0 0;
           0 uperf(1,1) uperf(1,2) uperf(1,3) uperf(1,4);
           0 uperf(2,1) h^2*f(2,2) uperf(2,3) uperf(2,4);
           0 uperf(3,1) uperf(3,2) uperf(3,3) uperf(3,4);
           0 uperf(4,1) uperf(4,2) uperf(4,3) uperf(4,4)];
    uperf(2,2) = rho1(:)'*m23(:)/r1;
    % bottom left corner
    m24 = [0 uperf(N-1,1) uperf(N-1,2) uperf(N-1,3) uperf(N-1,4);
           0 uperf(N,1)   uperf(N,2)   uperf(N,3)   uperf(N,4);
           0 uperf(N+1,1) h^2*f(N+1,2) uperf(N+1,3) uperf(N+1,4);
           0 uperf(N+2,1) uperf(N+2,2) uperf(N+2,3) uperf(N+2,4);
           0 0 0 0 0];
    uperf(N+1,2) = rho1(:)'*m24(:)/r1;
    % bottom right corner
    m25 = [uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) uperf(N-1,N+2) 0;
           uperf(N,N-1)   uperf(N,N)   uperf(N,N+1)   uperf(N,N+2)   0;
           uperf(N+1,N-1) uperf(N+1,N) h^2*f(N+1,N+1) uperf(N+1,N+2) 0;
           uperf(N+2,N-1) uperf(N+2,N) uperf(N+2,N+1) uperf(N+2,N+2) 0;
           0 0 0 0 0];
    uperf(N+1,N+1) = rho1(:)'*m25(:)/r1;
end
end
