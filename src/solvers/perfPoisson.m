function uperf = perfPoisson(f,h,maxiter)
%PERFPOISSON Solves the Poisson equation on a square.
%   UPERF = PERFPOISSON(F,H,MAXITER) solves the Poisson equation on the
%   square with homogeneous Dirichlet boundary conditions. The argument F
%   is the right-hand side of the Poisson equation. The argument H is
%   the step size of the x-interval and MAXITER is the maximum number of
%   iterations. The algorithm uses Jacobi's method and utilizes the
%   parameterized normalized 25-point perfect Laplacian. This consists of  
%   six operators that vary depending on the distance of the central point 
%   from the boundary.
%
%   Example: Solve Poisson's equation on the square with homogeneous
%            boundary conditions with a right-hand side of 
%            f(x,y) = peaks(x-a,y-a), for 0 <= x,y <= 2*a
%            for a = 1 using the perfPoisson function with 150 iterations.
%
%       N = 25;
%       maxiter = 150;
%       a = 1;
%       h = a/(N + 1);
%       [X,Y] = meshgrid(0:h:a,0:h:a);
%       f = peaks(X-a,Y-a);
%       u = perfPoisson(f,h,maxiter);
%       mesh(X,Y,f)
%
%   See also MYPOISSON1, MYPOISSON2, MYPOISSON3.

%  Copyright 2020 Latimer De'Shone Harris-Ward.

N = length(f(1,:)) - 2;
uperf = zeros(N+2,N+2);

% constructing the rho matrices of the perfect laplacian
rho1 = zeros(5);
rho2top = rho1;
rho3top = rho1;
rho4tr = rho1;
rho5trr = rho1;
rho6tr = rho1;

% rho1
rho1(:,1) = -[.00162;-.00068;-.002;-.00068;.00162];
rho1(:,2) = -[-.00068;-.19024;-.61773;-.19024;-.00068];
rho1(:,3) = -[-.002;-.61773;1;-.61773;-.002];
rho1(:,4) = rho1(:,2);
rho1(:,5) = rho1(:,1);

% rho2 top (uses rho1_mat)
rho2top(:,1) = -[0;-.0023;-.002;-.00068;.00162];
rho2top(:,2) = -[0;-.18956;-.61773;-.19024;-.00068];
rho2top(:,3) = -[0;-.61573;1;-.61773;-.002];
rho2top(:,4) = rho2top(:,2);
rho2top(:,5) = rho2top(:,1);

% rho2 bottom (uses rho1_mat)
rho2bot = flip(rho2top);

% rho2 left (uses rho1_mat)
rho2l = rot90(rho2top);

% rho2 right (uses rho1_mat)
rho2r = rot90(rho2bot);

% rho3 top
rho3top(:,1) = -[0;0;-.00132;-.0023;.00162];
rho3top(:,2) = -[0;0;-.42749;-.18956;-.00068];
rho3top(:,3) = -[0;0;1;-.61573;-.002];
rho3top(:,4) = rho3top(:,2);
rho3top(:,5) = rho3top(:,1);

% rho3 bottom
rho3bot = flip(rho3top);

% rho3 left
rho3l = rot90(rho3top);

% rho3 right
rho3r = rot90(rho3bot);

% rho4 top right (uses rho1_mat)
rho4tr(:,1) = -[0;-.0023;-.002;-.00068;.00162];
rho4tr(:,2) = -[0;-.18956;-.61773;-.19024;-.00068];
rho4tr(:,3) = -[0;-.61573;1;-.61773;-.002];
rho4tr(:,4) = -[0;-.18726;-.61573;-.18956;-.0023];

% rho4 top left (uses rho1_mat)
rho4tl = flip(rho4tr,2);

% rho4 bottom left (uses rho1_mat)
rho4bl = flip(rho4tl);

% rho4 bottom right (uses rho1_mat)
rho4br = flip(rho4tr);

% rho5 top right row (uses rho6tr_mat)
rho5trr(:,1) = -[0;0;-.00132;-.0023;.00162];
rho5trr(:,2) = -[0;0;-.42749;-.18956;-.00068];
rho5trr(:,3) = -[0;0;1;-.61573;-.002];
rho5trr(:,4) = -[0;0;-.42617;-.18726;-.0023];

% rho5 top right col (uses rho6tr_mat)
rho5trc = rot90(flip(rho5trr));

% rho5 top left row (uses rho6tl_mat)
rho5tlr = flip(rho5trr,2); 

% rho5 top left col (uses rho6tl_mat)
rho5tlc = rot90(flip(rho5tlr),-1);

% rho5 bottom left col (uses rho6bl_mat)
rho5blc = flip(rho5tlc);

% rho5 bottom left row (uses rho6bl_mat)
rho5blr = flip(rho5tlr);

% rho5 bottom right row (uses rho6br_mat)
rho5brr = flip(rho5trr);

% rho5 bottom right col (uses rho6br_mat)
rho5brc = flip(rho5trc);

% rho6 top right
rho6tr(:,1) = -[0;0;-.00132;-.00230;.00162];
rho6tr(:,2) = -[0;0;-.42617;-.18726;-.00230];
rho6tr(:,3) = -[0;0;1;-.42617;-.00132];

% rho6 top left
rho6tl = flip(rho6tr,2);

% rho6 bottom left
rho6bl = flip(rho6tl);

% rho6 bottom right
rho6br = flip(rho6tr);           

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
        uperf(3,jj) = rho2top(:)'*m2(:)/r1;
        % left column
        m3 = [uperf(jj-2,1) uperf(jj-2,2) uperf(jj-2,3) uperf(jj-2,4) uperf(jj-2,5);
              uperf(jj-1,1) uperf(jj-1,2) uperf(jj-1,3) uperf(jj-1,4) uperf(jj-1,5);
              uperf(jj,1)   uperf(jj,2)   h^2*f(jj,3)   uperf(jj,4)   uperf(jj,5);
              uperf(jj+1,1) uperf(jj+1,2) uperf(jj+1,3) uperf(jj+1,4) uperf(jj+1,5);
              uperf(jj+2,1) uperf(jj+2,2) uperf(jj+2,3) uperf(jj+2,4) uperf(jj+2,5)];
        uperf(jj,3) = rho2l(:)'*m3(:)/r1;
        % bottom row
        m4 = [uperf(N-2,jj-2) uperf(N-2,jj-1) uperf(N-2,jj) uperf(N-2,jj+1) uperf(N-2,jj+2);
              uperf(N-1,jj-2) uperf(N-1,jj-1) uperf(N-1,jj) uperf(N-1,jj+1) uperf(N-1,jj+2);
              uperf(N,jj-2)   uperf(N,jj-1)   h^2*f(N,jj)   uperf(N,jj+1)   uperf(N,jj+2);
              uperf(N+1,jj-2) uperf(N+1,jj-1) uperf(N+1,jj) uperf(N+1,jj+1) uperf(N+1,jj+2);
              uperf(N+2,jj-2) uperf(N+2,jj-1) uperf(N+2,jj) uperf(N+2,jj+1) uperf(N+2,jj+2)];
        uperf(N,jj) = rho2bot(:)'*m4(:)/r1;
        % right column
        m5 = [uperf(jj-2,N-2) uperf(jj-2,N-1) uperf(jj-2,N) uperf(jj-2,N+1) uperf(jj-2,N+2);
              uperf(jj-1,N-2) uperf(jj-1,N-1) uperf(jj-1,N) uperf(jj-1,N+1) uperf(jj-1,N+2);
              uperf(jj,N-2)   uperf(jj,N-1)   h^2*f(jj,N)   uperf(jj,N+1)   uperf(jj,N+2);
              uperf(jj+1,N-2) uperf(jj+1,N-1) uperf(jj+1,N) uperf(jj+1,N+1) uperf(jj+1,N+2);
              uperf(jj+2,N-2) uperf(jj+2,N-1) uperf(jj+2,N) uperf(jj+2,N+1) uperf(jj+2,N+2)];
        uperf(jj,N) = rho2r(:)'*m5(:)/r1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho3 laplacian
    for jj = 3:N
        r3 = 3.85654;
        % top row
        m6 = [0 0 0 0 0;
              0 0 0 0 0;
              uperf(2,jj-2) uperf(2,jj-1) h^2*f(2,jj) uperf(2,jj+1) uperf(2,jj+2);
              uperf(3,jj-2) uperf(3,jj-1) uperf(3,jj) uperf(3,jj+1) uperf(3,jj+2);
              uperf(4,jj-2) uperf(4,jj-1) uperf(4,jj) uperf(4,jj+1) uperf(4,jj+2)];
        uperf(2,jj) = rho3top(:)'*m6(:)/r3;
        % left column
        m7 = [0 0 uperf(jj-2,2) uperf(jj-2,3) uperf(jj-2,4);
              0 0 uperf(jj-1,2) uperf(jj-1,3) uperf(jj-1,4);
              0 0 h^2*f(jj,2)   uperf(jj,3)   uperf(jj,4);
              0 0 uperf(jj+1,2) uperf(jj+1,3) uperf(jj+1,4);
              0 0 uperf(jj+2,2) uperf(jj+2,3) uperf(jj+2,4)];
        uperf(jj,2) = rho3l(:)'*m7(:)/r3;
        % bottom row 
        m8 = [uperf(N-1,jj-2) uperf(N-1,jj-1) uperf(N-1,jj) uperf(N-1,jj+1) uperf(N-1,jj+2);
              uperf(N,jj-2)   uperf(N,jj-1)   uperf(N,jj)   uperf(N,jj+1)   uperf(N,jj+2);
              uperf(N+1,jj-2) uperf(N+1,jj-1) h^2*f(N+1,jj) uperf(N+1,jj+1) uperf(N+1,jj+2);
              0 0 0 0 0;
              0 0 0 0 0];
        uperf(N+1,jj) = rho3bot(:)'*m8(:)/r3;
        % right column
        m9 = [uperf(jj-2,N-1) uperf(jj-2,N) uperf(jj-2,N+1) 0 0;
              uperf(jj-1,N-1) uperf(jj-1,N) uperf(jj-1,N+1) 0 0;
              uperf(jj,N-1)   uperf(jj,N)   h^2*f(jj,N+1)   0 0;
              uperf(jj+1,N-1) uperf(jj+1,N) uperf(jj+1,N+1) 0 0;
              uperf(jj+2,N-1) uperf(jj+2,N) uperf(jj+2,N+1) 0 0];
        uperf(jj,N+1) = rho3r(:)'*m9(:)/r3;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho4 laplacian
    % top left sub-corner
    m10 = [uperf(1,1) uperf(1,2) uperf(1,3) uperf(1,4) uperf(1,5);
           uperf(2,1) uperf(2,2) uperf(2,3) uperf(2,4) uperf(2,5);
           uperf(3,1) uperf(3,2) h^2*f(3,3) uperf(3,4) uperf(3,5);
           uperf(4,1) uperf(4,2) uperf(4,3) uperf(4,4) uperf(4,5);
           uperf(5,1) uperf(5,2) uperf(5,3) uperf(5,4) uperf(5,5)];
    uperf(3,3) = rho4tl(:)'*m10(:)/r1;
    % bottom left sub-corner
    m11 = [uperf(N-2,1) uperf(N-2,2) uperf(N-2,3) uperf(N-2,4) uperf(N-2,5);
           uperf(N-1,1) uperf(N-1,2) uperf(N-1,3) uperf(N-1,4) uperf(N-1,5);
           uperf(N,1)   uperf(N,2)   h^2*f(N,3)   uperf(N,4)   uperf(N,5);
           uperf(N+1,1) uperf(N+1,2) uperf(N+1,3) uperf(N+1,4) uperf(N+1,5);
           uperf(N+2,1) uperf(N+2,2) uperf(N+2,3) uperf(N+2,4) uperf(N+2,5)];
    uperf(N,3) = rho4bl(:)'*m11(:)/r1;
    % top right sub-corner
    m12 = [uperf(1,N-2) uperf(1,N-1) uperf(1,N) uperf(1,N+1) uperf(1,N+2);
           uperf(2,N-2) uperf(2,N-1) uperf(2,N) uperf(2,N+1) uperf(2,N+2);
           uperf(3,N-2) uperf(3,N-1)  h^2*f(3,N) uperf(3,N+1) uperf(3,N+2);
           uperf(4,N-2) uperf(4,N-1) uperf(4,N) uperf(4,N+1) uperf(4,N+2);
           uperf(5,N-2) uperf(5,N-1) uperf(5,N) uperf(5,N+1) uperf(5,N+2)];
    uperf(3,N) = rho4tr(:)'*m12(:)/r1;
    % bottom right sub-corner
    m13 = [uperf(N-2,N-2) uperf(N-2,N-1) uperf(N-2,N) uperf(N-2,N+1) uperf(N-2,N+2);
           uperf(N-1,N-2) uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) uperf(N-1,N+2);
           uperf(N,N-2)   uperf(N,N-1)   h^2*f(N,N)   uperf(N,N+1)   uperf(N,N+2);
           uperf(N+1,N-2) uperf(N+1,N-1) uperf(N+1,N) uperf(N+1,N+1) uperf(N+1,N+2);
           uperf(N+2,N-2) uperf(N+2,N-1) uperf(N+2,N) uperf(N+2,N+1) uperf(N+2,N+2)];
    uperf(N,N) = rho4br(:)'*m13(:)/r1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho5 laplacian
    % top row %
    % top right row
    m14 = [0 0 0 0 0;
           uperf(1,N-2) uperf(1,N-1) uperf(1,N) uperf(1,N+1) 0;
           uperf(2,N-2) uperf(2,N-1) h^2*f(2,N) uperf(2,N+1) 0;
           uperf(3,N-2) uperf(3,N-1) uperf(3,N) uperf(3,N+1) 0;
           uperf(4,N-2) uperf(4,N-1) uperf(4,N) uperf(4,N+1) 0];
    uperf(2,N) = rho5trr(:)'*m14(:)/r3;
    % top left row
    m15 = [0 0 0 0 0;
           0 uperf(1,2) uperf(1,3) uperf(1,4) uperf(1,5);
           0 uperf(2,2) h^2*f(2,3) uperf(2,4) uperf(2,5);
           0 uperf(3,2) uperf(3,3) uperf(3,4) uperf(3,5);
           0 uperf(4,2) uperf(4,3) uperf(4,4) uperf(4,5)];
    uperf(2,3) = rho5tlr(:)'*m15(:)/r3;
    % bottom row %
    % bottom right row
    m16 = [uperf(N-1,N-2) uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) 0;
           uperf(N,N-2)   uperf(N,N-1)   uperf(N,N)   uperf(N,N+1) 0;
           uperf(N+1,N-2) uperf(N+1,N-1) h^2*f(N+1,N) uperf(N+1,N+1) 0;
           uperf(N+2,N-2) uperf(N+2,N-1) uperf(N+2,N) uperf(N+2,N+1) 0;
           0 0 0 0 0];
    uperf(N+1,N) = rho5brr(:)'*m16(:)/r3;
    % bottom left row
    m17 = [0 uperf(N-1,2) uperf(N-1,3) uperf(N-1,4) uperf(N-1,5);
           0 uperf(N,2)   uperf(N,3)   uperf(N,4)   uperf(N,5);
           0 uperf(N+1,2) h^2*f(N+1,3) uperf(N+1,4) uperf(N+1,5);
           0 uperf(N+2,2) uperf(N+2,3) uperf(N+2,4) uperf(N+2,5);
           0 0 0 0 0];
    uperf(N+1,3) = rho5blr(:)'*m17(:)/r3;
    % left column %
    % top left column
    m18 = [0 0 0 0 0;
           0 uperf(2,1) uperf(2,2) uperf(2,3) uperf(2,4);
           0 uperf(3,1) h^2*f(3,2) uperf(3,3) uperf(3,4);
           0 uperf(4,1) uperf(4,2) uperf(4,3) uperf(4,4);
           0 uperf(5,1) uperf(5,2) uperf(5,3) uperf(5,4)];
    uperf(3,2) = rho5tlc(:)'*m18(:)/r3;
    % bottom left column
    m19 = [0 uperf(N-2,1) uperf(N-2,2) uperf(N-2,3) uperf(N-2,4);
           0 uperf(N-1,1) uperf(N-1,2) uperf(N-1,3) uperf(N-1,4);
           0 uperf(N,1)   h^2*f(N,2)   uperf(N,3)   uperf(N,4);
           0 uperf(N+1,1) uperf(N+1,2) uperf(N+1,3) uperf(N+1,4);
           0 0 0 0 0];
    uperf(N,2) = rho5blc(:)'*m19(:)/r3;
    % right column %
    % top right column
    m20 = [0 0 0 0 0;
           uperf(2,N-1) uperf(2,N) uperf(2,N+1) uperf(2,N+2) 0;
           uperf(3,N-1) uperf(3,N) h^2*f(3,N+1) uperf(3,N+2) 0;
           uperf(4,N-1) uperf(4,N) uperf(4,N+1) uperf(4,N+2) 0;
           uperf(5,N-1) uperf(5,N) uperf(5,N+1) uperf(5,N+2) 0];
    uperf(3,N+1) = rho5trc(:)'*m20(:)/r3;
    % bottom right column
    m21 = [uperf(N-2,N-1) uperf(N-2,N) uperf(N-2,N+1) uperf(N-2,N+2) 0;
           uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) uperf(N-1,N+2) 0;
           uperf(N,N-1)   uperf(N,N)   h^2*f(N,N+1)   uperf(N,N+2)   0;
           uperf(N+1,N-1) uperf(N+1,N) uperf(N+1,N+1) uperf(N+1,N+2) 0;
           0 0 0 0 0];
    uperf(N,N+1) = rho5brc(:)'*m21(:)/r3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho 6 laplacian
    r6 = 4.28403;
    % top right corner
    m22 = [0 0 0 0 0;
           uperf(1,N-1) uperf(1,N) uperf(1,N+1) uperf(1,N+2) 0;
           uperf(2,N-1) uperf(2,N) h^2*f(2,N+1) uperf(2,N+2)  0;
           uperf(3,N-1) uperf(3,N) uperf(3,N+1) uperf(3,N+2) 0;
           uperf(4,N-1) uperf(4,N) uperf(4,N+1) uperf(4,N+2) 0];
    uperf(2,N+1) = rho6tr(:)'*m22(:)/r6;
    % top left corner
    m23 = [0 0 0 0 0;
           0 uperf(1,1) uperf(1,2) uperf(1,3) uperf(1,4);
           0 uperf(2,1) h^2*f(2,2) uperf(2,3) uperf(2,4);
           0 uperf(3,1) uperf(3,2) uperf(3,3) uperf(3,4);
           0 uperf(4,1) uperf(4,2) uperf(4,3) uperf(4,4)];
    uperf(2,2) = rho6tl(:)'*m23(:)/r6;
    % bottom left corner
    m24 = [0 uperf(N-1,1) uperf(N-1,2) uperf(N-1,3) uperf(N-1,4);
           0 uperf(N,1)   uperf(N,2)   uperf(N,3)   uperf(N,4);
           0 uperf(N+1,1) h^2*f(N+1,2) uperf(N+1,3) uperf(N+1,4);
           0 uperf(N+2,1) uperf(N+2,2) uperf(N+2,3) uperf(N+2,4);
           0 0 0 0 0];
    uperf(N+1,2) = rho6bl(:)'*m24(:)/r6;
    % bottom right corner
    m25 = [uperf(N-1,N-1) uperf(N-1,N) uperf(N-1,N+1) uperf(N-1,N+2) 0;
           uperf(N,N-1)   uperf(N,N)   uperf(N,N+1)   uperf(N,N+2)   0;
           uperf(N+1,N-1) uperf(N+1,N) h^2*f(N+1,N+1) uperf(N+1,N+2) 0;
           uperf(N+2,N-1) uperf(N+2,N) uperf(N+2,N+1) uperf(N+2,N+2) 0;
           0 0 0 0 0];
    uperf(N+1,N+1) = rho6br(:)'*m25(:)/r6;
end
end
