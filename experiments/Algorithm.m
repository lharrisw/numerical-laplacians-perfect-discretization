%% Algorithm Part 1

clear;
clc;
close all;

s = 25; % sample size
n = 6; % number of columns for the erros
hvec = zeros(s,1); % initialization for step-sizes
Nvec = hvec; % initialization for matrix indices
evec = zeros(s,6); % initialization for errors

for j = n:(s+n-1)
% initialization
maxiter = 1000; % number of iterations
L = 2*pi; % length of intervals
h = L/(j+1); % step-size
hvec(j-n+1) = h; % store step-sizes
Nvec(j-n+1) = j; % store matrix indices

% forming the grid
[x,y] = meshgrid(0:h:L,0:h:L);

% right-hand sides and solutions 
% system 1
f = 2*sin(x).*sin(y); % rhs
soln = -sin(x).*sin(y); % solution

% % system 2
% f = 16.*((x-L/2).^2+(y-L/2).^2-.5).*exp(-2*((x-L/2).^2+(y-L/2).^2)); % rhs
% soln = exp(-2*((x-L/2).^2+(y-L/2).^2)); % solution

% % system 3
% f = 2*y.*(2*pi-y).*(sin(2*x)+x.*cos(2*x))-2*x.*sin(x).^2;
% soln = x.*y.*(2*pi-y).*sin(x).^2;

% poisson solvers
ustd5 = myPoisson1(f,h,'maxiter',maxiter,'stencil',5); % STD 5-pt operator
u9m = myPoisson1(f,h,'maxiter',maxiter,'stencil','9m'); % M 9-pt operator
u9pm = myPoisson1(f,h,'maxiter',maxiter,'stencil','9pm'); % PM 9-pt operator
up1 = myPoisson3(f,h,maxiter); % SN 25-point operator
up2 = myPoisson2(f,h,maxiter); % PNN 25-pt operator
up3 = perfPoisson(f,h,maxiter); % PN 25-point operator

% maximum error
eustd5 = norm(soln-ustd5,'inf'); % STD 5-pt error
eu9m = norm(soln-u9m,'inf'); % M 9-pt error
eu9pm = norm(soln-u9pm,'inf'); % PM 9-pt error
eup1 = norm(soln-up1,'inf'); % SN 25-pt error
eup2 = norm(soln-up2,'inf'); % PNN 25-pt error
eup3 = norm(soln-up3,'inf'); % PN 25-pt error

evec(j-n+1,:) = [eustd5 eu9m eu9pm eup1 eup2 eup3]; % store errors
% disp(j)
end

% % Visualization
% error for plots
eustd5 = abs(soln-ustd5);
eu9m = abs(soln-u9m);
eu9pm = abs(soln-u9pm);
eup1 = abs(soln-up1);
eup2 = abs(soln-up2);
eup3 = abs(soln-up3);

% error plots
figure
subplot(2,3,1)
mesh(x,y,eustd5)
title('STD 5-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,2)
mesh(x,y,eu9m)
title('M 9-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,3)
mesh(x,y,eu9pm)
title('PM 9-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,4)
mesh(x,y,eup1)
title('SN 25-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,5)
mesh(x,y,eup2)
title('PNN 25-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,6)
mesh(x,y,eup3)
title('PN 25-pt Error','interpreter','latex')
set(gca,'fontsize',20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution plots
figure
subplot(2,3,1)
mesh(x,y,ustd5)
title('STD 5-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,2)
mesh(x,y,u9m)
title('M 9-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,3)
mesh(x,y,u9pm)
title('PM 9-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,4)
mesh(x,y,up1)
title('SN 25-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,5)
mesh(x,y,up2)
title('PNN 25-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,6)
mesh(x,y,up3)
title('PN 25-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

%% Algorithm Part 2

clear;
clc;
close all;

s = 25; % sample size
n = 6; % number of columns for the erros
hvec = zeros(s,1); % initialization for step-sizes
Nvec = hvec; % initialization for matrix indices
evec = zeros(s,6); % initialization for errors

for j = n:(s+n-1)
% initialization
maxiter = 1000; % number of iterations
L = 2*pi; % length of intervals
h = L/(j+1); % step-size
hvec(j-n+1) = h; % store step-sizes
Nvec(j-n+1) = j; % store matrix indices

% forming the grid
[x,y] = meshgrid(0:h:L,0:h:L);

% right-hand sides and solutions
% Region 1
R = tetragon(x,y); % square region

% % Region 2
% R = ellipse(x,y); % circular region

% % Region 3
% R = triangle(x,y); % triangular region

% System 4
f = 2*sin(x).*sin(y).*R; % rhs
soln = -sin(x).*sin(y).*R; % solution

% poisson solvers
ustd5 = myPoisson1(f,h,'maxiter',maxiter,'stencil',5); % STD 5-pt operator
u9m = myPoisson1(f,h,'maxiter',maxiter,'stencil','9m'); % M 9-pt operator
u9pm = myPoisson1(f,h,'maxiter',maxiter,'stencil','9pm'); % PM 9-pt operator
up1 = myPoisson3(f,h,maxiter); % SN 25-point operator
up2 = myPoisson2(f,h,maxiter); % PNN 25-pt operator
up3 = perfPoisson(f,h,maxiter); % PN 25-point operator

% maximum error
eustd5 = norm(soln-ustd5,'inf'); % STD 5-pt error
eu9m = norm(soln-u9m,'inf'); % M 9-pt error
eu9pm = norm(soln-u9pm,'inf'); % PM 9-pt error
eup1 = norm(soln-up1,'inf'); % SN 25-pt error
eup2 = norm(soln-up2,'inf'); % PNN 25-pt error
eup3 = norm(soln-up3,'inf'); % PN 25-pt error

evec(j-n+1,:) = [eustd5 eu9m eu9pm eup1 eup2 eup3]; % store errors
% disp(j)
end

% % Visualization
% errors for plots using absolute values
eustd5 = abs(soln-ustd5); 
eu9m = abs(soln-u9m);
eu9pm = abs(soln-u9pm);
eup1 = abs(soln-up1);
eup2 = abs(soln-up2);
eup3 = abs(soln-up3);

% error plots
figure
subplot(2,3,1)
mesh(x,y,eustd5)
title('STD 5-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,2)
mesh(x,y,eu9m)
title('M 9-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,3)
mesh(x,y,eu9pm)
title('PM 9-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,4)
mesh(x,y,eup1)
title('SN 25-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,5)
mesh(x,y,eup2)
title('PNN 25-pt Error','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,6)
mesh(x,y,eup3)
title('PN 25-pt Error','interpreter','latex')
set(gca,'fontsize',20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution plots
figure
subplot(2,3,1)
mesh(x,y,ustd5)
title('STD 5-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,2)
mesh(x,y,u9m)
title('M 9-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,3)
mesh(x,y,u9pm)
title('PM 9-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,4)
mesh(x,y,up1)
title('SN 25-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,5)
mesh(x,y,up2)
title('PNN 25-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

subplot(2,3,6)
mesh(x,y,up3)
title('PN 25-pt Solution','interpreter','latex')
set(gca,'fontsize',20)

%% Algorithm 1 & 2 Error and OrderAnalysis
clear;
clc;
close all;

% % load data
load('datasys1.mat')
% load('datasys2.mat')
% load('datasys3.mat')
% load('datasys4_1.mat')
% load('datasys4_2.mat')
% load('datasys4_3.mat')

% labels for the graphs
strvec = {'STD 5-pt Error'...
          'M 9-pt Error'...
          'PM 9-pt Error'...
          'SN 25-pt Error'...
          'PNN 25-pt Error'...
          'PN 25-pt Error'};

% plot of the error against the step-size
figure
for j = 1:6
    subplot(2,3,j)
    plot(hvec,evec(:,j,1),'linewidth',2)
    xlabel('$\Delta x$','interpreter','latex')
    ylabel('$\varepsilon$','interpreter','latex')
    title(sprintf('%s',strvec{j}),'interpreter','latex')
    set(gca,'xdir','reverse','fontsize',20)
end

% log-log plots for order determination

strveclog = {'STD 5-pt Order'...
             'M 9-pt Order'...
             'PM 9-pt Order'...
             'SN 25-pt Order'...
             'PNN 25-pt Order'...
             'PN 25-pt Order'};

% determination of the order
figure
% idx = 10:25; % datasys1
% idx = 1:25; % datasys2, and all datasys4 mat files
idx = 5:20; % datasys3
idxlen = length(idx);
logevec = zeros(idxlen,6);
for j = 1:6
    % logarithmic data
    logevec(idx,j) = log(evec(idx,j,1));
    logNvec = log(hvec(idx));
    % fitting the data
    p = polyfit(logNvec,logevec(idx,j),1); 
    y = polyval(p,logNvec);
    % plotting
    subplot(2,3,j)
    hold on
    plot(logNvec,logevec(idx,j),'k.','markersize',20);
    plot(logNvec,y,'b-','linewidth',2);
    xlabel('$\log(\Delta x)$','interpreter','latex')
    ylabel('$\log(\varepsilon)$','interpreter','latex')
    title(sprintf('%s: %1.3f',strveclog{j},round(p(1),3)),'interpreter','latex')
    set(gca,'fontsize',20)
end
