
%% Brochistochrone Bonus 2
clear all;close all;clc

% Given
x0 = 0;
y0 = 0;
v0 = 0;
t0 = 0;
xf = 2;
yf = 2;
n = 3;% No. of states
lambdaxguess = 0;
lambdayguess = 0;
lambdavguess = 0;
tfguess = 10;
g = 10;
k =2; % 2 stage shooting
tau = linspace(-1,+1,k+1);
p0guess = zeros(2*n,k-1);
z0g = [lambdaxguess;lambdayguess;lambdavguess;tfguess;];
z0 = [z0g;p0guess(:)]; % Initial variables for shooting
options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
z = fsolve(@BSError,z0,options,t0,x0,y0,v0,xf,yf,g,n,k,tau);
[E,t,p] = BSError(z,t0,x0,y0,v0,xf,yf,g,n,k,tau);

%%Plots 

figure(1)
plot(t,p(:,1),'r*-',t,p(:,2),'b*-',t,p(:,3),'b--')
title('Trajectory and velocity components plot in tau grid');
legend('X(tau)','Y(tau)','V(tau)');
grid on;

% Evaluating the control plot

for j = 1: length(t)
vtj = p(j,3);
lambdaxtj = p(j,4);
lambdaytj = p(j,5);
lambdavtj = p(j,6);
thetatj(j) = fsolve(@solvecontrol,0,options,lambdaxtj,lambdaytj,lambdavtj,vtj,g);
end

% control plot
figure(2)
plot(t,thetatj)
hold on
title('Control vs Tau')
xlabel('Tau')
ylabel('u(tau)')
legend('Control u(tau)');

%% Optimal trajectory plot
figure(3)
plot(p(:,1),p(:,2))
set(gca,'YDir','reverse')
hold on
title('Optimal trajectory')
xlabel('X(tau)')
ylabel('Y(tau)')
xlim([0 2.00])
ylim([0 2.00])
legend('optimal Trajectory');
%% Error function
function [E,t,p] = BSError(z0,t0,x0,y0,v0,xf,yf,g,n,k,tau)
% p0 = [x0;y0;v0;z0(1);z0(2);z0(3)];
tfguess = z0(4);
premint = z0(5:end);
premint = reshape(premint,2*n,k-1);
options = odeset('RelTol',1e-8);
options1 = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
E = [];
p = [];
t = [];
for i = 1:k
    if i == 1
        p0 = [x0;y0;v0;z0(1);z0(2);z0(3)];
    else
        p0 = premint(:,k-1);
    end

    tauspan = [tau(i),tau(i+1)];
    [tout,pout] = ode113(@BSODE,tauspan,p0,options,t0,tfguess,g);
    ptf = pout(end,:)';

    if i < k
        E = [E;ptf - premint(:,i)];
    end
    t = [t; tout];
    p = [p; pout];
end
xtf = ptf(1);
ytf = ptf(2);
vtf = ptf(3);
lambdaxtf = ptf(4);
lambdaytf = ptf(5);
lambdavtf = ptf(6);
theta0guess = 0;
thetaf = fsolve(@solvecontrol,theta0guess,options1,lambdaxtf,lambdaytf,lambdavtf,vtf,g);
htfout = lambdaxtf*vtf*sin(thetaf)+(lambdaytf*vtf+lambdavtf*g)*cos(thetaf);
E = [E;ptf(1)-xf;ptf(2)-yf;lambdavtf;htfout+1];
end

%% ODE Function
function pdot = BSODE(t,p,t0,tf,g)
thetaGuess = 0;
v = p(3);
lambdax = p(4);
lambday = p(5);
lambdav = p(6);
options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
theta = fsolve(@solvecontrol,thetaGuess,options,lambdax,lambday,lambdav,v,g);
xdot = v*sin(theta);
ydot = v*cos(theta);
vdot = g*cos(theta);
lambdaxdot = 0;
lambdaydot = 0;
lambdavdot = -lambdax*sin(theta)-lambday*cos(theta);
pdot = [xdot;ydot;vdot;lambdaxdot;lambdaydot;lambdavdot];
pdot = ((tf-t0)/2)*pdot;
end

%% Solve control function
function htheta = solvecontrol(theta,lambdax,lambday,lambdav,v,g)
htheta = lambdax*v*cos(theta) - (lambday*v + lambdav*g)*sin(theta);
end

