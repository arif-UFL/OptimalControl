clear all; close all; clc;
%% Problem2 Main function
%% Given Variable;
t0 = 0;
x0 = 0;
y0=0;
v0=0;
xf=2;
yf=2;
g=10;
% lambdavf = 0;
%Hf = -1;
lambdax0Guess=.1;
lambday0Guess=.1;
lambdav0Guess=.1;
tfGuess=1;

l0 = [lambdax0Guess; lambday0Guess; lambdav0Guess; tfGuess];
options = optimset('Display','Iter','Tolx',1e-8,'TolFun',1e-8);
[p,theta] = fsolve(@error2,l0,options,x0,y0,v0,xf,yf,t0,g);

%% Numerically evaluated values
lambdaxeval=p(1);
lambdayeval=p(2);
lambdaveval=p(3);
tfeval=p(4);

L0 = [x0;y0;v0;lambdaxeval;lambdayeval;lambdaveval];
tspannew = [t0,tfeval];
options = odeset('RelTol',1e-8);
[t,L]= ode113(@myode2,tspannew,L0,options,g);


%% Plotting values
%% position and velocity plots 
figure(1)
plot(t,L(:,1:3))
hold on
xline(t(end))
title('Brochistochrome problem position and velocity components')
xlabel('time(seconds)')
ylabel('X(t),Y(t) & V(t)')
legend('X(t)','Y(t)','V(t)','optimal Time (0.8165)');


% %% Control Plot
% figure(2)
% plot(t(:,end),theta)
% hold on
% xline(t(end))
% title('control vs time')
% xlabel('time(seconds)')
% ylabel('theta')
% legend('theta');