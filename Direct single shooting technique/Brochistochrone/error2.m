function [E,theta] = error2(l0,x0,y0,v0,xf,yf,t0,g)
% lambdax0Guess = l0(1,1);
% lambday0Guess = l0(2,1);
% lambdav0Guess = l0(3,1);
tf = l0(end);
% x0 = 0;
% y0 = 0;
% v0 = 0;
% xf = 2;
% yf = 2;
% g=10;

l0 = [x0;y0;v0;l0(1:3)];
options = odeset('RelTol',1e-8);
tspan = [t0,tf];
[t,p] = ode113(@myode2,tspan,l0,options,g);
peval=p(end,:);
lambdaxtf = peval(1,4);
lambdaytf = peval(1,5);
lambdavtf = peval(1,6);
xtf = peval(1,1);
ytf = peval(1,2);
vtf = peval(1,3);
thetaguess = 0;
options = optimset('Display','off','Tolx',1e-8,'TolFun',1e-8);
theta = fsolve(@SolveControl,thetaguess,options,lambdaxtf,lambdaytf,lambdavtf,vtf,g);
Htf = lambdaxtf*vtf*sin(theta) + (lambdaytf*vtf + lambdavtf*g)*cos(theta);
%% lambdavf = 0 implies p(3,1) = 0;E = P(3,1)
E = [xtf - xf; ytf - yf;lambdavtf; Htf+1];

end