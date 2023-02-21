%% Orbital transfer problem
%-----------------------------------------------------------------%
% Minimize t_f
% subjected to differential equation constraints
%  r_dot = vr;
% th_dot = vth / r;
% vr_dot = vth^2 / r + sin(beta)*T/m - mu/r^2;
% vth_dot = cos(beta)*T/m - vth*vr / r;
% m_dot = -T/ve;
%-----------------------------------------------------------------%
global mu psStuff nstates ncontrols ve
global iGfun jGvar
%-----------------------------------------------------------------%
%             Define the constants for the problem                %
%-----------------------------------------------------------------%
mu      = 1;ve   = 1.8758344;
%-----------------------------------------------------------------%

% Define the bounds on the variables in the optimal controls problem
r0 = 1; th0  = 0; vr0 = 0; vth0 = sqrt(mu/r0); m0 = 1; t0      = 0;
rf = 1.5; vrf = 0; vthf    = sqrt(mu/rf);

rmin    = 0.1;      rmax    = 2;
thmin   = -4*pi;    thmax   = 4*pi;
vrmin   = -50;      vrmax   = +50;
vthmin  = -50;      vthmax  = +50;
mmin    =  0;       mmax    = 1.1;
Tmin    = 0;        Tmax    = 0.1405;
betamin = -4*pi;    betamax = +4*pi;
t0min   = 0;        t0max   = 0;
tfmin   = 0;        tfmax   = 100;
%-----------------------------------------------------------------%
%-----------------------------------------------------------------%
% Define the sizes of quantities in the optimal control problem
nstates     = 5;
ncontrols   = 2;
%-----------------------------------------------------------------%
%-----------------------------------------------------------------%
%      Compute Points, Weights, and Differentiation Matrix        %
%-----------------------------------------------------------------%
%-----------------------------------------------------------------%
% Choose Polynomial Degree and Number of Mesh Intervals           %
%-----------------------------------------------------------------%
N = 4;
numIntervals = 10;
%-----------------------------------------------------------------%
% DO NOT ALTER THE LINE OF CODE SHOWN BELOW!                      %
%-----------------------------------------------------------------%
meshPoints = linspace(-1,1,numIntervals+1).';  
polyDegrees = N*ones(numIntervals,1);
[tau,w,D] = lgrPS(meshPoints,polyDegrees);
psStuff.tau = tau; psStuff.w = w; psStuff.D = D; NLGR = length(w);
%-----------------------------------------------------------------%
% DO NOT ALTER THE LINES OF CODE SHOWN ABOVE!                     %
%-----------------------------------------------------------------%
%-----------------------------------------------------------------%
% Set the bounds on the variables in the NLP.                     %
%-----------------------------------------------------------------%
zrmin       = rmin * ones(length(tau),1);
zrmax       = rmax * ones(length(tau),1);
zrmin(1)    = r0; zrmax(1) = r0;
zrmin(NLGR+1) = rf;  zrmax(NLGR+1) = rf;

zthmin       = thmin * ones(length(tau),1);
zthmax       = thmax * ones(length(tau),1);
zthmin(1)    = th0; zthmax(1) = th0;

zvrmin      = vrmin * ones(length(tau),1);
zvrmax      = vrmax * ones(length(tau),1);
zvrmin(1)   = vr0;  zvrmax(1) = vr0;
zvrmin(NLGR + 1) = vrf; zvrmax(NLGR +1) = vrf;

zvthmin     = vthmin* ones(length(tau),1);
zvthmax     = vthmax* ones(length(tau),1);
zvthmin(1)  = vth0; zvthmax(1) = vth0;
zvthmin(NLGR + 1) = vthf; zvrmax(NLGR +1) = vthf;

zmmin       = mmin * ones(length(tau),1);
zmmax       = mmax* ones(length(tau),1);
zmmin(1)    = m0; zmmax(1) = m0;

zbetamin    = betamin * ones(length(tau)-1,1);
zbetamax    = betamax * ones(length(tau)-1,1);

zTmin       = Tmin * ones(length(tau)-1,1);
zTmax       = Tmax * ones(length(tau)-1,1);

zmin = [zrmin; zthmin; zvrmin ; zvthmin; zmmin; zbetamin; zTmin; t0min; tfmin];
zmax = [zrmax; zthmax; zvrmax; zvthmax; zmmax; zbetamax; zTmax; t0max; tfmax];

%-----------------------------------------------------------------%
% Set the bounds on the constraints in the NLP.                   %
%-----------------------------------------------------------------%
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
pathMin = []; pathMax = [];
eventMin = []; eventMax = [];
objMin = 0; objMax = inf;
Fmin = [objMin; defectMin; pathMin; eventMin];
Fmax = [objMax; defectMax; pathMax; eventMax];

%-----------------------------------------------------------------%
% Supply an initial guess for the NLP.                            %
%-----------------------------------------------------------------%
rguess  = linspace(r0,rf,NLGR+1).';
thguess  = linspace(th0,2*pi,NLGR+1).';
vrguess  = linspace(vr0,vrf,NLGR+1).';
vthguess  = linspace(vth0,vthf,NLGR+1).';
mguess  = linspace(m0,0.1,NLGR+1).';
betaguess = linspace(0,2*pi,NLGR).';
Tguess  = linspace(0,0.2,NLGR).';
t0guess = 0;
tfguess = 10;
z0  = [rguess; thguess; vrguess; vthguess; mguess; betaguess; Tguess; t0guess; tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint Funtction Derivatives
xsize = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('OrbitTransferFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective Funtcion Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('OrbitTransferObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% Set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)OrbitalTransferObj(Z);
funcs.gradient    = @(Z)OrbitalTransferGrd(Z);
funcs.constraints = @(Z)OrbitalTransferCon(Z);
funcs.jacobian    = @(Z)OrbitalTransferJac(Z);
funcs.jacobianstructure = @()OrbitalTransferJacPat(S_jac);
options.ipopt.hessian_approximation = 'limited-memory';

%-----------------------------------------------------------------%
% Set IPOPT Options %
%-----------------------------------------------------------------%
options.ipopt.tol = 1e-5;
options.ipopt.linear_solver = 'ma57';
options.ipopt.max_iter = 2000;
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.ma57_automatic_scaling = 'yes';
options.ipopt.print_user_options = 'yes';
options.ipopt.output_file = ['OrbitalTransfer','IPOPTinfo.txt']; % print output file
options.ipopt.print_level = 5; % set print level default

options.lb = zmin; % Lower bound on the variables.
options.ub = zmax; % Upper bound on the variables.
options.cl = Fmin; % Lower bounds on the constraint functions.
options.cu = Fmax; % Upper bounds on the constraint functions.

%-----------------------------------------------------------------%
% Call IPOPT
%-----------------------------------------------------------------%
[z, info] = ipopt(z0,funcs,options);

%-----------------------------------------------------------------%
% extract lagrange multipliers from ipopt output, info
%-----------------------------------------------------------------%
Fmul = info.lambda;

%-----------------------------------------------------------------%
% Extract the state and control from the decision vector z.
%-----------------------------------------------------------------%
r       = z(1:NLGR+1);
th      = z(NLGR+2:2*(NLGR+1));
vr      = z(2*(NLGR+1)+1:3*(NLGR+1));
vth     = z(3*(NLGR+1)+1:4*(NLGR+1));
m       = z(4*(NLGR+1)+1:5*(NLGR+1));
beta    = z(5*(NLGR+1)+1:6*(NLGR+1)-1);
T       = z(6*(NLGR+1):7*(NLGR+1)-1);
t0      = z(end-1);
tf      = z(end);
t       = (tf-t0)*(tau+1)/2+t0;
tLGR    = t(1:end-1);

%-----------------------------------------------------------------%
% Extract the Lagrange multipliers corresponding to the defect constraint %
%-----------------------------------------------------------------%
multipliersDefects = Fmul(2:nstates*NLGR+1);
multipliersDefects = reshape(multipliersDefects,NLGR,nstates);
%-----------------------------------------------------------------%
% Compute the costates at the LGR points via transformation       %
%-----------------------------------------------------------------%
costateLGR = inv(diag(w))*multipliersDefects;
%-----------------------------------------------------------------%
% Compute the costate at the tau=+1 via transformation            %
%-----------------------------------------------------------------%
costateF = D(:,end).'*multipliersDefects;
%-----------------------------------------------------------------%
% Now assemble the costates into a single matrix                  %
%-----------------------------------------------------------------%
costate = [costateLGR; costateF];
lamr    = costate(:,1);
lamth   = costate(:,2);
lamvr   = costate(:,3);
lamvth  = costate(:,4);
lamm    = costate(:,5);
%-----------------------------------------------------------------%
% plot results
%-----------------------------------------------------------------%
figure(1)
plot(t,r,'-bs',t,th,'-ro',t,vr,'-gd',t,vth,'-m*',t,m,'-c*');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(r(t),th(t),vr(t),vth(t),m(t))$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontName','Times','FontSize',16);
grid on;

% Evaluating the control beta
for i = 1:length(t)
    betai(i) = atan2(lamvr(i),lamvth(i));
end
figure(2)
betai = 180*unwrap(betai)/pi;
plot(t,betai);
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$beta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontName','Times','FontSize',16);
grid on;

figure(3);
plot(t,T,'r*-');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$T(t)$','Interpreter','LaTeX');
set(xl,'FontSize',16);
set(yl,'FontSize',16);
set(gca,'FontName','Times','FontSize',16);
grid on;

figure(4)
polarplot(th,r,'r--')
hold on
polarplot(0:.1:3*pi,ones(length(0:.1:3*pi)))
hold on
polarplot(0:.1:3*pi,1.5*ones(length(0:.1:3*pi)))
title('Optimal polar trajectory for Helio centric Orbital transfer');








