%% Main Function-oribital transfer problem-Indirect Single Shooting
clear all; close all; clc
tstart = tic;
% Defining structural array and assigning fields.
s = struct;
s.r0      = 1;              r0      = s.r0;
s.rf      = 1.5;            rf      = s.rf;
s.th0     = 0;              th0     = s.th0;    % thf = free;   
s.m0      = 1;              m0      = s.m0;     % mf = free;
s.mu      = 1;              mu      = s.mu;
s.T       = 0.1405;         T       = s.T;
s.ve      = 1.8758344;       ve      = s.ve;
s.vr0     = 0;              vr0     = s.vr0;
s.vrtf    = 0;              vrtf    = s.vrtf;
s.vth0    = sqrt(mu/r0);    vth0    = s.vth0;
s.vthtf   = sqrt(mu/rf);    vthtf   = s.vthtf;
s.t0      = 0;              t0      = s.t0;
% Initial Co-state and final time guesses
lam_r0_guess        = -2;
lam_th0_guess       = 0;
lam_vr0_guess       = 2;
lam_vth0_guess      = 2;
lam_m0_guess        = 2;
tfguess             = 3.5;
z0                 = [lam_r0_guess;lam_th0_guess;lam_vr0_guess;lam_vth0_guess;lam_m0_guess;tfguess];
options            = optimset ('MaxFunEvals',100000,'MaxIter',10000,'Display','Iter','TolX',1e-8,'TolFun',1e-8);
z                  = fsolve(@errorIDSS,z0,options,s);
[E,t,p]            = errorIDSS(z,s);
tend = toc(tstart);
% Plotting the results.
figure(1)
plot(t,p(:,1),'b-',t,p(:,2),'g-',t,p(:,3),'r-',t,p(:,4),'c-',t,p(:,5),'m-');
hold on
title('Optimal states for Helio centric orbital transfer using Indirect single shooting');
legend('r','theta','v_r','v_theta','m');
% Computing Control
for j = 1:length(t)
    lam_vrj = p(j,8);
    lam_vthj = p(j,9);
    betaj(j) = atan2(lam_vrj,lam_vthj);
end
figure(2)
betaj = unwrap(betaj);
plot(t,betaj);
title('Optimal control for Helio centric Orbital transfer');
figure(3)
polarplot(p(:,2),p(:,1),'r--')
hold on
polarplot(0:.1:3*pi,ones(length(0:.1:3*pi)))
hold on
polarplot(0:.1:3*pi,1.5*ones(length(0:.1:3*pi)))
title('Optimal polar trajectory for Helio centric Orbital transfer');

%% Errror Function
function [E,t,p] = errorIDSS(z0,s)
%% Extracting data from structural array
r0      = s.r0; 
rf	    = s.rf;
th0     = s.th0; % thf = free;
m0      = s.m0; % mf = free;
mu      = s.mu;
T       = s.T;
ve      = s.ve;
vr0     = s.vr0;
vrtf    = s.vrtf;
vth0    = s.vth0;
vthtf   = s.vthtf;
t0      = s.t0;
tfguess = z0(end);
tspan   = [t0,tfguess] ;
options = odeset('RelTol',1e-8); 
p0      = [r0,th0,vr0,vth0,m0,z0(1),z0(2),z0(3),z0(4),z0(5)].';
[t,p]   = ode113(@ODEIDSS,tspan,p0,options,s);
peval           = p(end,:);
r_tf            = peval(1);
th_tf           = peval(2);
vr_tf           = peval(3);
vth_tf          = peval(4);
m_tf            = peval(5);
lam_r_tf        = peval(6);
lam_th_tf       = peval(7);
lam_vr_tf       = peval(8);
lam_vth_tf      = peval(9);
lam_m_tf        = peval(10);
beta_tf         = atan2(lam_vr_tf,lam_vth_tf);
h_tf            =-1 + lam_r_tf * vr_tf +lam_th_tf*vth_tf / r_tf + lam_vr_tf * (((vth_tf^2)/r_tf) - mu/r_tf^2 + sin(beta_tf)*T/m_tf) + lam_vth_tf*(cos(beta_tf)*T/m_tf  - vth_tf*vr_tf / r_tf) + lam_m_tf*(-T/ve);
% Terminal state matrix
E = [r_tf-rf; vr_tf-vrtf; vth_tf - vthtf; h_tf; lam_m_tf-1; lam_th_tf];
end

%% ODE Function
function [pdot] = ODEIDSS(~,p,s)
t0  = s.t0;
T   = s.T;
mu  = s.mu;
ve  = s.ve;
r   =p(1);
th  =p(2);
vr  =p(3);
vth =p(4);
m   =p(5);
lam_r   = p(6);
lam_th  = p(7);
lam_vr  = p(8);
lam_vth =p(9);
lam_m   = p(10);
beta    = atan2(lam_vr,lam_vth);
%%Dynamics
r_dot = vr;
th_dot = vth / r;
vr_dot = vth^2 / r + sin(beta)*T/m - mu/r^2;
vth_dot = cos(beta)*T/m - vth*vr / r;
m_dot = -T/ve;
lam_r_dot = lam_th*vth/r^2 - lam_vr*(2*mu/r^3 - vth^2/r^2)-lam_vth*(vr*vth/r^2);
lam_th_dot = 0;
lam_vr_dot = -lam_r + lam_vth*vth/r;
lam_vth_dot = -lam_th/r - 2*lam_vr*vth/r + lam_vth*vr/r;
lam_m_dot = T/m^2 * (sin(beta)*lam_vr + cos(beta)*lam_vth); %% sqrt(lam_vr^2 + lam_vth^2);
pdot = [r_dot; th_dot; vr_dot; vth_dot; m_dot; lam_r_dot; lam_th_dot; lam_vr_dot; lam_vth_dot; lam_m_dot;];
end