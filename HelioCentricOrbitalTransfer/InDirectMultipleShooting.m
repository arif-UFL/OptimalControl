%% Orbital transfer optimization using In-direct Multiple shooting 
clear all;close all ; clc
tstart = tic;
s = struct; % Defining a Structural array 
% Defining fields and assigning the values to the fields in structural array.
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
s.k       = 16;              k       = s.k; % no of intervals
s.n       = 5;              n       =s.n;   % No of states 
% Initial Costate and final time guesses;
tfguess         = 3.3;
lam_r0_guess    = -2;
lam_th0_guess   = 0;
lam_vr0_guess   = 2;
lam_vth0_guess  = 2;
lam_m0_guess    = 2;
% scaling the time line
tau             = linspace(-1,+1,k+1);
s.tau           = tau;
% initializing and assigning the guess matrix.
p0_guess        = ones(2*n,k-1);
z0guess         = [lam_r0_guess;lam_th0_guess;lam_vr0_guess;lam_vth0_guess;lam_m0_guess;tfguess];
z0              = [z0guess;p0_guess(:)];
% solving the NLP
options         = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',100000);
z               = fsolve(@errorims,z0,options,s);
[E,t,p]         = errorims(z,s);
%% Plotting the result
figure(1)% States Plot
plot(t,p(:,1),'b-',t,p(:,2),'g-',t,p(:,3),'r-',t,p(:,4),'c-',t,p(:,5),'m-');
title ('Optimal states for helio-centric orbital transfer using multiple shooting ');
legend('r','theta','vr','vtheta','m');
grid on;
xlabel('tau');
ylabel('states');
% Evaluating control 
for j = 1:length(t)
    lam_vrj     = p(j,8);
    lam_vthj    = p(j,9);
    betaj(j) = atan2(lam_vrj,lam_vthj);
end
% Control plot
figure(2)
betaj = unwrap(betaj);
plot(t,betaj)
hold on
title ('Optimal control for helio-centric orbital transfer using multiple shooting');
xlabel ('Tau')
ylabel('beta (rad/sec)')
legend('Control beta(tau)');
figure(3)% polar trajectory plot
polarplot(p(:,2),p(:,1),'r--')
hold on
polarplot(0:.1:3*pi,ones(length(0:.1:3*pi)))
hold on
polarplot(0:.1:3*pi,1.5*ones(length(0:.1:3*pi)))
title ('Optimal polar trajectory for helio-centric orbital transfer using multiple shooting');
figure(4)% Costate plots
plot(t,p(:,6:10));
title ('Optimal co-states for helio-centric orbital transfer using multiple shooting ');
legend('lam_r','lam_theta','lam_vr','lam_vtheta','lam_m');
grid on;
xlabel('tau');
ylabel('Co-states');
tend = toc(tstart);
function [E,t,p] = errorims(z0,s)
%% Extracting data from structures
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
n       = s.n; % No. of states
k       = s.k; % No. of stages for shooting
t0      = s.t0;
tau     = s.tau;
tfguess = z0(6);
p0_guess = z0(7:end);
p0_guess = reshape(p0_guess,2*n,k-1);
options = odeset('RelTol',1e-8);
E = [];
p = [];
t = [];
for i = 1:k
    if i == 1
        p0 = [r0;th0;vr0;vth0;m0;z0(1);z0(2);z0(3);z0(4);z0(5)];
    else
        p0 = p0_guess(:,i-1);
    end
    tauspan = [tau(i),tau(i+1)];
    [tout,pout] = ode113(@ODEIDSS,tauspan,p0,options,tfguess,s);
    ptf = pout(end,:).';

    if i < k
        E = [E;ptf - p0_guess(:,i)];
    end

    t = [t; tout];
    p = [p ; pout];
end
r_tf        = ptf(1);
th_tf       = ptf(2);
vr_tf       = ptf(3);
vth_tf      = ptf(4);
m_tf        = ptf(5);
lam_r_tf    = ptf(6);
lam_th_tf   = ptf(7);
lam_vr_tf   = ptf(8);
lam_vth_tf  = ptf(9);
lam_m_tf    = ptf(10);
beta_tf     =atan2(lam_vr_tf,lam_vth_tf);
h_tf            =-1 + lam_r_tf * vr_tf +lam_th_tf*vth_tf / r_tf + lam_vr_tf * (((vth_tf^2)/r_tf) - mu/r_tf^2 + sin(beta_tf)*T/m_tf) + lam_vth_tf*(cos(beta_tf)*T/m_tf  - vth_tf*vr_tf / r_tf) + lam_m_tf*(-T/ve);
% Computing Final error.
E = [E;r_tf - rf; vr_tf - vrtf; vth_tf - vthtf; h_tf;lam_m_tf-1;lam_th_tf;];
end
%% ODE Function
function [pdot] = ODEIDSS(~,p,tf_guess,s)
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
pdot = (tf_guess - t0 /2) .* pdot;
end