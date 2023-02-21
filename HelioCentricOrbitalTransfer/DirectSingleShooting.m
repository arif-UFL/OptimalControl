%% Orbital transfer Main function
close all; clear all;clc
tstart = tic;
% Defining structure array and Assigning values to the fields
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
% Fmincon function inputs genration 
n        = 6;   % degree of polynomial used to parameterize the control
cguess  = 0.01*ones(n+1,1);% guess for polynomial coefficients
tfguess = 8;% Final time guess
zguess  = [cguess;tfguess];% decision vector
% setting the FMINCON solver inputs
tfmin    = 0;
tfmax    = 100;
zmin    = [-20*ones(n+1,1);tfmin];
zmax    = [20*ones(n+1,1);tfmax];
A       = [];
B       = [];
Aeq     = [];
Beq     = [];

options1=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5000);
z       = fmincon(@OBJDSS,zguess,A,B,Aeq,Beq,zmin,zmax,@ERRORDSS,options1,s);
c       = z(1:end-1);
[Eineq,Eeq,t,p] = ERRORDSS(z,s);
beta            = polyval(c,t); % evaluating the control beta(t)
%plotting the results
figure(1)
plot(t,beta,'r--')
legend('beta vs t')
% States Plot 
figure(2)
plot(t,p(:,1),'b-',t,p(:,2),'g-',t,p(:,3),'r-',t,p(:,4),'c-',t,p(:,5),'m-');
hold on
title('Optimal states for Helio centric orbital transfer using Direct single shooting');
xlabel('States')
ylabel('t(secs)')
legend('r(t)','theta(t)','vr(t)','vtheta(t)','m(t)')
%Polar trajectory plot
figure(3)
polarplot(p(:,2),p(:,1),'r--')
hold on
polarplot(0:.1:3*pi,ones(length(0:.1:3*pi)))
hold on
polarplot(0:.1:3*pi,1.5*ones(length(0:.1:3*pi)))
title('Optimal polar trajectory for Helio centric Orbital transfer using direct single shooting');
tend = toc(tstart);  % computational time evaluation
%% ODE Function
function [pdot] = ODEDSS(t,p,c,s)
t0  = s.t0;
T   = s.T;
mu  = s.mu;
ve  = s.ve;
% extracting the states from decision vector
r   =p(1);
th  =p(2);
vr  =p(3);
vth =p(4);
m   =p(5);
beta    = polyval(c,t);% evaluating the initial Control
% Dynamics
r_dot = vr;
th_dot = vth / r;
vr_dot = (vth^2 / r) + sin(beta)*T/m - mu/r^2;
vth_dot = cos(beta)*T/m - vth*vr / r;
m_dot = -T/ve;
pdot = [r_dot; th_dot; vr_dot; vth_dot; m_dot];
end

%% Objective function 
% optimizing the final time
function J = OBJDSS(z,s)
tf = z(end);
J = tf;
end

%% Error fn
function [Eineq,Eeq,t,p] = ERRORDSS(z,s)
c             = z(1:end-1);% extracting co-efficient guesses from decision vector 
tf_guess      = z(end);% final guess extraction
% Extracting data from Structure array
r0      = s.r0; 
rf	    = s.rf;
th0     = s.th0; % thf = free;
m0      = s.m0; % mf = free;
mu      = s.mu;
T       = s.T;
ve      = s.ve;
vr0     = s.vr0;
vrf    = s.vrtf;
vth0    = s.vth0;
vthf   = s.vthtf;
t0      = s.t0;
p0      =[r0;th0;vr0;vth0;m0];% initial states
tspan   =[t0,tf_guess];% Time span for integration
options =odeset('RelTol',1e-8);% ODE113 solver tolerence setup
[t,p]   =ode113(@ODEDSS,tspan,p0,options,c,s);% Integrating the system dynamics
% Extracting evaluated values
r_tf    =p(end,1);
th_tf   =p(end,2);
vr_tf   =p(end,3);
vth_tf  =p(end,4);
m_th    =p(end,5);
% Composing Error
Eeq     = [r_tf-rf;vr_tf - vrf;vth_tf - vthf];
Eineq   = [];
end
