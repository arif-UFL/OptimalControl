% Direct multiple shooting
close all; clear all; clc
tstart = tic;
% Definig the structral array, creating and assigning the fields
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
% Direct Multiple shooting
nx      =     5; % no. of states
s.nx    =     nx;
n       =     6; % degree of polynomial.
s.n     =     n ;
k       =     16; % no. of intervals
s.k     =     k;
tau     =     linspace(-1,+1,k+1);% scaling timeline
s.tau   =     tau;
cguess  =     zeros(k*(n+1),1);% co-efficient guesses
pguess  =     ones(nx,k-1);% initial state guesses
tf_guess  =   3.5;% final time guess
zguess  =     [cguess;pguess(:);tf_guess];% guesses coloumn matrix.
% Inputs for FMINCON
tfmin   =     0;
tfmax   =     20;
pgmin   =     -50*ones(nx,k-1);
pgmax   =     +50*ones(nx,k-1);
zmin    =     [-50*ones(k*(n+1),1);pgmin(:);tfmin];
zmax    =     [+50*ones(k*(n+1),1);pgmax(:);tfmax];

A       =     [];
B       =     [];
Aeq     =     [];
Beq     =     [];

options1 = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-4,'MaxFunEvals',100000,'MaxIter',1000);
z   =   fmincon(@OBJDMS,zguess,A,B,Aeq,Beq,zmin,zmax,@ERRORDMS,options1,s);
[Eineq,Eeq,t,p] = ERRORDMS(z,s);
% plotting the results
figure(1)
% evaluating control sub-arcs
c    = z(1:k*(n+1));
c    = reshape(c,n+1,k);
ti   = t(1);
j   =   1;
for i = 1:k
    while ti(end) < tau(i+1)
        ti = [ti,t(j)];
        j = j+1;
    end
    beta = polyval(c(:,i),ti);
    plot(ti,beta,'r--');
    hold on
    if j <= length(t)
        ti = t(j);
    end
title ('Optimal control of orbital transfer using direct multiple shooting ')
xlabel ('Tau')
ylabel('beta (rad/sec)')
end
figure(2) % state plots
plot(t,p(:,1),'b-',t,p(:,2),'g-',t,p(:,3),'r-',t,p(:,4),'c-',t,p(:,5),'m-');
hold on
title('Optimal states for Helio centric orbital transfer using Direct Multiple shooting');
legend('r','theta','v_r','v_theta','m');
xlabel('States')
ylabel('tau')
figure(3)% Polar trajectory plot
polarplot(p(:,2),p(:,1),'r--')
hold on
polarplot(0:.1:3*pi,ones(length(0:.1:3*pi)))
hold on
polarplot(0:.1:3*pi,1.5*ones(length(0:.1:3*pi)))
title ('Optimal polar trajectory for orbital transfer using direct multiple shooting ')
tend = toc(tstart);% Computational time evalustion
% Error Function
function [Eineq,Eeq,t,P] = ERRORDMS(z,s)
% extracting data from structural array
nx      = s.nx;
n       = s.n;
k       = s.k;
r0      = s.r0;
rf      = s.rf;
th0     = s.th0;
m0      = s.m0;
t0      = s.t0;
vrtf    = s.vrtf;
vr0     = s.vr0;
vth0    = s.vth0;
vthtf   = s.vthtf;
k       = s.k;
tau     = s.tau;
% Guess arrays generation and Matrix reshaping
cguess  =   z(1:end-nx*(k-1)-1);
cguess  =   reshape(cguess,n+1,k);
p0_guess =  z(end-nx*(k-1):end-1);
p0_guess = reshape(p0_guess,nx,k-1);
tf_guess = z(end);
% Generating empty arrays for E,t,p
E = [];
t = [];
P = [];
% Direct multiple shooting process
options = odeset('RelTol',1e-8);
for i = 1:k
    if i == 1
        P0  = [r0;th0;vr0;vth0;m0];
    else
        P0  = p0_guess(:,i-1);
    end
    c = cguess(:,i);
    tauspan     =   [tau(i),tau(i+1)];
    [tout,pout] =   ode113(@ODEDMS,tauspan,P0,options,c,s,tf_guess);
    ptf     =   pout(end,:).';
    if i < k
        E = [E;ptf-p0_guess(:,i)];
    end
    t   =   [t; tout];
    P   =   [P; pout];
end
r_tf    = ptf(1);
th_tf   = ptf(2);
vr_tf   = ptf(3);
vth_tf  = ptf(4);
m_tf    = ptf(5);
Eeq     = [E;r_tf-rf;vr_tf-vrtf;vth_tf-vthtf];% Terminal state matrix
Eineq  =   [];
end
% Objective function
% Objective is to minimize the final time which is same as maximizing mass
% of space craft
function J = OBJDMS(z,s)
tf = z(end);
J   = tf;
end
% ODE Function
function pdot = ODEDMS(t,p,c,s,tf_guess)
t0  = s.t0;
T   = s.T;
mu  = s.mu;
ve  = s.ve;
% Extracting initial state values
r   =p(1);
th  =p(2);
vr  =p(3);
vth =p(4);
m   =p(5);
% control evaluation
beta    = polyval(c,t);
% Dynamics
r_dot = vr;
th_dot = vth / r;
vr_dot = (vth^2 / r) + sin(beta)*T/m - mu/r^2;
vth_dot = cos(beta)*T/m - vth*vr / r;
m_dot = -T/ve;
pdot =(tf_guess-t0)/2 * [r_dot; th_dot; vr_dot; vth_dot; m_dot];
end
