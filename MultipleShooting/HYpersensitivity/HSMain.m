clear all; close all; clc
%% HS Main function bonus 2
% Given
x0 = 1;
xf=1.5;
t0 = 0;
tf = 50;
nx = 1; %(only lambda0 is variable)
k = 100; %(No. of intervals)
tau = linspace(-1,+1,k+1);
lambda0guess =0;
pRemint0guess = zeros(2*nx,k-1);
zguess = [lambda0guess;pRemint0guess(:)];
options = optimset('TolX',1e-8,'TolFun',1e-8,'display','Iter','MaxFunEvals',10000,'MaxIter',10000);
z = fsolve(@HSError1,zguess,options,x0,t0,xf,tf,nx,k,tau);
[E,t,y] = HSError1(z,x0,t0,xf,tf,nx,k,tau);


figure(1)
plot(t,y(:,1),'r--');
xlabel('tau');
ylabel('x(tau)');

figure(2)
plot(t,-y(:,2),'b--');
xlabel('tau');
ylabel('lambda(tau)')


%% Error Function
function [E,t,y] = HSError1(z,x0,t0,xf,tf,nx,k,tau)
lambda0 = z(1);
pRemint = z(2:end);
pRemint = reshape(pRemint,2*nx,k-1);% Reshaping      2 X k-1 matrix.
options = odeset('RelTol',1e-8);
E = [];% alocating memory for E
t = [];% alocating memory for t
y = [];% alocating memory for y
% assigning the guess values for x0 and lambda0guess at all the intervals
for i = 1:k
    if i == 1
        y0 = [x0;lambda0];
    else
        y0 = pRemint(:,i-1);
    end
    tauspan = [tau(i),tau(i+1)];
    [tout,yout] = ode113(@HSmyode1,tauspan,y0,options,t0,tf);
    ytf = yout(end,:)';

    
    if i < k
        E = [E;ytf-pRemint(:,i)];
    end
    t = [t;tout];
    y = [y ;yout];
end
E = [E;ytf(1) - xf];% Computing final error.
end

%% Dynamics Function
function dydt = HSmyode1(t,y,t0,tf)
dydt = zeros(2,1);
dydt(1) = -y(1)-y(2);
dydt(2) = -y(1)+y(2);
dydt = (tf-t0)/2 * dydt;%% For scaling it between -1 and +1 timespan
end