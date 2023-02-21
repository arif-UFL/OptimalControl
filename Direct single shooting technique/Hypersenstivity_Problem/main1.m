clear all, clc;
x0=1;
xf=1;
t0=0;
tf=[10:10:50];

for(i = 1:5)
options = optimset('display','Iter','TolX',1e-8,'TolFun',1e-8);
lambda0 = fsolve(@myerror1,0.1,options,x0,t0,xf,tf(i));
end

y0 = [x0;lambda0]
options = odeset('RelTol',1e-8);
tspan = [t0,tf];
[t,F] = ode113(@myode1,tspan,y0,options);

figure(1)
plot(t,F(:,1))
hold on
title('Trajectory x(t) Numerical optimisation' )
xlabel('Time(seconds)')
legend('Numerical solution trajectory');

 
figure(2)
plot(t,F(:,2))
title('control u(t) Numerical optimisation' )
xlabel('Time(seconds)')
legend('Numerical solution control');
 
