function E = myerror1(lambda0, x0, xf,t0,tf)
    p0 = [x0;lambda0];
    options = odeset('RelTol',1e-8);
    tspan = [t0,tf];
    [t,y] = ode113(@myode1,tspan,p0,options);
    ytf = y(end,:);
    xtf = ytf(1,:);
    E = xtf - xf;
end
    