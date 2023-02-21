%%myode2 function
function pdot = myode2(t,p,g)
thetaguess = 1;
v = p(3);
lambdax = p(4);
lambday = p(5);
lambdav = p(6);
options = optimset('Display','off','Tolx',1e-8,'TolFun',1e-8);
theta = fsolve(@SolveControl,thetaguess,options,lambdax,lambday,lambdav,v,g);
xdot = v*sin(theta);
ydot = v*cos(theta);
vdot = g*cos(theta);
lambdaxdot = 0;
lambdaydot = 0;
lambdavdot = -lambdax*sin(theta)-lambday*cos(theta);
pdot = [xdot;ydot;vdot;lambdaxdot;lambdaydot;lambdavdot];
end