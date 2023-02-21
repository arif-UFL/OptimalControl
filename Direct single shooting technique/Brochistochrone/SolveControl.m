
function Htheta = SolveControl(theta, lambdax, lambday, lambdav, v, g)
Htheta = lambdax*v*cos(theta)-(lambday*v + lambdav*g)*sin(theta);
% Htf = lambdax*v*sin(theta) + (lambday*v + lambdav*g)*cos(theta);
end