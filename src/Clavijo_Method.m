function [xClav, yClav] = Clavijo_Method(f,D0,We,V,Oh,lamhat,theta)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  Matlab Struct, "f", short for "fluid" contains:
%                 Fluid density, "rho"
%                 Fluid kinematic viscocity, "nu"
%                 Fluid dynamic viscocity, "mu"
%                 Fluid surface tension, "sigma
%         Drop diameter, "D0", in meters
%         Drop impact Weber number, "We"
%         Drop impact Velocity, V, in m/s
%         Drop impact Ohnesorge number, Oh
%         Non-dimensional slip length, lamhat
%         Static contact angle, "theta"
%
% Returns Time for Clavijo method, "xClav", in seconds
%         Radius for Clavijo method, "yClav", in meters
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

rhat0 = 0.39*2;

numer = (2*lamhat + rhat0^2/6)*sqrt(324*rhat0^10*We);
denom = sqrt(6*rhat0^6 + 108*lamhat*rhat0^8 + 648*lamhat^2*rhat0^10 + 0.2 ...
    + 6*lamhat*rhat0^2 + 48*lamhat^2*rhat0^7);
drdt = numer/denom;
tc = pi/(2*sqrt(2*(1-cos(theta))))*sqrt(f.rho*D0^3/f.sigma);

rhat0 = 0.39;
initialConditions = [rhat0 drdt];

tspan = [1e-5 tc]*V/(D0*sqrt(We));
[t,r] = ode45(@(t,r) radius_v_t(t,r,lamhat,theta,Oh)...
    ,tspan,initialConditions); 

xClav = t*D0*sqrt(We)/V;
yClav = r(:,1)*D0;
end
