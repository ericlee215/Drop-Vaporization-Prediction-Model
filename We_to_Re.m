function [Re, Oh, V] = We_to_Re(We, rho, mu, sigma, D)

V = sqrt(We*sigma/(rho*D));

Re = rho*V*D/mu;

Oh = sqrt(We)/Re;

end