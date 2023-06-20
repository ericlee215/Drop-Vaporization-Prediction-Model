function [Re, Oh, V] = We_to_Re(d,f)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  Matlab Struct, "d", short for "drop" contains:
%                 Drop impact Weber number, "We"
%                 Drop initial radius, "r0"
%         Matlab Struct, "f", short for "fluid" contains:
%                 Fluid density, "rho"
%                 Fluid kinematic viscocity, "nu"
%                 Fluid dynamic viscocity, "mu"
%                 Fluid surface tension, "sigma
%
%       (Verify consistent Units)
%                 
% Returns Impact Reynolds number, "Re" 
%         Impact Ohnesorge number, "Oh"
%         Impact Velocity, "V"
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

V = sqrt(d.We*f.sigma/(f.rho*2*d.r0));

Re = f.rho*V*2*d.r0/f.mu;

Oh = sqrt(d.We)/Re;

end
