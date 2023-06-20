function [Vstar, mdot] = vapor_from_heat_flux(t,Area,qvap,d,f,lamThat)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  time, "t", in seconds
%         Area, "Area", in m^2
%         Latent Heat flux, "qvap", in W/m^2
%         Matlab Struct, "d", short for "drop" contains:
%                 Drop initial radius, "r0", in meters
%         Matlab Struct, "f", short for "fluid" contains:
%                 Fluid density, "rho", in kg/m^3
%                 Fluid Heat of Vaporization, "Hv", in 
%         Non-Dimensional temperature jump length, "lamThat"
%
% Returns Vapor Produced normalized by the initial droplet mass, "Vstar"
%         mass flow rate of vapor generation "mdot"
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

Q = qvap.*Area;
mdot = Q/f.Hv;

if lamThat == 0
    [pks,locs] = findpeaks(mdot);
    loc = locs(pks==max(pks));
    p = polyfitB(t(loc-5:loc+5), mdot(loc-5:loc+5), 2,0);
    mdot(1:loc-5) = p(1)*t(1:loc-5).^2 + p(2)*t(1:loc-5);
end

initial_mass = 4/3*pi*d.r0^3*f.rho;
dt = [t(1); diff(t)];
Vstar = zeros(1,length(qvap));
for i = 1:length(qvap)  
    Vstar(i) = sum(mdot(1:i).*dt(1:i));
    Vstar(i) = Vstar(i)/initial_mass;
end

end
