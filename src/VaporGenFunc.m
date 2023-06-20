function [t,qvap,qsens,mdot,Vstar,z,T,Area] = VaporGenFunc(d,s)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  Matlab Struct, "d", short for "drop" contains:
%                 Drop impact Weber number, "We"
%                 Drop initial radius, "r0", in meters
%                 Drop initial temperature, "Tdr", in Celsius
%         Matlab Struct, "s", short for "surface" contains:
%                 surface Temperature, "Tw", in Celsius
%                 surface type ('post','rib',or 'smooth'),
%                 surface pitch, "pitch", in microns
%                 surface solid fraction, "Fs"
%                 surface height, "h", in microns
% 
% Returns time, "t", in seconds 
%         latent heat, "qvap", in W/m2
%         sensible heat, "qsens", in W/m2
%         mass flow rate of vapor generation, "mdot", in kg/s
%         Nondimensional mass of vapor normalized by initial drop mass,
%               "vapor" 
%         height into drop, "z", in meters
%         Temperature profile at various time steps, "T", in Celsius
%         Contact Area, "Area", in m^2
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

%  Find surface props
[lamhat, lamThat, Tmax, tmaxhat, theta] = surface_props(d,s);

% Fluid Parameters
f = struct('rho',998,'nu',1e-6,'mu',9.98e-4, ... 
           'sigma',0.07275,'Hv',2.26e6);

[Re, Oh, V] = We_to_Re(d,f);

[xClav, yClav] = Clavijo_Method(f,2*d.r0,d.We,V,Oh,lamhat,theta);

tc = 48*d.r0^2/(Re^(4/5)*f.nu^(4/5));

tmax = tmaxhat*tc;

[t,r] = Method_Combo(tmax,tc,xClav,yClav,d,V);

Area = pi*(r/2).^2;

[z, T, qw] = Temp_Profile(t,lamThat,d,s,Tmax);

[qvap, qsens] = separate_heat_flux(z,T,qw,d);

[Vstar, mdot] = vapor_from_heat_flux(t,Area,qvap,d,f,lamThat);

end
