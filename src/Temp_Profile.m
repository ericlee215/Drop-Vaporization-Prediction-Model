function [z, T, qw] = Temp_Profile(t,lamThat,d,s,Tmax)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  time, "t", in seconds
%         Non-Dimensional temperature jump length, "lamThat"
%         Matlab Struct, "d", short for "drop" contains:
%                 Drop initial radius, "r0", in meters
%                 Drop initial temperature, "Tdr", in Celsius
%         Matlab Struct, "s", short for "surface" contains:
%                 surface Temperature, "Tw", in Celsius
%         Temperature of maximum atomization, "Tmax"
%
% Returns Vapor Produced normalized by the initial droplet mass, "Vstar"
%         mass flow rate of vapor generation "mdot"
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

Nt = length(t);
Nz = 100;
zmax = (4*d.r0)*10^-2;
edr = sqrt(0.598*997*4184);          % Thermal Effusivity of Droplet
ew = sqrt(1300*2329*710);            % Thermal Effusivity of Wall
alpha = 1.4558e-7;                   % Thermal Diffusivity of Liquid
k = 0.598;                           % Thermal Conductivity of Water
kvap = .03;
Vap_thick = 20e-6;                           % Vapor layer thickness
z = linspace(0, zmax, Nz);           % Z Position
delT = sqrt(12/5*alpha*t);           % Instantaneous Thermal Boundary Layer

lamT = lamThat*2*d.r0;

Tc = (edr*d.Tdr + ew*s.Tw)/(edr + ew);
T = zeros(Nt,Nz);

for i = 1:Nt
    T(i,:) = Tc - 2*lamT*(Tc-d.Tdr)/(2*lamT+delT(i)) ...
        - 2*(Tc-d.Tdr)/(2*lamT+delT(i))*z ...
        + (Tc-d.Tdr)/(delT(i)*(2*lamT+delT(i)))*z.^2;
end

T(diff(T,1,2)>0) = d.Tdr; 
T(:,end) = d.Tdr;

vaporpercent = 1./(1+exp(-(s.Tw-Tmax)/15));
qw_insulated = (Tc-d.Tdr)./((2*lamT+delT)/(2*k)+Vap_thick/kvap);
qw_normal = 2*k*(Tc-d.Tdr)./(2*lamT+delT);
qw = qw_insulated*vaporpercent + qw_normal*(1-vaporpercent);

end
