function [lamhat, lamThat, Tmax, tmax, theta] = surface_props(surface,D0)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  initial diameter of the drop, "r0", in meters
%         Matlab Struct, "surface", contains surface type 
%                 ('post','rib',or 'smooth'),
%                 surface pitch in microns, surface solid fraction, and 
%                 surface height in microns
% 
% Returns Nondimensional slip length, "lamhat"
%         Nondimensional temperature jump length, "lamThat"
%         Temperature of max atomization, "Tmax", in Celsius
%         Nondimensional time of max atomization, "tmax", normalized by
%               Roisman contact time.
%         Static contact angle, "theta", in radians 
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

w = surface.pitch;
Fs = surface.fs;
Fc = 1-Fs;

A = .2066;
B = .2339;
C = .1756;
D = .1156;

rhocpSi = 2329*700;
kSi = 120;
rhocpAir = 1.29*1006;
kAir = 0.026;
eFluid = sqrt(998*4200*0.6);

if strcmp(surface.type,'smooth')
    Tmax = 240;
    tmax = 3.5/5.1 ;
    theta = 115*pi/180;
    lam = 0;
    lamT = 0;
else
    if strcmp(surface.type,'post')
        lam = w*(A/sqrt(Fs) - B);
        lamT = w*(C/sqrt(Fs) - D);
    elseif strcmp(surface.type,'rib')
        lam = w*(log(sec(Fc*pi/2))/(2*pi)*(0.19+(2.6e5)/ ...
            ((0+540)^2+(2.1e4))) + log(sec(Fc*pi/2))/pi)/2;
        lamT = w*(log(sec(Fc*pi/2))/pi ...
            + 1.285e+13*exp(-((Fc-4.382)/0.6191).^2) ...
            + 0.7415*exp(-((Fc-1.196)/0.6284).^2))/2;
    end
    rhocpWall = rhocpSi*Fs + rhocpAir*Fc;
    kWall = kSi*Fs + kAir*Fc;
    eWall = sqrt(rhocpWall*kWall);
    contactTmax = 80 + 120/(1+exp((lamT-4.155)/1.65)^1.677);
    Tmax = (contactTmax*(eWall+eFluid)-eFluid*25)/eWall;
    %tmax = 2.66 + 1/(0.5932+exp(-lamT+3.672)^2.919);
    tmax = 0.5227 + 1/(3.023+exp(-lamT+4.259)^2.799);
    tmax = tmax*1e-3;
    theta = 150*pi/180;
end

lamhat = lam*10^-6/D0;
lamThat = lamT*10^-6/D0;

end