function [t,qvap,qsens,mdot,vapor,z,T,Area] = VaporGenFunc(We,Tw,r0,surface)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  Weber number of impacting drop, "We"
%         Temperature of the surface, "Tw", in Celsius
%         initial radius of the drop, "r0", in meters
%         Matlab Struct, "surface", contains surface type 
%                 ('post','rib',or 'smooth'),
%                 surface pitch in microns, surface solid fraction, and 
%                 surface height in microns
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
[lamhat, lamThat, Tmax, tmaxhat, theta] = surface_props(surface,r0*2);

% Fluid Parameters
rho = 998;
nu = 1e-6;
mu = rho*nu;
sigma = 0.07275;
Hv  = 2.26e6;

D0 = 2*r0;
[Re, Oh, V] = We_to_Re(We, rho, mu, sigma, D0);

% Normalized Values
rhat0 = 0.39*2;
lamT = lamThat*D0;

numer = (2*lamhat + rhat0^2/6)*sqrt(324*rhat0^10*We);
denom = sqrt(6*rhat0^6 + 108*lamhat*rhat0^8 + 648*lamhat^2*rhat0^10 + 0.2 ...
    + 6*lamhat*rhat0^2 + 48*lamhat^2*rhat0^7);
drdt = numer/denom;
tc = pi/(2*sqrt(2*(1-cos(theta))))*sqrt(rho*D0^3/sigma);

rhat0 = 0.39;
initialConditions = [rhat0 drdt];
tspan = [1e-5 tc]*V/(D0*sqrt(We));
[t,r] = ode45(@(t,r) radius_v_t(t,r,lamhat,theta,Oh)...
    ,tspan,initialConditions); 

xClav = t*(D0*sqrt(We))/V;
yClav = r(:,1)*D0;

kRois = 1;
tc = .6*D0^2*20/(kRois^2*Re^(4/5)*nu^(4/5));
tmax = tmaxhat*tc;

if We <= 150
    Nt = sum(xClav>tmax & (xClav<tc));

    if tc < tmax + 4e-5
        tmax = tc-3.5e-4;
        Nt = 4;
    end

    t = linspace(tmax+1e-5,tc,Nt*2);
    % rmax = extrapolate_time(yClav,xClav, tmax+.04);
    tcFunc = @(x) x.^(5/2) + .5*x.^(3/2) + .125*x.^(1/2) ... 
        - sqrt(1000*(tc-tmax))*(1000*tc+.25)^2;
    b = fsolve(tcFunc, 0);

    RWet = sqrt(log((sqrt(b)*(b+.25)^2)./ ...
        (sqrt(1000*(t-tmax)).*(1000*t+.25).^2)));
    RWet = real(RWet);
    RWet = RWet*D0;

    xRois = t;
    yRois = RWet;

    P = InterX([xRois; yRois],[xClav'; yClav']);
    
    if isempty(P)
        yRois = yRois*1.1;
        P = InterX([xRois; yRois],[xClav'; yClav']);
    end
    
    yRois = yRois(xRois>P(1));
    xRois = xRois(xRois>P(1));
    yClav = yClav(xClav<P(1));
    xClav = xClav(xClav<P(1));

    
    t = [xClav; xRois'];
    r = [yClav; yRois'];    
else
    tcbar = tc*V/D0;
    tbar = xClav(xClav<=tc)*V/D0;
    Rwet = 0.6*(tbar+0.25).*sqrt(log((sqrt(tcbar)*(tcbar+.25)^2)./ ... 
        (sqrt(tbar).*(tbar+.25).^2)));
    
    t = tbar*D0/V;
    r = Rwet*D0;
end

Area = pi*(r/2).^2;

Nt = length(t);
Nz = 100;
zmax = 5.3e-05;
edr = sqrt(0.598*997*4184);          % Thermal Effusivity of Droplet
ew = sqrt(1300*2329*710);            % Thermal Effusivity of Wall
Tdr = 25;                            % Droplet Temperature
TBoiling = 95;
alpha = 1.4558e-7;                   % Thermal Diffusivity of Liquid
k = 0.598;                           % Thermal Conductivity of Water
kvap = .03;
s = 20e-6;                           % Vapor layer thickness
z = linspace(0, zmax, Nz);           % Z Position
delT = sqrt(12/5*alpha*t);           % Instantaneous Thermal Boundary Layer

Tc = (edr*Tdr + ew*Tw)/(edr + ew);
T = zeros(Nt,Nz);

for i = 1:Nt
    T(i,:) = Tc - 2*lamT*(Tc-Tdr)/(2*lamT+delT(i)) ...
        - 2*(Tc-Tdr)/(2*lamT+delT(i))*z ...
        + (Tc-Tdr)/(delT(i)*(2*lamT+delT(i)))*z.^2;
end

T(diff(T,1,2)>0) = Tdr; 
T(:,end) = Tdr;

vaporpercent = 1./(1+exp(-(Tw-Tmax)/15));
qw_insulated = (Tc-Tdr)./((2*lamT+delT)/(2*k)+s/kvap);
qw_normal = 2*k*(Tc-Tdr)./(2*lamT+delT);
qw = qw_insulated*vaporpercent + qw_normal*(1-vaporpercent);

percent = zeros(Nt,1);
for i = 1:Nt
    
    if i == 1
        t2 = 0;
    else 
        t2 = t1;
        if t2
            top2 = top;
            bottom2 = bottom;
        end
    end
    
    R = T(i,:)-Tdr;
    R(R==0) = [];
    p = polyfit(z(1:length(R)), R, 2);
    x = linspace(0,z(length(R)),500);
    y = p(1)*x.^2 + p(2)*x + p(3);
    yT = zeros(1,length(y))+(TBoiling-Tdr);
    t1 = InterX([x;y],[x;yT]);

    if t1 
        int = trapz(x,y);
        if length(x(x<t1(1))) == 1
            top = 0;
        else
            top = trapz(x(x<t1(1)),y(y>t1(2))-(TBoiling-Tdr));
        end
        bottom = int-top;
        if t2
            percent(i) = (top-top2)/(bottom-bottom2);
        else
            percent(i) = top/bottom;
        end
    else
        percent(i) = 0;
    end
    
end  

qvap = qw.*percent;
qsens = qw.*(1-percent);

Q = qvap.*Area;
mdot = Q/Hv;

if lamThat == 0
    [pks,locs] = findpeaks(mdot);
    loc = locs(pks==max(pks));
    p = polyfitB(t(loc-5:loc+5), mdot(loc-5:loc+5), 2,0);
    mdot(1:loc-5) = p(1)*t(1:loc-5).^2 + p(2)*t(1:loc-5);
end

initial_mass = 4/3*pi*r0^3*rho;
dt = [t(1); diff(t)];
vapor = zeros(1,Nt);
for i = 1:Nt  
    vapor(i) = sum(mdot(1:i).*dt(1:i));
    vapor(i) = vapor(i)/initial_mass;
end

end


