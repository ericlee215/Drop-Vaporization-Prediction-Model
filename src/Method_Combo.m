function [t,r] = Method_Combo(tmax,tc,xClav,yClav,d,V)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  time of maximum atomizatino, "tmax", in seconds
%         Contact time, "tc", in seconds
%         Time for Clavijo method, "xClav", in seconds
%         Radius for Clavijo method, "yClav", in meters
%         Matlab Struct, "d", short for "drop" contains:
%                 Drop impact Weber number, "We"
%                 Drop initial radius, "r0"
%         Drop impact Velocity, V, in m/s
%
% Returns Time for Combined method, "t", in seconds
%         Radius for Combined method, "r", in meters
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if d.We <= 150
    Nt = sum(xClav>tmax & (xClav<tc));

    if tc < tmax + 4e-5
        tmax = tc-3.5e-4;
        Nt = 4;
    end

    t = linspace(tmax+1e-5,tc,Nt*2);
    tcFunc = @(x) x.^(5/2) + .5*x.^(3/2) + .125*x.^(1/2) ... 
        - sqrt(1000*(tc-tmax))*(1000*tc+.25)^2;
    b = fsolve(tcFunc, 0,optimset('Display','off'));

    RWet = sqrt(log((sqrt(b)*(b+.25)^2)./ ...
        (sqrt(1000*(t-tmax)).*(1000*t+.25).^2)));
    RWet = real(RWet);
    RWet = RWet*d.r0*2;

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
    tcbar = tc*V/(2*d.r0);
    tbar = xClav(xClav<=tc)*V/(2*d.r0);
    Rwet = 0.6*(tbar+0.25).*sqrt(log((sqrt(tcbar)*(tcbar+.25)^2)./ ... 
        (sqrt(tbar).*(tbar+.25).^2)));
    
    t = tbar*2*d.r0/V;
    r = Rwet*2*d.r0;
end

end
