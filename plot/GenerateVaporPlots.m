clc
clear
%close all

We = 85;
r0 = 1.1e-3;
edr = 1.5794e+03;                       % Thermal Effusivity of Droplet
ew = 1.4662e+03;                        % Thermal Effusivity of Wall
Tw = 220;
Tdr = 25;
TBoiling = 95;
Tc = (edr*Tdr + ew*Tw)/(edr + ew);

surface = struct('type','','pitch',0,'fs',0,'height',0);
surface.type = 'smooth';
surface.fs = 0.1;
surface.pitch = 16;

[t,qvap,qsens,mdot,vapor,z,T,Area] = VaporGenFunc(We,Tw,r0,surface);
initmass = 4/3*pi*r0^3*998;
figure(1)
subplot(2,2,1)
plot(t*1000,qvap)
ylabel('$$\hat{q}"_{lat} (W/m^2)$$','Interpreter','Latex')
hold on
subplot(2,2,2)
plot(t*1000,qsens)
ylabel('$$\hat{q}"_{sens} (W/m^2)$$','Interpreter','Latex')
hold on
subplot(2,2,3)
xlabel('t (ms)')
plot(t*1000,Area*1000^2)
ylabel('Contact Area (mm^2)')
hold on
subplot(2,2,4)
plot(t*1000,mdot)
xlabel('t (ms)')
ylabel('$$\dot{m} (kg/s)$$','Interpreter','Latex')
hold on

figure(2)
plot(t*1000,vapor,'DisplayName','12')
xlabel('normalized time')
ylabel('Vapor*')
legend()
hold on

figure(3)
plotTemp = T;
plotTemp(2:2:end,:) = [];
plotTemp(2:2:end,:) = [];
plotTemp(2:2:end,:) = [];
figure(3)
plot([0 5.25e-2], [TBoiling TBoiling], '--b')
hold on
plot(z*1000, plotTemp, 'k')
legend('Boiling Temperature','Increasing Time')
xlim([0 5.25e-2])
ylim([Tdr Tc])
xlabel('z (mm)')
ylabel('Temperature (C)')
hold off