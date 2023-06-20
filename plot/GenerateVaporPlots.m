clc
clear
%close all


% Fluid Parameters
fluid = struct('rho',998,'nu',1e-6,'mu',9.98e-4, ... 
           'sigma',0.07275,'Hv',2.26e6);
% Drop Parameters
drop = struct('r0',1.1e-3,'We',85,'Tdr',25);
% Surface Parameters
surface = struct('Tw',220,'type','post',...
                'pitch',24,'fs',0.1,'height',4);

% Main Function
[t,qvap,qsens,mdot,Vstar,z,T,Area] = VaporGenFunc(drop,surface);

% Normalizing Components
[Re, Oh, V] = We_to_Re(drop,fluid);
tc = 48*drop.r0^2/(Re^(4/5)*(fluid.nu)^(4/5));
initmass = 4/3*pi*drop.r0^3*998;

%% Plots

figure(1)
subplot(2,2,1), hold on
plot(t*1000,qvap)
ylabel('$$\hat{q}"_{lat} (W/m^2)$$','Interpreter','Latex')
subplot(2,2,2), hold on
plot(t*1000,qsens)
ylabel('$$\hat{q}"_{sens} (W/m^2)$$','Interpreter','Latex')
subplot(2,2,3), hold on
xlabel('t (ms)')
plot(t*1000,Area*1000^2)
ylabel('Contact Area (mm^2)')
subplot(2,2,4), hold on
plot(t*1000,mdot)
xlabel('t (ms)')
ylabel('$$\dot{m} (kg/s)$$','Interpreter','Latex')

figure(2), hold on
plot(t/tc,Vstar)
xlabel('t/t_c')
ylabel('V^*')
legend()
hold on

figure(3)
subplot(1,2,1), hold on
plot(t/tc,zeros(1,length(t))+0.512)
plot(t/tc,qvap./(qvap+qsens))
xlabel('t/t_c')
ylabel('q"_{lat}/q"_w')
subplot(1,2,2), hold on
plot(t/tc,zeros(1,length(t))+(1-0.512))
plot(t/tc,qsens./(qvap+qsens))
xlabel('t/t_c')
ylabel('q"_{sens}/q"_w')

TBoiling = 95;
plotTemp = T;
plotTemp(2:2:end,:) = [];
plotTemp(2:2:end,:) = [];
plotTemp(2:2:end,:) = [];
figure(4), hold on
plot([0 5.25e-2], [TBoiling TBoiling], '--b')
plot(z*1000, plotTemp, 'k')
legend('Boiling Temperature','Increasing Time')
xlim([0 5.25e-2])
ylim([drop.Tdr max(max(plotTemp))])
xlabel('z (mm)')
ylabel('Temperature (C)')
hold off
