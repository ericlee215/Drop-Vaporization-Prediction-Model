clc
clear
close all

%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Contains sections that run loops for plotting maximum vapor generation
%   while varying three different parameters: 
%   surface temperature, drop radius, and surface temperature jump length
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

%% Fluid Parameters
% Fluid Parameters
fluid = struct('rho',998,'nu',1e-6,'mu',9.98e-4, ... 
           'sigma',0.07275,'Hv',2.26e6);

%% Vary Temperature

Tw = 120:5:320;
maxV = zeros(1,length(Tw));

% Drop Parameters
drop = struct('r0',1.1e-3,'We',85,'Tdr',25);
% Surface Parameters
surface = struct('Tw',0,'type','post',...
                'pitch',24,'fs',0.42,'height',4);

for i = 1:length(Tw)
    surface.Tw = Tw(i);
    [t,qvap,qsens,mdot,Vstar,z,T,Area] = VaporGenFunc(drop,surface);
    maxV(i) = max(Vstar);
end


figure(1), hold on
plot(Tw,maxV/.0072,'HandleVisibility','off')
plot(Tw(maxV == max(maxV)),max(maxV)/.0072)
xlabel('{\it T_w} (^\circC)')
ylabel('\it V^*_{max}')
xlim([Tw(1) Tw(end)])

%% Vary Diameter

r0 = linspace(0.5e-3,5e-3,25);

% Drop Parameters
drop = struct('r0',0,'We',85,'Tdr',25);
% Surface Parameters
surface = struct('Tw',220,'type','post',...
                'pitch',16,'fs',0.1,'height',4);

maxV = zeros(1,length(r0));
maxVap = maxV;
initial_mass = maxV;
Oh = maxV;
for i = 1:length(r0)
    drop.r0 = r0(i);
    [t,qvap,qsens,mdot,Vstar,z,T,Area] = VaporGenFunc(drop,surface);
    maxV(i) = max(Vstar);
    initial_mass(i) = 4/3*pi*r0(i)^3*998;
    maxVap(i) = maxV(i)*initial_mass(i);
    [~, Oh(i), ~] = We_to_Re(drop,fluid);
end

figure(2), hold on
plot(r0*1000,maxV)
xlabel('r_0 (mm)')
ylabel('V^*_{max}')
figure(3), hold on
subplot(1,2,1), hold on
plot(r0*1000,maxVap)
set(gca,'YScale','log')
xlabel('r_0 (mm)')
ylabel('V (kg)')
xlim([r0(1) r0(end)]*1000)
subplot(1,2,2), hold on
plot(1./Oh,maxV)
xlabel('Oh^{-1}')
ylabel('V^*_{max}')
xlim([Oh(1) Oh(end)].^-1)
figure(4), hold on
plot(1./Oh,maxV)
xlabel('Oh^{-1}')
ylabel('V^*_{max}')
xlim([Oh(1) Oh(end)].^-1)

%% Vary lamT
% This Section Runs multiple loops and will take a while. Decrease n for
%   faster run with lower resolution

n = 50;
w = linspace(1e-3,24,n);
C = .1756;
D = .1156;
Tw = 120:5:320;

% Drop Parameters
drop = struct('r0',1.1e-3,'We',100,'Tdr',25);
% Surface Parameters
surface = struct('Tw',220,'type','post',...
                'pitch',2,'fs',0.1,'height',4);

maxV = zeros(1,length(Tw));
maxmaxV = zeros(1,length(w));
lamT = maxmaxV;
for i = 1:length(w)
    surface.pitch = w(i);
    for j = 1:length(Tw)
        surface.Tw = Tw(j);
        [~,~,~,~,Vstar,~,~,~] = VaporGenFunc(drop,surface);
        maxV(j) = max(Vstar);
    end
    maxmaxV(i) = max(maxV);
    lamT(i) = w(i)*(C/sqrt(surface.fs) - D);
    % i       % I like to have this uncommented to see my progress
end

figure(4), hold on
plot(lamT,maxmaxV,'HandleVisibility','off')
plot(lamT(9),maxmaxV(9),'DisplayName',num2str(drop.We))
xlabel('{\it \lambda_T} (\mum)')
ylabel('{\it V^*_m}')
lgd = legend;
title(lgd,'\it We')
