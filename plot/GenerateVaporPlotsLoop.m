clc
clear
%close all

We = 85;
surface = struct('type','','pitch',0,'fs',0,'height',0);
surface.type = 'post';
surface.fs = 0.2;
surface.pitch = 8;

r0 = 1.1e-3;
Tw = 120:5:320;
maxV = zeros(1,length(Tw));
for i = 1:length(Tw)
    [t,qvap,qsens,mdot,vapor,z,T,Area] = VaporGenFunc(We,Tw(i),r0,surface);
    maxV(i) = max(vapor);
end


% plot(Tw,maxV,'DisplayName',Surface)
plot(Tw,maxV/.0072,'DisplayName',num2str(We))
hold on

%% Finish Plots

legend()
xlabel('T (^\circC)')
ylabel('Vapor*')
ylim([0 .04])

