function [qvap, qsens] = separate_heat_flux(z,T,qw,d)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Inputs  height into drop, "z", in meters
%         drop temperature profile, "T", in Celsius
%         heat flux, "qw", in W/m^2
%         Matlab Struct, "d", short for "drop" contains:
%                 Drop initial temperature, "Tdr", in Celsius
%
% Returns latent heat flux, "qvap", in W/m^2
%         sensible heat flux, "qsens", in W/m^2
%
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

TBoiling = 95;
percent = zeros(length(qw),1);
for i = 1:length(qw)
    
    if i == 1
        t2 = 0;
    else 
        t2 = t1;
        if t2
            top2 = top;
            bottom2 = bottom;
        end
    end
    
    R = T(i,:)-d.Tdr;
    R(R==0) = [];
    x = linspace(0,z(length(R)),500);
    yT = zeros(1,length(x))+(TBoiling-d.Tdr);
    if length(R)>2
        p = polyfit(z(1:length(R)), R, 2);   
        y = p(1)*x.^2 + p(2)*x + p(3);
    else
        y = zeros(1,length(x));
    end
    
    t1 = InterX([x;y],[x;yT]);

    if t1 
        int = trapz(x,y);
        if length(x(x<t1(1))) == 1
            top = 0;
        else
            top = trapz(x(x<t1(1)),y(y>t1(2))-(TBoiling-d.Tdr));
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

end
