function output = radius_v_t(t,r,lam,theta,Oh)
    output = zeros(2,1);
    s = 1.41*Oh^(-2/3);
    
    R = r(1);
    dr = r(2);
    
    first = 2*R*dr*(1-cos(theta)) - dr/(3*R^2);
    second = Oh/3*(lam + 1/(12*R^2))^-2*(1/(18*R^6) + lam/R^4 ... 
        + 6*lam^2/R^2 + 1/4 + s/(12*R^3))*dr^2;
    a = 1/R^10;
    b = (2*lam + 1/(6*R^2))^-2;
    c = dr^2;
    d = 6*R^6 + 108*lam*R^8 + 648*lam^2*R^10 + 1/5 ...
        + 6*lam*R^2 + 48*lam^2*R^7;
    da = -10*dr/R^11;
    db = -2*(2*lam + 1/(6*R^2))^-3 * (-dr/(3*R^3));
    dc = 2*dr;
    dd = 36*R^5*dr + 864*lam*R^7*dr + 6480*lam^2*R^9*dr ...
        + 12*lam*R*dr + 336*lam^2*R^6*dr;
    third = (a*b*c*dd + a*db*c*d + da*b*c*d)/3888;
    fourth = -a*b*dc*d/3888;
    
    output(1) = dr;
    output(2) = (first + second + third)/fourth;


end
