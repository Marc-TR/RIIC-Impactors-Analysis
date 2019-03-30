function [r0,v0]=IC_from_KEP(a, e, i, RAAN, AOP, f0) % Semi major axis in m and angles in radians

% It provides a position vector [m] and a velocity vector [m/s] at the
% given keplerian elements

mu_Sun=1.32712428e20; %m^3/s^2

R0=(a*(1-e^2))/(1+e*cos(f0)); % R at f0

% P vector
Lat=asin(sin(i)*sin(AOP)); % Latitude
if sin(Lat)>=0
        Beta=acos(cos(AOP)/cos(Lat));
    else
        Beta=2*pi-acos(cos(AOP)/cos(Lat));
end % Beta (Auxiliar angle)

Px=R0*cos(Lat)*cos(Beta+RAAN);
Py=R0*cos(Lat)*sin(Beta+RAAN);
Pz=R0*sin(Lat);
P=[Px Py Pz]; P=P/norm(P); %Unit vector

% Q vector
Lat=asin(sin(i)*sin(AOP+pi/2)); % Latitude
if sin(Lat)>=0
        Beta=acos(cos(AOP+pi/2)/cos(Lat));
    else
        Beta=2*pi-acos(cos(AOP+pi/2)/cos(Lat));
end % Beta (Auxiliar angle)

R90=(a*(1-e^2))/(1+e*cos(f0));

Qx=R90*cos(Lat)*cos(Beta+RAAN);
Qy=R90*cos(Lat)*sin(Beta+RAAN);
Qz=R90*sin(Lat);
Q=[Qx Qy Qz]; Q=Q/norm(Q); % Unit vector

r0=R0*cos(f0)*P+R0*sin(f0)*Q; % Initial position (vector)
v0=sqrt(mu_Sun/(a*(1-e^2)))*(-sin(f0)*P+(e+cos(f0))*Q); % Initial velocity (vector)

end