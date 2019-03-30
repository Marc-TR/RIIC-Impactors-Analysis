function [V]=Apparent_Magnitude(r_Obj,r_L1O)

AU=149597870700; %m
Diam=0.07; % Diameter [km]
Albedo=0.154;

kappa=acos(dot(r_Obj,r_L1O)/(norm(r_Obj)*norm(r_L1O)));

R_Obj=norm(r_Obj)/AU; % Distance from the object to the Sun [AU]
R_L1O=norm(r_L1O)/AU; % Distance from the object to the L1 [AU]

H=5*(log10(1329)-0.5*log10(Albedo)-log10(Diam));
phi1=exp(-3.33*tan(kappa/2)^0.63);
phi2=exp(-1.87*tan(kappa/2)^1.22);
G=0.15;

V=H+5*log10(R_L1O*R_Obj)-2.5*log10((1-G)*phi1+G*phi2);

end