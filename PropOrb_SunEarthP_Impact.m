function [r_Earth,v_Earth,r_Obj,v_Obj]=PropOrb_SunEarthP_Impact(r0_Earth,v0_Earth,r0_Obj,v0_Obj,t) % In m, m/s and s

% t is a vector of time from propagation start time until the propagation
% end time or until the propagated object impacts the Earth

opts=odeset('RelTol',5e-10,'AbsTol',5e-10,'Refine',50,'Events',@(~,r)Event_Impact2(t,r));

[~,rv]=ode45(@rdot2,t,[r0_Earth v0_Earth r0_Obj v0_Obj],opts);

r_Earth=rv(:,1:3);
v_Earth=rv(:,4:6);
r_Obj=rv(:,7:9);
v_Obj=rv(:,10:12);

end