function [r,v]=PropOrb_Sun(r0,v0,t) % In m, m/s and s

% t is a vector of time from propagation start time until the propagation
% end time

opts=odeset('RelTol',5e-10,'AbsTol',5e-10);

[~,rv]=ode45(@rdot,t,[r0 v0],opts);

r=rv(:,1:3);
v=rv(:,4:6);

end