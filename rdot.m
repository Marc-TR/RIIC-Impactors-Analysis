%% Two Body Problem
% Function to set the differential Equation of Motion of a Satellite

function rdot=rdot(~,r)

mu_Sun = 1.32712428e20; % m^3/s^2

rdot(1)=r(4);
rdot(2)=r(5);
rdot(3)=r(6);

rdot(4)=-(mu_Sun/norm(r(1:3))^3)*r(1);
rdot(5)=-(mu_Sun/norm(r(1:3))^3)*r(2);
rdot(6)=-(mu_Sun/norm(r(1:3))^3)*r(3);

rdot=rdot'; % Vector in column form

