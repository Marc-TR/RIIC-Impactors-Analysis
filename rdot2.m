%% Two body problem with Earth perturbation

function rdot2=rdot2(~,r)

mu_Sun=1.32712428e20; % m^3/s^2
mu_Earth=3.986004415e14; % m^3/s^2

%% Earth
rdot2(1)=r(4);
rdot2(2)=r(5);
rdot2(3)=r(6);

rdot2(4)=-(mu_Sun/norm(r(1:3))^3)*r(1);
rdot2(5)=-(mu_Sun/norm(r(1:3))^3)*r(2);
rdot2(6)=-(mu_Sun/norm(r(1:3))^3)*r(3);

%% Object
rdot2(7)=r(10);
rdot2(8)=r(11);
rdot2(9)=r(12);

% Position vector from Earth to Object
rEO(1)=r(7)-r(1);
rEO(2)=r(8)-r(2);
rEO(3)=r(9)-r(3);

rdot2(10)=-(mu_Sun/norm(r(7:9))^3)*r(7) - mu_Earth*( (rEO(1)/norm(rEO)^3) - (-r(1)/norm(r(1:3))^3) );
rdot2(11)=-(mu_Sun/norm(r(7:9))^3)*r(8) - mu_Earth*( (rEO(2)/norm(rEO)^3) - (-r(2)/norm(r(1:3))^3) );
rdot2(12)=-(mu_Sun/norm(r(7:9))^3)*r(9) - mu_Earth*( (rEO(3)/norm(rEO)^3) - (-r(3)/norm(r(1:3))^3) );

rdot2=rdot2'; % Vector in column form

