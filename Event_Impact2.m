function [value,isterminal,direction] = Event_Impact2(~,r)

R_Earth=6378136.3; % Earth radius [m]

R_EO=norm(r(7:9)-r(1:3)); % Distance between Earth and object

if R_EO-R_Earth<=0
    value=0; % When the object reaches the Earth surface (impact) it will be zero
else
    value=1;
end

isterminal=1;
direction=0;

end