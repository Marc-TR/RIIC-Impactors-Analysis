clc; close all; clear all;

%% Constants
AU=149597870700; %m
R_Earth=6378136.3; %m
a_Earth=AU; %m
e_Earth=0;
mu_Sun = 1.32712428e20; %m^3/s^2
R_Sun=695508000; %m
i_Earth=0*(pi/180); %Earth inclination
RAAN_Earth=0*(pi/180); %Right ascencion of the ascending node (Earth)
AOP_Earth=0*(pi/180); %Argument of periapsis (Earth)
f0_Earth=0;
T_Earth=2*pi*sqrt(a_Earth^3/mu_Sun); %Earth period
N=500; %Step number

T=T_Earth; % Seconds until impact
t=linspace(0,T,N);
t_b=linspace(T,0,N); % For backwards propagation

%% Orbit Propagation (Earth)
[r0_Earth,v0_Earth]=IC_from_KEP(a_Earth, e_Earth, i_Earth, RAAN_Earth, AOP_Earth, f0_Earth); % Earth Initial State Vector
[r_Earth,v_Earth]=PropOrb_Sun(r0_Earth,v0_Earth,t_b); % Earth backwards propagation

r0_Earth=r_Earth(N,:);
v0_Earth=v_Earth(N,:);

%% Asteroid
load('AstList_newDef.mat');

index=630; % Number of object

a_Obj=ImpactingObjects{index}.astUpOut(1)*AU; %Semi-major axis [m]
e_Obj=ImpactingObjects{index}.astUpOut(2); %Asteroid's eccentricity
i_Obj=ImpactingObjects{index}.astUpOut(3); %Inclination [Rad]
RAAN_Obj=ImpactingObjects{index}.astUpOut(4); % RAAN [Rad]
AOP_Obj=ImpactingObjects{index}.astUpOut(5); %Argument of Periapsis [Rad]
T_Obj=2*pi*sqrt(a_Obj^3/mu_Sun); %Period [s]
f0_Obj=ImpactingObjects{index}.astUpOut(6);

[r0_Obj, v0_Obj]=IC_from_KEP(a_Obj, e_Obj, i_Obj, RAAN_Obj, AOP_Obj, f0_Obj);

[r_Obj, v_Obj]=PropOrb_Sun(r0_Obj, v0_Obj, t_b); % Backwards propagation
r0_Obj=r_Obj(N,:); % We only want the last value
v0_Obj=v_Obj(N,:);% We only want the last value

% Earth and Object propagation
[r_Earth,v_Earth,r_Obj, v_Obj]=PropOrb_SunEarthP_Impact(r0_Earth, v0_Earth, r0_Obj, v0_Obj, t);

N2=r_Earth';
N2=length(N2(1,:));

%% L1 Position
r_L1=zeros(N2,3);
v_L1=zeros(N2,3);

for i=1:N2
    r_L1(i,:)=r_Earth(i,:)-((r_Earth(i,:)/norm(r_Earth(i,:)))*1.5*10^9);
    v_L1(i,:)=norm(r_L1(i,:))*norm(v_Earth(i,:))/norm(r_Earth(i,:))*(v_Earth(i,:)/norm(v_Earth(i,:)));
end

%% Relative position (Earth Centered)
r_EO=r_Obj-r_Earth; % Position vector from earth (0,0,0) to Asteroid
r_EL1=r_L1-r_Earth; % Position vector from earth (0,0,0) to L1

%% Apparent Magnitude
V=zeros(1,length(r_Earth(:,1)));
r_L1O=r_Obj-r_L1; % Position vector from L1 to object

for i=1:length(r_Earth(:,1))
   [V(i)]=Apparent_Magnitude(r_Obj(i,:),r_L1O(i,:));
end

% {
%% Plot Orbits

figure(1)
grid on;
hold on;
daspect([1 1 1])
axis equal
axis([-1.1*AU 1.1*AU -1.1*AU 1.1*AU -0.2*AU 0.2*AU])
view(3)

% Plot sphere (Sun)
[A,B,C]=ellipsoid(0,0,0,R_Sun,R_Sun,R_Sun);
surf (A,B,C)

% Initiliaze animated lines with color and line width properties
Earth_an=animatedline('Color','b','LineWidth',0.2);
Obj_an=animatedline('Color','k','LineWidth',0.1);
L1_an=animatedline('Color','g','LineWidth',0.2);

% Earth trajectory
addpoints(Earth_an,r_Earth(1,1),r_Earth(1,2),r_Earth(1,3)); % Earth trajectory
Earth=scatter3(r_Earth(1,1),r_Earth(1,2),r_Earth(1,3), 15, 'b','filled');

% Object trajectory
addpoints(Obj_an,r_Obj(1,1),r_Obj(1,2),r_Obj(1,3));
OBJ=scatter3(r_Obj(1,1),r_Obj(1,2),r_Obj(1,3), 8, 'k','filled');

% L1 "point" trajectory
addpoints(L1_an,r_L1(1,1),r_L1(1,2),r_L1(1,3)); 
L1=scatter3(r_L1(1,1),r_L1(1,2),r_L1(1,3), 6, 'g','filled'); % Sun centered

legend([Earth OBJ L1],'Earth','Asteroid', 'L1')
%Maximise
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

figure(2)
grid on;
hold on;
daspect([1 1 1])
axis equal

% Initialize line and scatter
EO_an=animatedline('Color','k','LineWidth',0.2);
EO=scatter3(r_EO(1,1),r_EO(1,2),r_EO(1,3), 6, 'k','filled'); % Earth centered
EL1_an=animatedline('Color','g','LineWidth',0.05);
EL1=scatter3(r_EL1(1,1),r_EL1(1,2),r_EL1(1,3), 6, 'g','filled'); % Earth centered

axis([-10*R_Earth 10*R_Earth -10*R_Earth 10*R_Earth -2*R_Earth 2*R_Earth])
view(3)

%Plot sphere (Earth)
[X,Y,Z]=ellipsoid(0,0,0,R_Earth,R_Earth,R_Earth);
surf (X,Y,Z)

%Maximise
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);



for i=1:N2   %1:N
    
    % Earth
    addpoints(Earth_an,r_Earth(i,1),r_Earth(i,2),r_Earth(i,3));
    set(Earth,'xdata',r_Earth(i,1),'ydata',r_Earth(i,2),'zdata',r_Earth(i,3));
    
    % Obj
    addpoints(Obj_an,r_Obj(i,1),r_Obj(i,2),r_Obj(i,3));
    set(OBJ,'xdata',r_Obj(i,1),'ydata',r_Obj(i,2),'zdata',r_Obj(i,3));
    
    % L1
    addpoints(L1_an,r_L1(i,1),r_L1(i,2),r_L1(i,3));
    set(L1,'xdata',r_L1(i,1),'ydata',r_L1(i,2),'zdata',r_L1(i,3));
    
    % Object (Earth centered)
    addpoints(EO_an,r_EO(i,1),r_EO(i,2),r_EO(i,3));
    set(EO,'xdata',r_EO(i,1),'ydata',r_EO(i,2),'zdata',r_EO(i,3));
    
    % L1 (Earth centered)
    addpoints(EL1_an,r_EL1(i,1),r_EL1(i,2),r_EL1(i,3));
    set(EL1,'xdata',r_EL1(i,1),'ydata',r_EL1(i,2),'zdata',r_EL1(i,3));
       
    drawnow %Load the updates to the plot
    disp(i)
    
end

%}




