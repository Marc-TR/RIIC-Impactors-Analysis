%% Velocity and apparent motion from L1 of 17.518 objects that impact the Earth with Earth Perturbation
% The propagation stops when the object impacts the Earth

clear all; clc; close all;
%-----------------------------------------------------------------------------------------------------------------


%% Constants
AU=149597870700; %m
a_Earth=AU; %m
e_Earth=0;
mu_Sun = 1.32712428e20; %m^3/s^2
i_Earth=0*(pi/180); %Earth inclination
RAAN_Earth=0*(pi/180); %Right ascencion of the ascending node (Earth)
AOP_Earth=0*(pi/180); %Argument of periapsis (Earth)
f0_Earth=0; % Initial true anomaly of Earth
T_Earth=2*pi*sqrt(a_Earth^3/mu_Sun); %Earth period
N=2000; %Step number
ON=8759; % Number of objects set Out (1)
ON2=17518; % Number of objects set In (2)
load('AstList_newDef.mat'); % Load Asteroids
T=4 *24*3600; % Duration of the propagation [s]
t=linspace(0,T,N); % Forward propagation
t_b=linspace(T,0,N); % Backwards propagation
Lim_V=17; % Limiting (apparent) magnitude of our telescope

SUNL1EARTH=[AU 0 0]; % Sun-Earth vector at t=0

%% Orbit backwards propagation (Earth initial conditions)
[r0_Earth,v0_Earth]=IC_from_KEP(a_Earth, e_Earth, i_Earth, RAAN_Earth, AOP_Earth, f0_Earth); % Earth Initial State Vector
[r_Earth,v_Earth]=PropOrb_Sun(r0_Earth,v0_Earth,t_b); % Earth backwards propagation

r0_Earth=r_Earth(N,:); % Earth initial position
v0_Earth=v_Earth(N,:); % Earth initial velocity


%% Prealocated values

VP_Obj=zeros(N,ON2); % Velocity (scalar) of the object relative to the Sun WITH Earth pert. [m/s]
VP_ObjEarth=zeros(N,ON2); % Velocity (scalar) of the object relative to the Earth WITH Earth pert. [m/s]
AppMot_P=zeros(N,ON2); % Apparent motion of the object as seen from L1 WITH Earth pert. [deg/day]
V=zeros(N,ON2); % Apparent magnitude
r_Det=zeros(3,ON2); % Object position at detection (Non inertial reference frame)
t_Obs=zeros(1,ON2); % Time from detection to impact


%-------------------------------------------------------------------------------------


%% Object set 1 (Out)

for index=1:ON
    
    a_Obj=ImpactingObjects{index}.astUpOut(1)*AU; %Semi-major axis [m]
    e_Obj=ImpactingObjects{index}.astUpOut(2); %Asteroid's eccentricity
    i_Obj=ImpactingObjects{index}.astUpOut(3); %Inclination [Rad]
    RAAN_Obj=ImpactingObjects{index}.astUpOut(4); % RAAN [Rad]
    AOP_Obj=ImpactingObjects{index}.astUpOut(5); %Argument of Periapsis [Rad]
    f0_Obj=ImpactingObjects{index}.astUpOut(6); % Asteroid's True anomaly (initial) [Rads]
        
    % Orbit backwards propagation (Object)
    [r0_Obj,v0_Obj]=IC_from_KEP(a_Obj, e_Obj, i_Obj, RAAN_Obj, AOP_Obj, f0_Obj); % Initial state vector
    
    [r_Obj,v_Obj]=PropOrb_Sun(r0_Obj,v0_Obj,t_b); % Backwards propagation
    r0_Obj=r_Obj(N,:); % Object initial position
    v0_Obj=v_Obj(N,:); % Object initial velocity
    
    % Earth and object propagation
    [r_Earth,v_Earth,r_Obj,v_Obj]=PropOrb_SunEarthP_Impact(r0_Earth,v0_Earth,r0_Obj,v0_Obj,t); % Propagation with Earth pert
    
    N2=length(r_Earth(:,1));
    
    % L1
    r_L1=zeros(N2,3);
    v_L1=zeros(N2,3);
    
    for i=1:N2
        r_L1(i,:)=r_Earth(i,:)-((r_Earth(i,:)/norm(r_Earth(i,:)))*1.5*10^9);
        v_L1(i,:)=norm(r_L1(i,:))*norm(v_Earth(i,:))/norm(r_Earth(i,:))*(v_Earth(i,:)/norm(v_Earth(i,:)));
    end
    
    % Vector position from Earth to the object
    r_EO=r_Obj-r_Earth; 
    
    % Velocities
    vP_ObjEarth=v_Obj-v_Earth; % Vector
    
    % Position and velocity vectors from L1 to object
    v_L1O=v_Obj-v_L1; % Velocity vector of the object relative to L1
    r_L1O=r_Obj-r_L1; % Position vector from L1 to object
    
    for i=1:N2
        VP_Obj(i,index)=norm(v_Obj(i,:)); % Object velocity (relative to the Sun) [m/s]
        VP_ObjEarth(i,index)=norm(vP_ObjEarth(i,:)); % Object velocity (relative to the Earth)[m/s]
        AppMot_P(i,index)=norm(cross(r_L1O(i,:),v_L1O(i,:))/(norm(r_L1O(i,:)))^2)*4950355.35; % Apparent motion [deg/day]
        [V(i,index)]=Apparent_Magnitude(r_Obj(i,:),r_L1O(i,:));
        
        if t_Obs(index)==0
            
            Epsilon=acos(dot(r_L1O(i,:),(-r_L1(i,:)))/(norm(r_L1O(i,:))*norm(r_L1(i,:))));
            
            if V(i,index)<=Lim_V && Epsilon>40*(pi/180) % Limiting magnitude and Sun exclusion angle
                t_Obs(index)=t(N2)-t(i); % Time since the object is observable until impact
                
                % Save the position in a non inertial reference frame
                Beta=acos(dot(r_Earth(i,:),SUNL1EARTH)/(norm(r_Earth(i,:))*norm(SUNL1EARTH))); % Earth angular displacement
                r_Det(1,index)=r_EO(i,1)*cos(Beta)-r_EO(i,2)*sin(Beta);
                r_Det(2,index)=r_EO(i,1)*sin(Beta)+r_EO(i,2)*cos(Beta);
                r_Det(3,index)=r_EO(i,3);
                
            end
        end
        
    end
      
    disp(index) % To check the progress
    
end


%% Object set 2 (In)

index=1;

for index2=(ON+1):ON2
    
    a_Obj=ImpactingObjects{index}.astUpIn(1)*AU; %Semi-major axis [m]
    e_Obj=ImpactingObjects{index}.astUpIn(2); %Asteroid's eccentricity
    i_Obj=ImpactingObjects{index}.astUpIn(3); %Inclination [Rad]
    RAAN_Obj=ImpactingObjects{index}.astUpIn(4); % RAAN [Rad]
    AOP_Obj=ImpactingObjects{index}.astUpIn(5); %Argument of Periapsis [Rad]
    f0_Obj=ImpactingObjects{index}.astUpIn(6); % Asteroid's True anomaly (initial) [Rads]
      
    % Orbit backwards propagation (Object)
    [r0_Obj,v0_Obj]=IC_from_KEP(a_Obj, e_Obj, i_Obj, RAAN_Obj, AOP_Obj, f0_Obj); % Initial state vector
       
    [r_Obj,v_Obj]=PropOrb_Sun(r0_Obj,v0_Obj,t_b); % Backwards propagation
    r0_Obj=r_Obj(N,:); % Object initial position
    v0_Obj=v_Obj(N,:); % Object initial velocity
    
    % Earth and object propagation
    [r_Earth,v_Earth,r_Obj,v_Obj]=PropOrb_SunEarthP_Impact(r0_Earth,v0_Earth,r0_Obj,v0_Obj,t); % Propagation with Earth pert
    
    N2=length(r_Earth(:,1));
    
    % L1 Position
    r_L1=zeros(N2,3);
    
    for i=1:N2
        r_L1(i,:)=r_Earth(i,:)-((r_Earth(i,:)/norm(r_Earth(i,:)))*1.5*10^9);
    end
    
    % Vector position from Earth to the object
    r_EO=r_Obj-r_Earth; 
    
    % Velocities
    vP_ObjEarth=v_Obj-v_Earth; % Vector
    
    % Position and velocity vectors from L1 to object
    v_L1O=vP_ObjEarth; % Velocity vector of the object relative to L1 (v_Earth=v_L1)
    r_L1O=r_Obj-r_L1; % Position vector from L1 to object
    
    for i=1:N2
        VP_Obj(i,index2)=norm(v_Obj(i,:)); % Object velocity (relative to the Sun) [m/s]
        VP_ObjEarth(i,index2)=norm(vP_ObjEarth(i,:)); % Object velocity (relative to the Earth)[m/s]
        AppMot_P(i,index2)=norm(cross(r_L1O(i,:),v_L1O(i,:))/(norm(r_L1O(i,:)))^2)*4950355.35; % Apparent motion [deg/day]
        [V(i,index2)]=Apparent_Magnitude(r_Obj(i,:),r_L1O(i,:));
        
        if t_Obs(index2)==0
            
            Epsilon=acos(dot(r_L1O(i,:),(-r_L1(i,:)))/(norm(r_L1O(i,:))*norm(r_L1(i,:))));
                        
            if V(i,index2)<=Lim_V && Epsilon>40*(pi/180) % Limiting magnitude and Sun exclusion angle
                t_Obs(index2)=t(N2)-t(i); % Time since the object is observable until impact
                
                % Save the position in a non inertial reference frame
                Beta=acos(dot(r_Earth(i,:),SUNL1EARTH)/(norm(r_Earth(i,:))*norm(SUNL1EARTH))); % Earth angular displacement
                r_Det(1,index2)=r_EO(i,1)*cos(Beta)-r_EO(i,2)*sin(Beta);
                r_Det(2,index2)=r_EO(i,1)*sin(Beta)+r_EO(i,2)*cos(Beta);
                r_Det(3,index2)=r_EO(i,3);
                
            end
        end
        
    end
    
    
    
    disp(index2) % To check the progress
    index=index+1;
    
end