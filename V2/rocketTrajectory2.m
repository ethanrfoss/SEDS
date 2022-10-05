function rocketTrajectory2(R,Reco,H,LaunchLat,LaunchLong,LaunchAlt,compressibility,dt)

%% Environemntal Parameters:
global g
g = 9.81; % Gravitational Acceleration(m/s^2)

%% Rocket Initial Conditions:

% Iterator:
n = 1;

% Time:
t(n) = 0; %(sec)

% System State, Initial Conditions:
systemState(:,n) = [0;... % X - Easting(m)
                   0;... % Y = Northing(m)
                   LaunchAlt;... % Z - Altitude(m)
                   0;... % Vx - Velocity in X-dir(m/s)
                   0;... % Vy - Velocity in Y-dir(m/s)
                   0;... % Vz - Velocity in Z-dir(m/s)
                   0;... % alpha - 1-axis euler-angle(rad)
                   0;... % beta - 2-axis euler-angle(rad)
                   0;... % gamma - 3-axis euler-angle(rad)
                   0;... % w1 - Rotation Rate about rocket 1-axis(rad/s)
                   0;... % w2 - Rotation Rate about rocket 2-axis(rad/s)
                   0;... % w3 - Rotation Rate about rocket 3-axis(rad/s)
                   0;... % Xd - Easting of Drogue Chute(m)
                   0;... % Yd - Northing of Drogue Chute(m)
                   0;... % Zd - Altitude of Drogue Chute(m)
                   0;... % Vxd - Velocity of Drogue Chute in X-dir(m/s)
                   0;... % Vyd - Velocity of Drogue Chute in Y-dir(m/s)
                   0;... % Vzd - Velocity of Drogue Chute in Z-dir(m/s)
                   0;... % Xm - Easting of Main Chute(m)
                   0;... % Ym - Northing of Main Chute(m)
                   0;... % Zm - Altitude of Main Chute(m)
                   0;... % Vxm - Velocity of Main Chute in X-dir(m/s)
                   0;... % Vym - Velocity of Main Chute in Y-dir(m/s)
                   0]; % Vzm - Velocity of Main Chute in Z-dir(m/s)

% Forces:
FR(:,n) = [0;... % Fx - Force in X-dir(N)
         0;... % Fy - Force in Y-dir(N)
         0]; % Fz - Force in Z-dir(N)
Fd(n) = 0; % Drag Force(N)
Fn(n) = 0; % Normal Force(N)

% Moments:
MR(:,n) = [0;... % Mx - Moment in X-dir(Nm)
         0;... % My - Moment in Y-dir(Nm)
         0]; % Mz - Moment in Z-dir(Nm)

%% Other Global Variables:
global LaunchLat LaunchLong LaunchAlt H compressibility;

global drogueDeployed mainDeployed;
drogueDeployed = false;
mainDeployed = false;

%% Launch Rail Loops:

%% Ascent:
while systemState(6,n)>=0 || R.thrust(n-1) ~= 0
   
    % Check if Rocket Bounds Exceeded:
    if n > length(R.mass)
        R = extendBounds(R,n);
    end
    
    % Iterate time:
    t(n+1) = t(n) + dt; 
    
    % Apply Dynamics using RK4:
    systemState(:,n+1) = RK4(@systemN,systemState(:,n),@systemFandM,R,Reco,n,dt);
    
    % Iterator:
    n = n+1;
end

%% Drogue Chute Initial Conditions:

% Iterator:
nDrogueDeploy = n;

% Drogue State:
systemState(13:18,n) = systemState(1:6,n);

% Indicate Deployment:
drogueDeployed = true;

%% Drogue Chute Loop:

while systemState(3,n) >= Reco.MainDeployHeight+LaunchAlt
    
    % Check if Rocket Bounds Exceeded:
    if n > length(R.mass)
        R = extendBounds(R,n);
    end
    
    % Iterate time:
    t(n+1) = t(n) + dt; 
    
    % Apply Dynamics using RK4:
    systemState(:,n+1) = RK4(@systemN,systemState(:,n),@systemFandM,R,Reco,n,dt);
    
    % Iterator:
    n = n+1;
    
end

%% Main Chute Initial Conditions:

% Iterator:
nMainDeploy = n;

% Drogue State:
systemState(19:24,n) = systemState(1:6,n);

% Indicate Deployment:
mainDeployed = true;
%% Main Chute Loop:

while systemState(3,n) >= LaunchAlt
    
     % Check if Rocket Bounds Exceeded:
    if n > length(R.mass)
        R = extendBounds(R,n);
    end
    
    % Iterate time:
    t(n+1) = t(n) + dt; 
    
    % Apply Dynamics using RK4:
    systemState(:,n+1) = RK4(@systemN,systemState(:,n),@systemFandM,R,Reco,n,dt);
    
    % Iterator:
    n = n+1;
        
end

plot3(systemState(1,:),systemState(2,:),systemState(3,:));
end

%% Forces on Rocket System, Including Parachutes, from Current State:
function [FR,MR,FD,FM] = systemFandM(systemState,R,Reco,n)

% Load Global Variables:
global LaunchLat LaunchLong LaunchAlt H compressibility g drogueDeployed mainDeployed;

% Determine Atmospheric Conditions:
[Temp,~,rho] = standardAtmosphere(systemState(3));

% Get Wind Speed:
[wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,systemState(3));

% Determine Velocities, Reynolds number, Mach number:
VinAir = [systemState(4)-wlong;systemState(5)-wlat;systemState(6)];
VmaginAir = sqrt(sum(VinAir.^2));
Re = rho*VmaginAir*R.lTR/sutherland(Temp);
Ma = VmaginAir/(sqrt(1.4*287*Temp));

% Angle of Attack:
if systemState(3) < LaunchAlt + H
    a = 0;    
else
    a = real(acos(dot(VinAir/VmaginAir,Tbi(systemState(7),systemState(8),systemState(9))*[0;0;1])));    
end

% Drag Coefficients:
[CN, CA, Xcp] = dragCoefs(R,Re,a,Ma,compressibility);

% Determine Forces:
Fd = .5*CA*rho*R.Aref*VmaginAir^2;
Fn = .5*CN*rho*R.Aref*VmaginAir^2;
if a == 0
    % If AoA is zero, no normal force
    dirFn = [0;0;0];
else
    % Determine direction of Force(cross of moment direction with b3 axis)
    dirFn = cross(cross(Tib(systemState(7),systemState(8),systemState(9))*VinAir/VmaginAir,[0;0;1]),[0;0;1]);
    % Normalize Force Direction
    dirFn = dirFn/sqrt(sum(dirFn.^2));
end
FR = Tbi(systemState(7),systemState(8),systemState(9))*([dirFn(1)*Fn;dirFn(2)*Fn;R.thrust(n)-Fd])-[0;0;g*R.mass(n)];

% Determine Moments:
if a == 0 
    % If angle of attack is zero, no moment(avoids singularity)
    dirM = [0;0;0];
else
    % Calculate the direction of the Moment(cross of b3 axis with Velocity vector)
    dirM = cross([0;0;1],Tib(systemState(7),systemState(8),systemState(9))*VinAir/VmaginAir); % Moment Due to Normal Force
    % Normalize Moment Direction
    dirM = dirM/sqrt(sum(dirM.^2));
end
if R.thrust(n) == 0
    % If thrust is zero, no thrust damping(avoids singularity)
    Td = [0;0;0];
else
    % Calculate vector of thrust damping(in same direction as rotational velocity)
    Td = -R.massflow(n)*abs(R.rocketCentroid{n}(3)^2-((R.LOX.cent{n}(3)*R.LOX.mass(n)+R.CH4.cent{n}(3)*R.CH4.mass(n))/(R.LOX.mass(n)+R.CH4.mass(n)))^2)*systemState(10:12); % Torque from Torque Damping
end
MR = [Fn*(R.rocketCentroid{n}(3) + Xcp)*dirM(1);Fn*(R.rocketCentroid{n}(3) + Xcp)*dirM(2);0] + Td;

% Drogue:
if drogueDeployed
    
    % Rope Parameters:
    RopeStart = Tbi(systemState(7),systemState(8),systemState(9))*[0;0;1]*(-R.ln-R.rocketCentroid{n}(3))+systemState(1:3);
    Ldrogue = sqrt(sum((systemState(13:15)-systemState(1:3)).^2));
    Vldrogue = sum((systemState(1:3)-systemState(13:15)).*(systemState(4:6)-systemState(16:18)))/Ldrogue;
   
    % Drogue Velocities:
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(systemState(15));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,systemState(15)); %wind
    
    % Determine Velocities:
    VdrogueinAir = systemState(16:18) - [wlong;wlat;0];
    VdroguemaginAir = sqrt(sum(VdrogueinAir.^2));
    
    % Rope Tension: 
    if Ldrogue > Reco.Drogue.L
        Tdrogue = Reco.Drogue.k*(Ldrogue-Reco.Drogue.L)+Vldrogue*Reco.Drogue.c;
    else
        Tdrogue = 0;
    end
    
    % Rope Forces and Moments:
    if Tdrogue == 0
        dirFt = [0;0;0];
        Ft = [0;0;0];
        Mt = [0;0;0];
    else
        dirFt = (systemState(13:15)-RopeStart)/Ldrogue;
        Ft = dirFt*Tdrogue;
        Mt = Tdrogue*cross(RopeStart-systemState(1:3),dirFt);
    end
    
    % Drogue Forces:
    Fpdrogue = .5*Reco.Drogue.Cd*Reco.Drogue.A*rho*VdroguemaginAir^2;
    FD = -Fpdrogue*VdrogueinAir/VdroguemaginAir + [0;0;-Reco.Drogue.m*g] - Tdrogue*dirFt;
    
    % Drogue Rope Forces on Rocket:
    FR = FR + Ft;
    MR = MR + Mt;
else
   
    FD = [0,0,0];
    
end

% Main:
if mainDeployed
    
    % Rope Parameters:
    RopeStart = Tbi(systemState(7),systemState(8),systemState(9))*[0;0;1]*(-R.ln-R.rocketCentroid{n}(3))+systemState(1:3);
    Lmain = sqrt(sum((systemState(19:21)-systemState(1:3)).^2));
    Vlmain = sum((systemState(1:3)-systemState(19:21)).*(systemState(4:6)-systemState(22:24)))/Lmain;
    
    % Drogue Velocities:
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(systemState(15));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,systemState(15)); %wind
    
    % Determine Velocities:
    VmaininAir = systemState(22:24) - [wlong;wlat;0];
    VmainmaginAir = sqrt(sum(VmaininAir.^2));
    
    % Rope Tension: 
    if Lmain > Reco.Main.L
        Tmain = Reco.Main.k*(Lmain-Reco.Main.L)+Vlmain*Reco.Main.c;
    else
        Tmain = 0;
    end
    
    % Rope Forces and Moments:
    if Tmain == 0
        dirFt = [0;0;0];
        Ft = [0;0;0];
        Mt = [0;0;0];
    else
        dirFt = (systemState(19:21)-RopeStart)/Lmain;
        Ft = dirFt*Tmain;
        Mt = Tmain*cross(RopeStart-systemState(1:3),dirFt);
    end
    
    % Drogue Forces:
    Fpmain = .5*Reco.Main.Cd*Reco.Main.A*rho*VmainmaginAir^2;
    FM = -Fpmain*VmaininAir/VmainmaginAir + [0;0;-Reco.Main.m*g] - Tmain*dirFt;
    
    % Drogue Rope Forces on Rocket:
    FR = FR + Ft;
    MR = MR + Mt;
else
    
    FM = [0,0,0];
    
end

end

%% Derivative of State Vector, Function of Current State:
function N = systemN(systemState,FandM,R,Reco,n)

[FR,MR,FD,FM] = FandM(systemState,R,Reco,n);

N = [systemState(4);...
    systemState(5);...
    systemState(6);...
    FR(1)/R.mass(n);...%+rocketState(12)*rocketState(5)-rocketState(11)*rocketState(6);...
    FR(2)/R.mass(n);...%+rocketState(10)*rocketState(6)-rocketState(12)*rocketState(4);...
    FR(3)/R.mass(n);...%+rocketState(11)*rocketState(4)-rocketState(10)*rocketState(5);...
    sin(systemState(9))/cos(systemState(8))*systemState(11)+cos(systemState(9))/cos(systemState(8))*systemState(12);...
    cos(systemState(9))*systemState(11)-sin(systemState(9))*systemState(12);...
    systemState(10)+sin(systemState(9))*tan(systemState(8))*systemState(11)+cos(systemState(9))*tan(systemState(8))*systemState(12);...
    (MR(1)+(R.rocketInertia{n}(2,2)-R.rocketInertia{n}(3,3))*systemState(11)*systemState(12))/R.rocketInertia{n}(1,1);...
    (MR(2)+(R.rocketInertia{n}(3,3)-R.rocketInertia{n}(1,1))*systemState(10)*systemState(12))/R.rocketInertia{n}(2,2);...
    (MR(3)+(R.rocketInertia{n}(1,1)-R.rocketInertia{n}(2,2))*systemState(10)*systemState(11))/R.rocketInertia{n}(3,3);...
    systemState(16);...
    systemState(17);...
    systemState(18);...
    FD(1)/Reco.Drogue.m;...
    FD(2)/Reco.Drogue.m;...
    FD(3)/Reco.Drogue.m;...
    systemState(22);...
    systemState(23);...
    systemState(24);...
    FM(1)/Reco.Main.m;...
    FM(2)/Reco.Main.m;...
    FM(3)/Reco.Main.m];
    
end

%% RK4 Method to Integrate to Next State:
function nextState = RK4(computeN,currentState,FandM,R,Reco,n,dt)

f1 = computeN(currentState,FandM,R,Reco,n);
f2 = computeN(currentState+dt*f1/2,FandM,R,Reco,n);
f3 = computeN(currentState+dt*f2/2,FandM,R,Reco,n);
f4 = computeN(currentState+dt*f3,FandM,R,Reco,n);

nextState = currentState+dt*(f1/6+(f2+f3)/3+f4/6);

end

%% Extend Bounds of Rocket Variables:
function R = extendBounds(R,n)

R.thrust(n) = 0;
R.mass(n) = R.mass(n-1);
R.massflow(n) = 0;
R.LOX.mass(n) = 0;
R.CH4.mass(n) = 0;
R.LOX.cent{n} = R.LOX.cent{n-1};
R.CH4.cent{n} = R.CH4.cent{n-1};
R.LOX.I{n} = R.LOX.I{n-1};
R.CH4.I{n} = R.CH4.I{n-1};
R.rocketCentroid{n} = R.rocketCentroid{n-1};
R.rocketInertia{n} = R.rocketInertia{n-1};
        
end

%% Inertial to Body Frame Transformation:
function T = Tib(a,b,g)

T = [cos(a)*cos(b),cos(b)*sin(a),-sin(b); cos(a)*sin(b)*sin(g)-cos(g)*sin(a),cos(a)*cos(g)+sin(a)*sin(b)*sin(g),cos(b)*sin(g); sin(a)*sin(g)+cos(a)*cos(g)*sin(b),cos(g)*sin(a)*sin(b)-cos(a)*sin(g),cos(b)*cos(g)]; 

end

%% Body to Inertial Frame Transformation:
function T = Tbi(a,b,g)

T = Tib(a,b,g)';

end