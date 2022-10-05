function [t,R,Reco,rocketState,drogueState,mainState,a,Ma,FR,nDrogueDeploy,Ldrogue,Tdrogue,Vldrogue,nMainDeploy,Lmain,Tmain,Vlmain,RopeStart] = rocketTrajectory(R,Reco,H,LaunchLat,LaunchLong,LaunchAlt,compressibility,dt)

%% Environemntal Parameters:
g = 9.81; % Gravitational Acceleration(m/s^2)

%% Transformation matrices:
% Inertial to Body Frame:
Tib = @(a,b,g) [cos(a)*cos(b),cos(b)*sin(a),-sin(b); cos(a)*sin(b)*sin(g)-cos(g)*sin(a),cos(a)*cos(g)+sin(a)*sin(b)*sin(g),cos(b)*sin(g); sin(a)*sin(g)+cos(a)*cos(g)*sin(b),cos(g)*sin(a)*sin(b)-cos(a)*sin(g),cos(b)*cos(g)]; 
% Body to Inertial Frame:
Tbi = @(a,b,g) Tib(a,b,g)';

%% Rocket Initial Conditions:

% Iterator:
n = 1;

% Time:
t(n) = 0; %(sec)

% Rocket State:
rocketState(:,n) = [0;... % X - Easting(m)
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
                   0]; % w3 - Rotation Rate about rocket 3-axis(rad/s)

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


%% Launch Rail Loop:

while rocketState(3,n) <= LaunchAlt+H
    
    % Check if Rocket Bounds Exceeded:
    if n > length(R.mass)
        R = extendBounds(R,n);
    end
    
    % Iterate time:
    t(n+1) = t(n) + dt;  
    
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(rocketState(3,n));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,rocketState(3,n));

    % Determine Velocities, Reynolds number, Mach number:
    VinAir(:,n) = [rocketState(4,n)-wlong;rocketState(5,n)-wlat;rocketState(6,n)];
    VmaginAir = sqrt(sum(VinAir(:,n).^2));
    Re(n) = rho*VmaginAir*R.lTR/sutherland(Temp);
    Ma(n) = VmaginAir/(sqrt(1.4*287*Temp));
    
    % Assume Zero Angle of Attack on Launch Rail:
    a(n) = 0;    
    
    % Drag Coefficients:
    [~, CA, Xcp(n)] = dragCoefs(R,Re(n),a(n),Ma(n),compressibility);

    % Determine Forces:
    Fd(n) = .5*CA*rho*R.Aref*VmaginAir^2;
    Fn(n) = 0;
    FR(:,n) = [0;0;R.thrust(n)-Fd(n)-R.mass(n)*g];
    
    % Determine Moments:
    MR(:,n) = [0;0;0]; 
    
    % Apply Dynamics using RK4:
    rocketState(:,n+1) = RK4(@rocketN,rocketState(:,n),FR(:,n),MR(:,n),R,n,dt);
    
    % Iterator:
    n = n+1;
end

%% Ascent Loop

while rocketState(6,n)>=0 || R.thrust(n-1) ~= 0
    
    % Check if Rocket Bounds Exceeded:
    if n > length(R.mass)
        R = extendBounds(R,n);
    end
    
    % Iterate time:
    t(n+1) = t(n) + dt;  
    
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(rocketState(3,n));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,rocketState(3,n));

    % Determine Velocities, Reynolds number, Mach number:
    VinAir(:,n) = [rocketState(4,n)-wlong;rocketState(5,n)-wlat;rocketState(6,n)];
    VmaginAir = sqrt(sum(VinAir(:,n).^2));
    Re(n) = rho*VmaginAir*R.lTR/sutherland(Temp);
    Ma(n) = VmaginAir/(sqrt(1.4*287*Temp));
    
    % Assume Zero Angle of Attack on Launch Rail:
    a(n) = real(acos(dot(VinAir(:,n)/VmaginAir,Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*[0;0;1])));    
    
    % Drag Coefficients:
    [CN, CA, Xcp(n)] = dragCoefs(R,Re(n),a(n),Ma(n),compressibility);

    % Determine Forces:
    Fd(n) = .5*CA*rho*R.Aref*VmaginAir^2;
    Fn(n) = .5*CN*rho*R.Aref*VmaginAir^2;
    if a(n) == 0
        % If AoA is zero, no normal force
        dirFn = [0;0;0];
    else
        % Determine direction of Force(cross of moment direction with b3 axis)
        dirFn = cross(cross(Tib(rocketState(7,n),rocketState(8,n),rocketState(9,n))*VinAir(:,n)/VmaginAir,[0;0;1]),[0;0;1]);
        % Normalize Force Direction
        dirFn = dirFn/sqrt(sum(dirFn.^2));
    end
    FR(:,n) = Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*([dirFn(1)*Fn(n);dirFn(2)*Fn(n);R.thrust(n)-Fd(n)])-[0;0;g*R.mass(n)];
    
    % Determine Moments:
    if a(n) == 0 
        % If angle of attack is zero, no moment(avoids singularity)
        dirM = [0;0;0];
    else
        % Calculate the direction of the Moment(cross of b3 axis with Velocity vector)
        dirM = cross([0;0;1],Tib(rocketState(7,n),rocketState(8,n),rocketState(9,n))*VinAir(:,n)/VmaginAir); % Moment Due to Normal Force
        % Normalize Moment Direction
        dirM = dirM/sqrt(sum(dirM.^2));
    end
    if R.thrust(n) == 0
        % If thrust is zero, no thrust damping(avoids singularity)
        Td = [0;0;0];
    else
        % Calculate vector of thrust damping(in same direction as rotational velocity)
        Td = -R.massflow(n)*abs(R.rocketCentroid{n}(3)^2-((R.LOX.cent{n}(3)*R.LOX.mass(n)+R.CH4.cent{n}(3)*R.CH4.mass(n))/(R.LOX.mass(n)+R.CH4.mass(n)))^2)*rocketState(10:12,n); % Torque from Torque Damping
    end
    MR(:,n) = [Fn(n)*(R.rocketCentroid{n}(3) + Xcp(n))*dirM(1);Fn(n)*(R.rocketCentroid{n}(3) + Xcp(n))*dirM(2);0] + Td; 

    % Apply Dynamics using RK4:
    rocketState(:,n+1) = RK4(@rocketN,rocketState(:,n),FR(:,n),MR(:,n),R,n,dt);
    
    % Iterator:
    n = n+1;
    
end

%% Drogue Chute Initial Conditions:

% Iterator:
nDrogueDeploy = n;

% Drogue State:
drogueState(:,n) = rocketState(1:6,n);

%% Drogue Chute Loop:

while rocketState(3,n) >= Reco.MainDeployHeight+LaunchAlt
    
    %Extend Rocket Bounds:
    R = extendBounds(R,n);
    
    % Rope Parameters:
    RopeStart(:,n) = Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*[0;0;1]*(-R.ln-R.rocketCentroid{n-1}(3))+rocketState(1:3,n);
    Ldrogue(n) = sqrt(sum((drogueState(1:3,n)-rocketState(1:3,n)).^2));
    Vldrogue(n) = sum((rocketState(1:3,n)-drogueState(1:3,n)).*(rocketState(4:6,n)-drogueState(4:6,n)))/Ldrogue(n);

    % Iterate time:
    t(n+1) = t(n) + dt;  
    
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(rocketState(3,n));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,rocketState(3,n));

    % Determine Velocities, Reynolds number, Mach number:
    VinAir(:,n) = [rocketState(4,n)-wlong;rocketState(5,n)-wlat;rocketState(6,n)];
    VmaginAir = sqrt(sum(VinAir(:,n).^2));
    Re(n) = rho*VmaginAir*R.lTR/sutherland(Temp);
    Ma(n) = VmaginAir/(sqrt(1.4*287*Temp));
    
    % Assume Zero Angle of Attack on Launch Rail:
    a(n) = real(acos(dot(VinAir(:,n)/VmaginAir,Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*[0;0;1])));    
    
    % Drag Coefficients(IP):
    %[CN, CA, Xcp(n)] = dragCoefs(R,Re(n),a(n),Ma(n),compressibility);
    CN = 0;
    CA = 0;
    
    % Rope Tension: 
    if Ldrogue(n) > Reco.Drogue.L
        Tdrogue(n) = Reco.Drogue.k*(Ldrogue(n)-Reco.Drogue.L)+Vldrogue(n)*Reco.Drogue.c;
    else
        Tdrogue(n) = 0;
    end
    
    % Rope Forces and Moments:
    if Tdrogue(n) == 0
        dirFt = [0;0;0];
        Ft = [0;0;0];
        Mt = [0;0;0];
    else
        dirFt = (drogueState(1:3,n)-RopeStart(:,n))/Ldrogue(n);
        Ft = dirFt*Tdrogue(n);
        Mt = Tdrogue(n)*cross(RopeStart(:,n)-rocketState(1:3,n),dirFt);
    end
    
    % Determine Forces:
    Fd(n) = .5*CA*rho*R.Aref*VmaginAir^2;
    Fn(n) = .5*CN*rho*R.Aref*VmaginAir^2;
    if a(n) == 0
        % If AoA is zero, no normal force
        dirFn = [0;0;0];
    else
        % Determine direction of Force(cross of moment direction with b3 axis)
        dirFn = cross(cross(Tib(rocketState(7,n),rocketState(8,n),rocketState(9,n))*VinAir(:,n)/VmaginAir,[0;0;1]),[0;0;1]);
        % Normalize Force Direction
        dirFn = dirFn/sqrt(sum(dirFn.^2));
    end
    FR(:,n) = Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*([dirFn(1)*Fn(n);dirFn(2)*Fn(n);R.thrust(n)-Fd(n)])-[0;0;g*R.mass(n)]+Ft;
    
    % Determine Moments:
    if a(n) == 0 
        % If angle of attack is zero, no moment(avoids singularity)
        dirM = [0;0;0];
    else
        % Calculate the direction of the Moment(cross of b3 axis with Velocity vector)
        dirM = cross([0;0;1],Tib(rocketState(7,n),rocketState(8,n),rocketState(9,n))*VinAir(:,n)/VmaginAir); % Moment Due to Normal Force
        % Normalize Moment Direction
        dirM = dirM/sqrt(sum(dirM.^2));
    end
    MR(:,n) = [Fn(n)*(R.rocketCentroid{n}(3) + Xcp(end))*dirM(1);Fn(n)*(R.rocketCentroid{n}(3) + Xcp(end))*dirM(2);0] + Mt; 

    % Apply Rocket Dynamics using RK4:
    rocketState(:,n+1) = RK4(@rocketN,rocketState(:,n),FR(:,n),MR(:,n),R,n,dt);
    
    % Parachute Velocities:
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(drogueState(3,n));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,drogueState(3,n)); %wind
    
    % Determine Velocities:
    VdrogueinAir(:,n) = drogueState(4:6,n) - [wlong;wlat;0];
    VdroguemaginAir = sqrt(sum(VdrogueinAir(:,n).^2));
    
    % Parachute Forces:
    Fpdrogue(n) = .5*Reco.Drogue.Cd*Reco.Drogue.A*rho*VdroguemaginAir^2;
    FDR(:,n) = -Fpdrogue(n)*VdrogueinAir(:,n)/VdroguemaginAir + [0;0;-Reco.Drogue.m*g] - Tdrogue(n)*dirFt;
    
    % Apply Drogue Dynamics using RK4:
    drogueState(:,n+1) = RK4(@drogueN,drogueState(:,n),FDR(:,n),1,Reco,n,dt);
    
    % Iterator:
    n = n+1;
    
end

%% Main Chute Initial Conditions:

% Iterator:
nMainDeploy = n;

% Drogue State:
mainState(:,n) = rocketState(1:6,n);

%% Drogue Chute Loop:

while rocketState(3,n) >= LaunchAlt
    
    %Extend Rocket Bounds:
    R = extendBounds(R,n);
    
    % Rope Parameters:
    RopeStart(:,n) = Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*[0;0;1]*(-R.ln-R.rocketCentroid{n-1}(3))+rocketState(1:3,n);
    Ldrogue(n) = sqrt(sum((drogueState(1:3,n)-rocketState(1:3,n)).^2));
    Vldrogue(n) = sum((rocketState(1:3,n)-drogueState(1:3,n)).*(rocketState(4:6,n)-drogueState(4:6,n)))/Ldrogue(n);
    Lmain(n) = sqrt(sum((mainState(1:3,n)-rocketState(1:3,n)).^2));
    Vlmain(n) = sum((rocketState(1:3,n)-mainState(1:3,n)).*(rocketState(4:6,n)-mainState(4:6,n)))/Lmain(n);

    % Iterate time:
    t(n+1) = t(n) + dt;  
    
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(rocketState(3,n));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,rocketState(3,n));

    % Determine Velocities, Reynolds number, Mach number:
    VinAir(:,n) = [rocketState(4,n)-wlong;rocketState(5,n)-wlat;rocketState(6,n)];
    VmaginAir = sqrt(sum(VinAir(:,n).^2));
    Re(n) = rho*VmaginAir*R.lTR/sutherland(Temp);
    Ma(n) = VmaginAir/(sqrt(1.4*287*Temp));
    
    % Assume Zero Angle of Attack on Launch Rail:
    a(n) = real(acos(dot(VinAir(:,n)/VmaginAir,Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*[0;0;1])));    
    
    % Drag Coefficients(IP):
    %[CN, CA, Xcp(n)] = dragCoefs(R,Re(n),a(n),Ma(n),compressibility);
    CN = 0;
    CA = 0;
    
    % Rope Tension: 
    if Ldrogue(n) > Reco.Drogue.L
        Tdrogue(n) = Reco.Drogue.k*(Ldrogue(n)-Reco.Drogue.L)+Vldrogue(n)*Reco.Drogue.c;
    else
        Tdrogue(n) = 0;
    end
    if Lmain(n) > Reco.Main.L
        Tmain(n) = Reco.Main.k*(Lmain(n)-Reco.Main.L)+Vlmain(n)*Reco.Main.c;
    else
        Tmain(n) = 0;
    end
    
    % Rope Forces and Moments:
    if Tdrogue(n) == 0
        dirFt = [0;0;0];
        Ft = [0;0;0];
        Mt = [0;0;0];
    else
        dirFt = (drogueState(1:3,n)-RopeStart(:,n))/Ldrogue(n);
        Ft = dirFt*Tdrogue(n);
        Mt = Tdrogue(n)*cross(RopeStart(:,n)-rocketState(1:3,n),dirFt);
    end
    if Tmain(n) == 0
        dirFtm = [0;0;0];
        Ftm = [0;0;0];
        Mtm = [0;0;0];
    else
        dirFtm = (mainState(1:3,n)-RopeStart(:,n))/Lmain(n);
        Ftm = dirFt*Tmain(n);
        Mtm = Tmain(n)*cross(RopeStart(:,n)-rocketState(1:3,n),dirFtm);
    end
    
    % Determine Forces:
    Fd(n) = .5*CA*rho*R.Aref*VmaginAir^2;
    Fn(n) = .5*CN*rho*R.Aref*VmaginAir^2;
    if a(n) == 0
        % If AoA is zero, no normal force
        dirFn = [0;0;0];
    else
        % Determine direction of Force(cross of moment direction with b3 axis)
        dirFn = cross(cross(Tib(rocketState(7,n),rocketState(8,n),rocketState(9,n))*VinAir(:,n)/VmaginAir,[0;0;1]),[0;0;1]);
        % Normalize Force Direction
        dirFn = dirFn/sqrt(sum(dirFn.^2));
    end
    FR(:,n) = Tbi(rocketState(7,n),rocketState(8,n),rocketState(9,n))*([dirFn(1)*Fn(n);dirFn(2)*Fn(n);R.thrust(n)-Fd(n)])-[0;0;g*R.mass(n)]+Ft+Ftm;
    
    % Determine Moments:
    if a(n) == 0 
        % If angle of attack is zero, no moment(avoids singularity)
        dirM = [0;0;0];
    else
        % Calculate the direction of the Moment(cross of b3 axis with Velocity vector)
        dirM = cross([0;0;1],Tib(rocketState(7,n),rocketState(8,n),rocketState(9,n))*VinAir(:,n)/VmaginAir); % Moment Due to Normal Force
        % Normalize Moment Direction
        dirM = dirM/sqrt(sum(dirM.^2));
    end
    MR(:,n) = [Fn(n)*(R.rocketCentroid{n}(3) + Xcp(end))*dirM(1);Fn(n)*(R.rocketCentroid{n}(3) + Xcp(end))*dirM(2);0] + Mt + Mtm; 

    % Apply Rocket Dynamics using RK4:
    rocketState(:,n+1) = RK4(@rocketN,rocketState(:,n),FR(:,n),MR(:,n),R,n,dt);
    
    % Drogue Velocities:
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(drogueState(3,n));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,drogueState(3,n)); %wind
    
    % Determine Velocities:
    VdrogueinAir(:,n) = drogueState(4:6,n) - [wlong;wlat;0];
    VdroguemaginAir = sqrt(sum(VdrogueinAir(:,n).^2));
    
    % Drogue Forces:
    Fpdrogue(n) = .5*Reco.Drogue.Cd*Reco.Drogue.A*rho*VdroguemaginAir^2;
    FDR(:,n) = -Fpdrogue(n)*VdrogueinAir(:,n)/VdroguemaginAir + [0;0;-Reco.Drogue.m*g] - Tdrogue(n)*dirFt;
    
    % Apply Drogue Dynamics using RK4:
    drogueState(:,n+1) = RK4(@drogueN,drogueState(:,n),FDR(:,n),1,Reco,n,dt);
    
    % Main Velocities:
    % Determine Atmospheric Conditions:
    [Temp,~,rho] = standardAtmosphere(mainState(3,n));
    
    % Get Wind Speed:
    [wlat,wlong] = atmoshwm(LaunchLat,LaunchLong,mainState(3,n)); %wind
    
    % Determine Velocities:
    VmaininAir(:,n) = mainState(4:6,n) - [wlong;wlat;0];
    VmainmaginAir = sqrt(sum(VmaininAir(:,n).^2));
    
    % Main Forces:
    Fpmain(n) = .5*Reco.Main.Cd*Reco.Main.A*rho*VmainmaginAir^2;
    FMR(:,n) = -Fpmain(n)*VmaininAir(:,n)/VmainmaginAir + [0;0;-Reco.Main.m*g] - Tmain(n)*dirFtm;
    
    % Apply Drogue Dynamics using RK4:
    mainState(:,n+1) = RK4(@drogueN,mainState(:,n),FMR(:,n),1,Reco,n,dt);
    
    % Iterator:
    n = n+1;
    
end
t(end) = [];
rocketState(:,end) = [];
drogueState(:,end) = [];
mainState(:,end) = [];

end

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

function N = rocketN(rocketState,F,M,R,n)

N = [rocketState(4);...
    rocketState(5);...
    rocketState(6);...
    F(1)/R.mass(n);...%+rocketState(12)*rocketState(5)-rocketState(11)*rocketState(6);...
    F(2)/R.mass(n);...%+rocketState(10)*rocketState(6)-rocketState(12)*rocketState(4);...
    F(3)/R.mass(n);...%+rocketState(11)*rocketState(4)-rocketState(10)*rocketState(5);...
    sin(rocketState(9))/cos(rocketState(8))*rocketState(11)+cos(rocketState(9))/cos(rocketState(8))*rocketState(12);...
    cos(rocketState(9))*rocketState(11)-sin(rocketState(9))*rocketState(12);...
    rocketState(10)+sin(rocketState(9))*tan(rocketState(8))*rocketState(11)+cos(rocketState(9))*tan(rocketState(8))*rocketState(12);...
    (M(1)+(R.rocketInertia{n}(2,2)-R.rocketInertia{n}(3,3))*rocketState(11)*rocketState(12))/R.rocketInertia{n}(1,1);...
    (M(2)+(R.rocketInertia{n}(3,3)-R.rocketInertia{n}(1,1))*rocketState(10)*rocketState(12))/R.rocketInertia{n}(2,2);...
    (M(3)+(R.rocketInertia{n}(1,1)-R.rocketInertia{n}(2,2))*rocketState(10)*rocketState(11))/R.rocketInertia{n}(3,3)];

end

function N = drogueN(drogueState,F,~,Reco,~)

N = [drogueState(4);...
    drogueState(5);...
    drogueState(6);...
    F(1)/Reco.Drogue.m;...
    F(2)/Reco.Drogue.m;...
    F(3)/Reco.Drogue.m];

end

function N = mainN(drogueState,F,~,Reco,~)

N = [drogueState(4);...
    drogueState(5);...
    drogueState(6);...
    F(1)/Reco.Main.m;...
    F(2)/Reco.Main.m;...
    F(3)/Reco.Main.m];

end

function nextState = RK4(computeN,currentState,F,M,R,n,dt)

f1 = computeN(currentState,F,M,R,n);
f2 = computeN(currentState+dt*f1/2,F,M,R,n);
f3 = computeN(currentState+dt*f2/2,F,M,R,n);
f4 = computeN(currentState+dt*f3,F,M,R,n);

nextState = currentState+dt*(f1/6+(f2+f3)/3+f4/6);

end