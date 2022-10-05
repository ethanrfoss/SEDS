
function [t,S] = LanderSim

% load([cd '\GimbalData\' 'U' '.mat'],'T','P','U1','U2');
% global T P U1 U2;

global g;
g = 9.81;

global flightParams;
flightParams.F(:,1) = [0;0;0];
flightParams.M(:,1) = [0;0;0];
flightParams.T(1) = 0;

L.m = 385; % Lander Mass Wet
L.I = L.m*[1/12*(3*.8^2+3^2),0,0;0,1/12*(3*.8^2+3^2),0;0,0,.8^2/4]; % Lander Inertia
L.RMETH = .3; % Radius of Methane Tank
L.RLOX = .3; % Radius of Lox Tank
L.LU = 2.7; % Length to Universal Joint
L.LE = .2; % Length of Engine
L.LMETH = 1; % Length to Methane Tank Center
L.LLOX = 2.1; % Length to Lox Tank Center
L.rhoLOX = 1141; % LOX Density
L.rhoMETH = 424; % METH Density
L.mLOX = 99.154; % Mass of LOX
L.mMETH = 36.846; % Mass of METH
L.cTLOX = 0.000247885; % LOX Relation Between Thrust and Mass Flow 
L.cTMETH = 0.000092115; % METH Relation Between Thrust and Mass Flow 

S(:,1) = [0 % X[m] 1
          0 % Y[m] 2
          0 % Z[m] 3
          0 % Vx[m] 4
          0 % Vy[m] 5
          0 % Vz[m] 6
          0 % alpha[rad] 7
          0 % beta[rad] 8 
          0 % gamma[rad] 9
          0 % w1[rad/s] 10
          0 % w2[rad/s] 11
          0 % w3[rad/s] 12
          L.m % mass[kg] 13
          L.I(1,1) % XX Inertia[kgm^2] 14
          L.I(2,2) % YY Inertia[kgm^2] 15
          L.I(3,3) % ZZ Inertia[kgm^2] 16 
          .187357 % Height Displaced in LOX Tank 17
          .187357 % Height Displaced in METH Tank 18
          1.7]; % COM of Rocket 19

dt = 1/1000;
C = 0;
i = 1;
t(1) = 0;

while t<200
    S(:,i+1) = RK4(S(:,i),C,L,dt,t);
    t(i+1) = t(i) + dt;
    i = i + 1;
end

figure;
plot3(S(1,:),S(2,:),S(3,:)); axis equal;

end

function [Sdot,F,M,T] = EOM(S,L,C,t)

[F,M,T] = Forces(S,L,C,t);

mdotLOX = -L.cTLOX*T;
mdotMETH = -L.cTMETH*T;

hdotLOX = 3*mdotLOX/(pi*L.rhoLOX*(-6*L.RLOX*S(17)+3*S(17)^2)); %-6*mdotLOX/(pi*L.rhoLOX*(3*L.RLOX^2+3*S(17)^2));
hdotMETH = 3*mdotMETH/(pi*L.rhoMETH*(-6*L.RMETH*S(18)+3*S(18)^2));

PIX = [0.4701,-2.3314,4.1022,-3.0241,-0.0686,1.6761]; % Polynomial Describing X-axis Inertia of Cutoff Sphere
PIY = [0.4701,-2.3314,4.1022,-3.0241,-0.0686,1.6761]; % Polynomial Describing Y-axis Inertia of Cutoff Sphere
PIZ = [-.3104,1.4885,-1.8111,-.2819,.0308,1.6761]; % Polynomial Describing Z-axis Inertia of Cutoff Sphere

Sdot = [S(4)
        S(5)
        S(6)
        F(1)/S(13)
        F(2)/S(13)
        F(3)/S(13)
        sin(S(9))/cos(S(8))*S(11)+cos(S(9))/cos(S(8))*S(12)
        cos(S(9))*S(11)-sin(S(9))*S(12)
        S(10)+sin(S(9))*tan(S(8))*S(11)+cos(S(9))*tan(S(8))*S(12)
        (M(1)+(S(15)-S(16))*S(11)*S(12))/S(14)
        (M(2)+(S(16)-S(14))*S(10)*S(12))/S(15)
        (M(3)+(S(14)-S(15))*S(10)*S(11))/S(16)
        mdotLOX+mdotMETH
        1/39.3483*(-hdotLOX*L.rhoLOX*(PIX(1)*L.RLOX^5*0+PIX(2)*L.RLOX^4*1+PIX(3)*L.RLOX^3*2*S(17)+PIX(4)*L.RLOX^2*3*S(17)^2+PIX(5)*L.RLOX^1*4*S(17)^3+PIX(6)*L.RLOX^0*5*S(17)^4)-hdotMETH*L.rhoMETH*(PIX(1)*L.RMETH^5*0+PIX(2)*L.RMETH^4*1+PIX(3)*L.RMETH^3*2*S(18)+PIX(4)*L.RMETH^2*3*S(18)^2+PIX(5)*L.RMETH^1*4*S(18)^3+PIX(6)*L.RMETH^0*5*S(18)^4))
        1/39.3483*(-hdotLOX*L.rhoLOX*(PIY(1)*L.RLOX^5*0+PIY(2)*L.RLOX^4*1+PIY(3)*L.RLOX^3*2*S(17)+PIY(4)*L.RLOX^2*3*S(17)^2+PIY(5)*L.RLOX^1*4*S(17)^3+PIY(6)*L.RLOX^0*5*S(17)^4)-hdotMETH*L.rhoMETH*(PIY(1)*L.RMETH^5*0+PIY(2)*L.RMETH^4*1+PIY(3)*L.RMETH^3*2*S(18)+PIY(4)*L.RMETH^2*3*S(18)^2+PIY(5)*L.RMETH^1*4*S(18)^3+PIY(6)*L.RMETH^0*5*S(18)^4))
        1/33.74*(-hdotLOX*L.rhoLOX*(PIZ(1)*L.RLOX^5*0+PIZ(2)*L.RLOX^4*1+PIZ(3)*L.RLOX^3*2*S(17)+PIZ(4)*L.RLOX^2*3*S(17)^2+PIZ(5)*L.RLOX^1*4*S(17)^3+PIZ(6)*L.RLOX^0*5*S(17)^4)-hdotMETH*L.rhoMETH*(PIZ(1)*L.RMETH^5*0+PIZ(2)*L.RMETH^4*1+PIZ(3)*L.RMETH^3*2*S(18)+PIZ(4)*L.RMETH^2*3*S(18)^2+PIZ(5)*L.RMETH^1*4*S(18)^3+PIZ(6)*L.RMETH^0*5*S(18)^4))
        hdotLOX
        hdotMETH
        3/4*hdotLOX*(1-L.RLOX^2/(L.RLOX^2+S(17)^2)) + 3/4*hdotMETH*(1-L.RMETH^2/(L.RMETH^2+S(18)^2))];
  
end

function [F,M,T] = Forces(S,L,C,t)

global g;

% if mod(t,C.dt) == 0
%     [u1,u2,valveAngle] = Control(S,L,C);
% end
if S(17) < 2*L.RLOX && S(18) < 2*L.RMETH
    T = 4000;
else
    T = 0;
end

phi = 0; theta = 0;

% [phi,theta] = GimbalKinematics(u1,u2);
% T = EnginePhysics(valveAngle);

Fx = -cos(theta)*sin(phi)*T;
Fy = -sin(theta)*sin(phi)*T;
Fz = cos(phi)*T;

F = Tb2i(S(7),S(8),S(9))*[Fx;Fy;Fz] + [0;0;-g]; % Body Forces

M = [Fy*(L.LU-S(19)+L.LE*cos(phi));-Fx*(L.LU-S(19)+L.LE*cos(phi));0];



end

function Snext = RK4(S,C,L,dt,t)

[f1,~,~,~] = EOM(S,L,C,t);
[f2,~,~,~] = EOM(S+dt*f1/2,L,C,t);
[f3,~,~,~] = EOM(S+dt*f2/2,L,C,t);
[f4,F,M,T] = EOM(S+dt*f3,L,C,t);

Snext = S+dt*(f1/6+(f2+f3)/3+f4/6);

global flightParams;
flightParams.F(:,end+1) = F;
flightParams.M(:,end+1) = M;
flightParams.T(end+1) = T;

end

function T = Ti2b(a,b,g)

T = [cos(a)*cos(b),cos(b)*sin(a),-sin(b); cos(a)*sin(b)*sin(g)-cos(g)*sin(a),cos(a)*cos(g)+sin(a)*sin(b)*sin(g),cos(b)*sin(g); sin(a)*sin(g)+cos(a)*cos(g)*sin(b),cos(g)*sin(a)*sin(b)-cos(a)*sin(g),cos(b)*cos(g)]; 

end

%% Body to Inertial Frame Transformation:
function T = Tb2i(a,b,g)

T = Ti2b(a,b,g)';

end

% function Re = GimbalKinematics(u1,u2)
% 
% global T P U1 U2;
% 
% if u1 == 0
%     
%     
% end
% 
% function T = EngineSystem(valveMETH,valveLOX)
