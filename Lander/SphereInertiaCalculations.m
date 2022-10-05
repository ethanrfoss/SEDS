
[X,Y,Z] = sphere(100);
[r,c] = size(Z);
 warning off
for i = 1:r-2
    Xt = X(1:r-i,:);
    Yt = Y(1:r-i,:);
    Zt = Z(1:r-i,:);
    
    Xt = reshape(Xt,numel(Xt),1);
    Yt = reshape(Yt,numel(Yt),1);
    Zt = reshape(Zt,numel(Zt),1);
    
    Ct = convhull(Xt,Yt,Zt);
    [RBP,~] = RigidBodyParams(triangulation(Ct,Xt,Yt,Zt));
    
    I1(i) = RBP.inertia_tensor(1,1)+RBP.centroid(3)^2*RBP.volume;
    I2(i) = RBP.inertia_tensor(2,2)+RBP.centroid(3)^2*RBP.volume;
    I3(i) = RBP.inertia_tensor(3,3);
    i
end
warning on
figure(1); hold on;
plot(Z(1:r-2,1)+1,I1);
plot(Z(1:r-2,1)+1,I2);
plot(Z(1:r-2,1)+1,I3);

ZPs = polyfit(Z(1:r-2,1)+1,I3,5)
XPs = polyfit(Z(1:r-2,1)+1,I2,5)

plot(Z(1:r-2,1)+1,polyval(ZPs,Z(1:r-2,1)+1))
plot(Z(1:r-2,1)+1,polyval(XPs,Z(1:r-2,1)+1))

%% Inertia Test:

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

PIX = [0.4701,-2.3314,4.1022,-3.0241,-0.0686,1.6761]; % Polynomial Describing X-axis Inertia of Cutoff Sphere

T = 4000;

hLOX(1) = .01;
Iz(1) = 0;
dt = .001;
i = 1;

while hLOX(i) < 2*L.RLOX

    mdotLOX = -L.cTLOX*T;
    hdotLOX = 3*mdotLOX/(pi*L.rhoLOX*(-6*L.RLOX*hLOX(i)+3*hLOX(i)^2));
    Izdot = -hdotLOX*L.rhoLOX*(PIX(1)*L.RLOX^5*0+PIX(2)*L.RLOX^4*1+PIX(3)*L.RLOX^3*2*hLOX(i)+PIX(4)*L.RLOX^2*3*hLOX(i)^2+PIX(5)*L.RLOX^1*4*hLOX(i)^3+PIX(6)*L.RLOX^0*5*hLOX(i)^4);
    
    hLOX(i+1) = hLOX(i) + dt*hdotLOX;
    Iz(i+1) = Iz(i) + Izdot*dt;
    i = i + 1;
end