% rocketTrajectoryTester

%% Time step:
dt = 1/1000;

%% Rocket Object:
% Test:
%R = rocket(1,.25,.25,.25,.25,.125,4,1,3.5,1.4,.5,1,1,.01,2,1,.2,16.3231,3.5,1,.2,16.3231,10,30,dt);

% Purdue:
%R = rocket(1,.16256,.16256,.16256,.16256,.16256/2,4.5,.3,5.1,.4,.2,.4,.2,.01,2,1,.16256,14,3.5,1,.16256,14,6,65,dt);

% Halya:
r = .2;
R = rocket(1,r,r,r,r,r/2,4.5,.3,5.1,.4,.2,.4,.2,.01,2,1,r,9.3461,3.5,1,r,3.38674,4.6,55,dt);

% Vulcan
% r = .2286;
% R = rocket(1,r,r,r,r,r/2,4.5,.3,5.1,.4,.2,.4,.2,.01,2,1,r,11.55,3.5,1,r,4.84,11,113,dt);

%% Launch Rail Height:
H = 15; %(m)

%% Recovery System:
Reco = recovery(50,100*10^2,.5,.97,.25,.35,60,100*10^2,.5,1.75,1.5,1,1000);

%% Rocket Dimensions:
R.displayRocketDimensions(R);

%% Run Simulation:
rocketTrajectory2(R,Reco,H,35.35053,-117.80814,1000,'Laitone-Ackeret',dt);

generatePlots(t,rocketState,drogueState,mainState,a,Ma,R,nDrogueDeploy,Reco,Ldrogue,Tdrogue,Vldrogue,nMainDeploy,Lmain,Tmain,Vlmain,true,true);

rocketParachuteAnimation(t,R,Reco,rocketState,drogueState,mainState,a,FR,nDrogueDeploy,nMainDeploy,RopeStart,H,createGif,saveName)