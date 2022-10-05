%% Aerodynamics Tester

%% Boomerang Parameters:
R = .15; %length of blade (m)
C = .04; %chord of blade (m)
A = 3*R*C; %area of boomn (m^2)
l = .01; %Boomerang Thickness
m = .04; %mass (kg)
rho = m/(A*l); %density (kg/m^3)
Roffset = .01; %Distance from center of lift to axis of symettry of blade(m)
I = 8e-04; % I33
I = [I/2,0,0;0,I/2,0;0,0,I]; %Inertia Tensor(kgm^2)
n = 3; % 3 bladed boomerang

%% Boomerang Airfoil:
Foil.CL0p = 0.39; 
Foil.CL0m = -0.39; 
Foil.CLap = 4.5; 
Foil.CLam = 4.5;
Foil.CD0p = 0.05; 
Foil.CD0m = 0.05;
Foil.CDap = 0.1;
Foil.CDam = 0.1;

%% New Boomerang Object:
B = boomerang(Foil,R,C,A,Roffset,l,m,rho,I,n);

%% Derive Coefficients:
Xi = .5:.01:1.5;
alpha = -10*pi/180:.01:40*pi/180;

for i = 1:length(Xi)
    for j = 1:length(alpha)
        [CL(i,j), CD(i,j), CMX(i,j), CMY(i,j), CMZ(i,j)] = boomerangCoefs(Xi(i),alpha(j),B);
    end
end

figure;
subplot(5,2,1);
plot(Xi,CL(:,20));
xlabel('Xi'); ylabel('CL');
subplot(5,2,3);
plot(Xi,CD(:,20));
xlabel('Xi'); ylabel('CD');
subplot(5,2,5);
plot(Xi,CMX(:,20));
xlabel('Xi'); ylabel('CMX');
subplot(5,2,7);
plot(Xi,CMY(:,20));
xlabel('Xi'); ylabel('CMY');
subplot(5,2,9);
plot(Xi,CMZ(:,20));
xlabel('Xi'); ylabel('CMZ');

subplot(5,2,2);
plot(alpha,CL(1,:)');
xlabel('\alpha'); ylabel('CL');
subplot(5,2,4);
plot(alpha,CD(1,:)');
xlabel('\alpha'); ylabel('CD');
subplot(5,2,6);
plot(alpha,CMX(1,:)');
xlabel('\alpha'); ylabel('CMX');
subplot(5,2,8);
plot(alpha,CMY(1,:)');
xlabel('\alpha'); ylabel('CMY');
subplot(5,2,10);
plot(alpha,CMZ(1,:)');
xlabel('\alpha'); ylabel('CMZ');