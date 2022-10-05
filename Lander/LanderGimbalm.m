
%% Transformation of engine to body frame from paramaters phi, theta(from cylindrical coordinates)
Te2b = @(phi,theta) [cos(phi) - sin(theta)^2*(cos(phi) - 1),sin(2*theta)*(cos(phi)/2 - 1/2),cos(theta)*sin(phi);
                     sin(2*theta)*(cos(phi)/2 - 1/2),(cos(phi) - 1)*sin(theta)^2 + 1,sin(phi)*sin(theta);
                     -cos(theta)*sin(phi),-sin(phi)*sin(theta),cos(phi)];          

%% Linkage System Parameters:
rb = 1; %Radius of Rocket Spherical joints from center
lb = .2; %Distance of Universal Joint from spherical joint plane
le = 1.5; %Length of Universal Joint to Engine spherical joint plane
re = .3; %Radius of Engine Spherical Joints from center

%% Actuator Vectors as a function of phi,theta:
u1 = @(phi,theta) [-rb;0;lb] + Te2b(phi,theta)*[re;0;le];
u2 = @(phi,theta) [0;-rb;lb] + Te2b(phi,theta)*[0;re;le];

%% Determine COntrol Input Lengths at Different phi,theta
theta = [0:.01:2*pi];
phi = [0:.01:pi/4];
u1nom = norm(u1(0,0));
u2nom = norm(u2(0,0));


for i = 1:length(theta)
    for j = 1:length(phi)
        
        U1(i,j) = norm(u1(phi(j),theta(i)))-u1nom;
        U2(i,j) = norm(u2(phi(j),theta(i)))-u2nom; 
    end
end

save([cd '\GimbalData\' 'U' '.mat'],'T','P','U1','U2');
    
%% Plot Theta,Phi vs. U1,U2
[T,P] = meshgrid(theta,phi);
T = T'; P = P';
figure(1);
subplot(2,1,1);
SU1 = surf(T*180/pi,P*180/pi,U1); shading interp;
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Control Input Length of First Actuator[m]');
subplot(2,1,2);
SU2 = surf(T*180/pi,P*180/pi,U2); shading interp;
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Control Input Length of Second Actuator[m]');

%% Fit U1, U2 to Polynomial:
U1fit = fit([T(:),P(:)],U1(:),'poly55');
U2fit = fit([T(:),P(:)],U2(:),'poly55');

%% Plot Fit Errors:
figure(2);
subplot(2,1,1);
SU1e = surf(T*180/pi,P*180/pi,U1-U1fit(T,P)); shading interp;
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Fit Difference[m]');
subplot(2,1,2);
SU2e = surf(T*180/pi,P*180/pi,U2-U2fit(T,P)); shading interp;
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Fit Difference[m]');

%% Fit U1,U2 to Sin:
U1sin = @(phi,theta) -phi/(pi/4).*(-(min(U1(:))-max(U1(:)))/2*cos(theta)-(max(U1(:))+min(U1(:)))/2);
U2sin = @(phi,theta) -phi/(pi/4).*(-(min(U2(:))-max(U2(:)))/2*sin(theta)-(max(U2(:))+min(U2(:)))/2);

%% Plot Fit Errors:
figure(3);
subplot(2,1,1);
SU1e = surf(T*180/pi,P*180/pi,U1-U1sin(P,T)); shading interp;
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Fit Difference[m]');
subplot(2,1,2);
SU2e = surf(T*180/pi,P*180/pi,U2-U2fit(P,T)); shading interp;
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Fit Difference[m]');



%% Determine phi and theta values from u1 and u2 control inputs:

u1p = 0;
u2p = -0;

C1 = contours(T,P,U1-u1p,[0 0]);
x1 = C1(1,2:end); y1 = C1(2,2:end);
C2 = contours(T,P,U2-u2p,[0 0]);
x2 = C2(1,2:end); y2 = C2(2,2:end);

[xi,yi] = intersections(x1,y1,x2,y2);

figure(4);
subplot(2,2,[1 3]); hold on;
plot(x1*180/pi,y1*180/pi,'r');
plot(x2*180/pi,y2*180/pi,'r'); 
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); axis([0 360 0 45]);
title("\Theta = " + num2str(xi*180/pi) + "[deg] || \Phi = " + num2str(yi*180/pi) + "[deg]");

subplot(2,2,2); hold on;
surf(T*180/pi,P*180/pi,U1); shading interp;
surf(T*180/pi,P*180/pi,u1p*ones(size(U1)));
plot3(x1*180/pi,y1*180/pi,zeros(1,length(x1)),'r');
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Control Input Length of First Actuator[m]'); axis([0 360 0 45 min(U1(:)),max(U1(:))]);

subplot(2,2,4); hold on;
surf(T*180/pi,P*180/pi,U2); shading interp;
surf(T*180/pi,P*180/pi,u2p*ones(size(U1)));
plot3(x2*180/pi,y2*180/pi,zeros(1,length(x2)),'r');
xlabel('\Theta[deg]'); ylabel('\Phi[deg]'); zlabel('Control Input Length of Second Actuator[m]'); axis([0 360 0 45 min(U2(:)),max(U2(:))]);

%% Determine x and y displacements from u1 and u2 control inputs:

u1p = 0;
u2p = 0;

C1 = contours(cos(T).*sin(P),sin(T).*sin(P),U1-u1p,[0 0]);
x1 = C1(1,2:end); y1 = C1(2,2:end);
C2 = contours(cos(T).*sin(P),sin(T).*sin(P),U2-u2p,[0 0]);
x2 = C2(1,2:end); y2 = C2(2,2:end);

[xi,yi] = intersections(x1,y1,x2,y2);

figure(4);
subplot(2,2,[1 3]); hold on;
plot(x1,y1,'r');
plot(x2,y2,'r'); 
xlabel('X[m]'); ylabel('Y[m]'); axis([min(cos(T(:)).*sin(P(:))) max(cos(T(:)).*sin(P(:))) min(sin(T(:)).*sin(P(:))) max(sin(T(:)).*sin(P(:))) min(U1(:)),max(U1(:))]);
title("X = " + num2str(xi) + "[m] || Y = " + num2str(yi) + "[m]");

subplot(2,2,2); hold on;
surf(cos(T).*sin(P),sin(T).*sin(P),U1); shading interp;
surf(cos(T).*sin(P),sin(T).*sin(P),u1p*ones(size(U1)));
plot3(x1,y1,zeros(1,length(x1)),'r');
xlabel('X[m]'); ylabel('Y[m]'); axis([min(cos(T(:)).*sin(P(:))) max(cos(T(:)).*sin(P(:))) min(sin(T(:)).*sin(P(:))) max(sin(T(:)).*sin(P(:))) min(U2(:)),max(U2(:))]);

subplot(2,2,4); hold on;
surf(cos(T).*sin(P),sin(T).*sin(P),U2); shading interp;
surf(cos(T).*sin(P),sin(T).*sin(P),u2p*ones(size(U1)));
plot3(x2,y2,zeros(1,length(x2)),'r');
xlabel('X[m]'); ylabel('Y[m]'); axis([min(cos(T(:)).*sin(P(:))) max(cos(T(:)).*sin(P(:))) min(sin(T(:)).*sin(P(:))) max(sin(T(:)).*sin(P(:)))]);



      