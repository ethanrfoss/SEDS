function [Eu,Ev,Ew] = generateWind(rocketState,Eu,Ev,Ew,dt,seeverity)

%% Velocity and Height of Rocket:
V = sqrt(sum(rocketState(4:6).^2))*3.28084; %Feet/sec
h = rocketState(3)*3.28084; %Feet

%% Length Parameters and Sigmas of Wind:
if h < 1000
    Lu = h;
    Lv = h./(.177+.000823*h).^1.2;
    Lw = h./(.177+.000823*h).^1.2;
elseif h<2000
    Lu = 1.5*h-1000;
    Lv = 1.5*h-1000;
    Lw = 1.5*h-1000;
else
    Lu = 2500;
    Lv = 2500;
    Lw = 2500;
end


W20 = sqrt(sum((atmoshwm(LaunchLat,LaunchLong,6)).^2));
sigu = .1*W20/(.177+.000823*h)^.4;
tauu = Lu/V;

Eu(1) = 0;
Ev(1) = 0;
Ew(1) = 0;
t = [0:dt:200];
for i = 1:length(t)
    Eu(i+1) = (1-dt/tauu)*Eu(i) + sigu*sqrt(2*dt/tauu)*normrnd(0,1);

    Ev(i+1) = (1-2*dt/tauu)*Ev(i) + sigu*sqrt(4*dt/tauu)*normrnd(0,1);

    Ew(i+1) = (1-2*dt/tauu)*Ew(i) + sigu*sqrt(4*dt/tauu)*normrnd(0,1);
    V = V + .001;
    tauu = Lu/V;
end

figure(1); hold on;
plot(t,Eu(1:end-1));
% plot(t,Ev(1:end-1));
% plot(t,Ew(1:end-1));
