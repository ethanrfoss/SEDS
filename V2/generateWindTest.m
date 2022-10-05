%% Generate Wind Tests:

% syms w
% 
% Lu = 762; Lv = 381; Lw = 381;
% Wind Spectral Densities:
% 
% SDu = sigu^2*2*Lu/pi*(1+8/3*(2.678*Lu*w/V)^2)/(1+(2.678*Lu*w/V)^2)^(11/6);
% 
% SDv = sigv^2*2*Lv/pi*(1+8/3*(2.678*Lv*w/V)^2)/(1+(2.678*Lv*w/V)^2)^(11/6);
% 
% SDw = sigw^2*2*Lw/pi*1/(1+(1.339*Lw*w/V)^2)^(5/6);

% % Test
% 
% Tnu = .0125;
% V = 100;
% bwing = 37.42;
% Lu = 1750;
% Lwv = 1750;
% Lp = sqrt(Lwv*bwing)/2.6;
% sigma = 5;
% 
% ft = logspace(-2,2,200)';
% w = 2*pi*ft;
% 
% tauu = Lu/V;
% tauwv = Lwv/V;
% taup = Lp/V;
% tauq = 4.*bwing/(pi*V);
% taur = 3.*bwing/(pi*V);
% 
% Ku = sigma^2*tauu/pi;
% 
% Snu = 10*log10(sin(w*Tnu/2).^2/(w*Tnu/2).^2);
% 
% Su = 10*log10(Ku./(1+(tauu*w).^2));
% 
% Laplace:
syms s

Hu = sigma*sqrt(2/pi*Lu/V)*(1+.25*Lu/V*s)/(1+1.357*Lu/V*s+.1987*(Lu/V)^2*s^2);


for i = 1:10000
    v(i) = normrnd(0,1);
end

%% Difference Equations

V = 100;
dt = 1/100;
h = 100;
Lu = h/(.177+.000823*h)^1.2;
W20 = sqrt(sum((atmoshwm(LaunchLat,LaunchLong,6)*3.2).^2));
sigu = .1*W20/(.177+.000823*h)^.4;
tauu = Lu/V;

Eu(1) = 0; %sqrt(sum((atmoshwm(LaunchLat,LaunchLong,h)*3.2).^2));
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

%% Calculated Difference Equation:
dt = 1/10;
Lu = 1750;
V = 100;
sigu = .9354;

t = [0:dt:30];
ran = randn(length(t),1);

Eu(1) = 0;
Eu(2) = 0;
for i = 3:length(t)

Eu(i) = 1/(dt^2+2.714*Lu/V*dt+.7948*(Lu/V)^2)*(-Eu(i-1)*(dt^2-1.5896*(Lu/V)^2)-Eu(i-2)*(dt^2-2.714*Lu/V*dt+.7958*(Lu/V)^2)+sigu*sqrt(2/pi*Lu/V)*(ran(i)*(dt^2+.5*Lu/V*dt)+ran(i-1)*2*dt^2+ran(i-2)*(dt^2-.5*Lu/V*dt)));

end

plot(t,Eu)


%% Atmoshwm Test

day = 1:365;
h = 0:100:10000;
FarLat = 35.35;
FarLong = -117.81;
for i = 1:length(day)
    for j = 1:length(h)
        w(i,j) = sqrt(sum(atmoshwm(FarLat,FarLong,h(j),'day',day(i)).^2));
    end
    i
end

figure;
[dayp,hp] = meshgrid(day,h);
surf(dayp',hp',w);