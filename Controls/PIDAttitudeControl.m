%% PID for Attitude Control:

%% System:

syms s
syms L real
syms m real
syms I real
syms g real

A = [0 1 0 0;m*g*L/I 0 0 0;0 0 0 1;-m*g*L^2/I 0 0 0];
B = [0 -L/I 0 (m*L^2+I)/(I*m)]';
D = 0;

utophi = [1 0 0 0]*(s*eye(4,4)-A)^-1*B;
utox = [0 0 1 0]*(s*eye(4,4)-A)^-1*B;

%% Numerical Substitution and Transfer Function:
L = .5; m = .8; I = 1.2; g = 9.81;

PU = minreal(syms2tf(subs(utophi)));
XU = minreal(syms2tf(subs(utox)));

G1 = PU;
G2 = minreal(XU/PU);

%% Control Design:
% Inner Loop:
%sisotool(G1);
D1 = -50*tf([1 1],[.005 1]);
P = 1/1.18;

% Outer Loop:
%sisotool(G2);
D2 = .042106*tf([2 1],[.59 2]);

% T1:
y = tf([-50 -50],1);
x = tf([.005 1],1);
[b a] = tfdata(G1, 'v');
b = tf(b,1);
a = tf(a,1);
T1 = b*y/(a*x+b*y);

% T2:
[num dem] = tfdata(D2*P*T1*G2, 'v');
num = tf(num,1);
dem = tf(dem,1);
T2 = num/(dem+num);

%% Step Response of Theta and X:
t = 0:.01:10;
theta = step(D2*P*D1*G1/(1+G2*D2*P*D1*G1+D1*G1),t);
x = step(T2,t);

figure(1); sgtitle('Rocket 2 DOF Response');
subplot(2,1,1);
plot(t,theta*180/pi);
xlabel('Time[s]'); ylabel('Rocket Angle[deg]');
subplot(2,1,2);
plot(t,x);
xlabel('Time[s]'); ylabel('Position[m]');

%% Impulse Response of Theta and X:
t = 0:.01:10;
theta = impulse(1/10*D2*P*D1*G1/(1+G2*D2*P*D1*G1+D1*G1),t);
x = impulse(1/10*T2,t);

figure(1); sgtitle('Rocket 2 DOF Response');
subplot(2,1,1);
plot(t,theta*180/pi);
xlabel('Time[s]'); ylabel('Rocket Angle[deg]');
subplot(2,1,2);
plot(t,x);
xlabel('Time[s]'); ylabel('Position[m]');

%% Symbolic to TF:
function[TF] = syms2tf(G)
[symNum,symDen] = numden(G); %Get num and den of Symbolic TF
TFnum = sym2poly(symNum);    %Convert Symbolic num to polynomial
TFden = sym2poly(symDen);    %Convert Symbolic den to polynomial
TF =tf(TFnum,TFden);
end