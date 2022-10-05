syms phir phim thetam s f1 h1 f2 h2 z A2 A1 b2 a2 b1 a1 P;

D2 = A2*(a2*s+1)/(b2*s+1);
D2 = subs(D2,{s},{2/(f1*h1)*(z-1)/(z+1)});
D1 = A1*(a1*s+1)/(b1*s+1);
D1 = subs(D1,{s},{2/(f2*h2)*(z-1)/(z+1)});

u = ((phir-phim)*D2*P-thetam)*D1;

[N, D] = numden(uz);
Nc = eval(coeffs(N)); %Get coeffs and evaluatle symbolic variable, i.e. make real matrix
Dc = eval(coeffs(D));
Nc = Nc./(Dc(1)); %Turn into proper polynomial, first coeffs of a is 1
Dc = Dc./(Dc(1));
M = idpoly(Dc, Nc, 'NoiseVariance',0)

%%

D2 = .0071465*tf([1 1],[.1 1]);
D1 = -100*tf([.18065 1],[6.7 1]);
P = 1/1.21;

u = ((phir-phim)*D2*P-thetam)*D1;

%%

D2 = .0071465; D1 = -100;
a2 = 1; b2 = .1; a1 = .1806; b1 = 6.7;
P = 1/1.21;
h2 = 1/20; w2 = 1.21; f2 = 2*(1-cos(w2*h2))/(w2*h2*sin(w2*h2));
h1 = 1/100; w1 = 38.3; f1 = 2*(1-cos(w2*h2))/(w2*h2*sin(w2*h2));

syms u(n) phi(n) theta(n)
RHS = (-u(n-1)*(2*f1*h1*f2*h2-8*b2*b1)-u(n-2)*(f2*h2-2*b2)*(f1*h1-2*b1)+phi(n)*(2*a2+f2*h2)*(2*a1+f1*h1)*D2*D1*P+phi(n-1)*(2*f1*h1*f2*h2-8*a2*a1)*D2*D1*P+phi(n-2)*(f2*h2-2*a2)*(f1*h1-2*a1)*D2*D1*P-theta(n)*(2*a1+f1*h1)*(2*b2+f2*h2)*D1-theta(n-1)*((2*a1+f1*h1)*(f2*h2-2*b2)+(f1*h1-2*a1)*(2*b2+f2*h2))*D1-theta(n-2)*(f1*h1-2*a1)*(f2*h2-2*b2)*D1)/((2*b2+f2*h2)*(2*b1+f1*h1));