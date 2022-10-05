function [CN CA Xcp] = dragCoefs(R,Re,alpha,Ma,compressibility)

alphaX = [0 4 6 8 10 12 14 16 18 20];
deltaY = [0 .79 .86 .92 .935 .95 .96 .97 .98 .98];
etaY = [0 .6 .63 .66 .68 .71 .725 .74 .75 .76];

Afe = .5*(R.lr+R.lt)*R.ls; %Exposed area
Afp = Afe+.5*R.df*R.lr; %Planform Area

%% CN
% Nose Cone
CNan = 2;
Xcpn = 2/3*R.ln; % conical nose cone

% Conical Change in body diameter
CNac = 2*((R.dd/R.dn)^2-(R.du/R.dn)^2);
Xcpc = R.Xc + R.lc/3*(1+(1-R.du/R.dd)/(1-(R.du/R.dd)^2));

% Fins
Kfb = 1+(R.df/2)/(R.ls+R.df/2);
CNaf = Kfb*(4*R.n*(R.ls/R.dn)^2)/(1+sqrt(1+(2*R.lm/(R.lr+R.lt))^2));
Xcpf =R.Xf+R.lm*(R.lr+2*R.lt)/(3*(R.lr+R.lt))+1/6*(R.lr+R.lt-R.lr*R.lt/(R.lr+R.lt));

% Rocket Body Lift
K = 1; %Between 1 and 1.5
Ap = .5*R.ln*R.dn + R.lb*R.db + .5*(R.du+R.dd)*R.lc; %Planform area(conical nose)
Ar = pi/4*(R.dn)^2;
CNal = K*Ap/Ar*alpha;
Xcpl = ((.5*R.ln*R.dn)*(2/3*R.ln)+(R.lb*R.db)*(R.Xb+.5*R.lb)+(R.Xc+R.lc/2*((R.du+2*R.dd)/(R.du+R.dd)))*(.5*(R.du+R.dd)*R.lc));

% Total
CNa = CNan+CNac+CNaf+CNal;
CN = CNa*alpha;
Xcp = (CNan*Xcpn+CNac*Xcpc+CNaf*Xcpf+CNal*Xcpl)/CNa;

%% CD
% Zero Angle of Attack
% Viscous Friction(body and fin)
Rec = 5*10^5;
Ref = Re*R.lm/R.lTR; %Reynolds number of fins
if Re<Rec
    Cfb = 1.328/sqrt(Re);
else
    Bb = Rec*(.074/Re^.2-1.328/sqrt(Re));
    Cfb = .074/Re^.2-Bb/Re;
end
if Ref<Rec
    Cff = 1.328/sqrt(Ref);
else
    Bf = Rec*(.074/Ref^.2-1.328/sqrt(Ref));
    Cff = .074/Ref^.2-Bf/Ref;
end
Cffb = Cfb+Cff;

% Body Drag
CDfb = (1+60/(R.lTR/R.db)^3+.0025*R.lb/R.db)*(2.7*R.ln/R.db+4*R.lb/R.db+2*(1-R.dd/R.db)*R.lc/R.db)*Cffb;

% Base Drag
CDb = real(.029*(R.dd/R.db)^3/sqrt(CDfb));

% Fin Drag
CDf = 2*Cff*(1+2*R.Tf/R.lm)*4*R.n*Afp/(pi*R.df^2);

%Iterference Drag
CDi = 2*Cff*(1+2*R.Tf/R.lm)*4*R.n*(Afp-Afe)/(pi*R.df^2);

% Coefficient of Drag at zero angle of attack
CD0 = CDfb+CDb+CDf+CDi;

% Alpha Drag
% Body
CDba = 2*interp1(alphaX,deltaY,alpha)*alpha^2+3.6*interp1(alphaX,etaY,alpha)*(1.36*R.lTR-.55*R.ln)*alpha^3/(pi*R.db);

% Fins
Rs = R.lTS/R.df;
kfb = .8065*Rs^2+1.1553*Rs;
kbf = .1935*Rs^2+.8174*Rs+1;
CDfa = alpha^2*(1.2*Afp*4/(pi*R.df^2)+3.12*(kfb+kbf-1)*Afe*4/(pi*R.df^2));

% Total
CD = CD0 +CDba +CDfa; 

% Axial Drag coefficient(drag coefficient along rocket)

CA = (CD*cos(alpha)-.5*CN*sin(2*alpha))/(1-sin(alpha)^2);

%% Compressible flow Correction
switch compressibility

    case 'Prandtl-Glauert'
        %Prandtl-Glauert
        % Worst Case Scenario, probably the best to use
        if Ma<=.8
            CN = CN/sqrt(1-Ma^2);
            CA = CA/sqrt(1-Ma^2);
        elseif Ma<1.1
            % Coeficients Given by: [.8^2 .8 1;1.1^2 1.1 1;2*.8 1 0]^-1*[1/sqrt(1-.8^2);1/sqrt(1.1^2-1);.8/(1-.8^2)^(3/2)]
            % Parabolic Fit at Ma 1.1
            aN = -43.8295*CN;
            bN = 84.9943*CN;
            cN = -38.2779*CN;
            aA = -43.8295*CA;
            bA = 84.9943*CA;
            cA = -38.2779*CA;
    
            CN = aN*Ma^2+bN*Ma+cN;
            CA = aA*Ma^2+bA*Ma+cA;
        else
            CN = CN/sqrt(Ma^2-1);
            CA = CA/sqrt(Ma^2-1);
        end
        
        
    case 'Karman-Tsien'
        % Karman-Tsien
        % Pretty much garbage
        Mcr = .8;
        if Ma<=Mcr
            CN = CN/(sqrt(1-Ma^2)+(Ma^2/(1+sqrt(1-Ma^2)))*CN/2);
            CA = CA/(sqrt(1-Ma^2)+(Ma^2/(1+sqrt(1-Ma^2)))*CA/2);
        else
            CN = CN/(sqrt(1-Mcr^2)+(Mcr^2/(1+sqrt(1-Mcr^2)))*CN/2);
            CA = CA/(sqrt(1-Mcr^2)+(Mcr^2/(1+sqrt(1-Mcr^2)))*CA/2);
        end
        
    case 'Laitone-Ackeret'
        
        % Laitone + Ackeret
        % Supposed to be more accurate maybe? But better case scenario than
        % prandtl
        gamma = 1.4;
        Laitone = @(M,C) C/(sqrt(1-M^2)+(M^2*(1+(gamma-1)/2*M^2))/(2*(1+sqrt(1-M^2)))*C);
        if Ma<.8
            CN = Laitone(Ma,CN);
            CA = Laitone(Ma,CA);
        elseif Ma < 1.1
            % Coeficients Given by: [.9^2 .9 1;1.1^2 1.1 1;2*1.1 1 0]^-1*[1/sqrt(1-.9^2);1/sqrt(1.1^2-1);-1.1/(1.1^2-1)^(3/2)]
            % Parabolic Fit
            aN = -43.8295*Laitone(.8,CN)*sqrt(1-.8^2);
            bN = 84.9943*Laitone(.8,CN)*sqrt(1-.8^2);
            cN = -38.2779*Laitone(.8,CN)*sqrt(1-.8^2);
            aA = -43.8295*Laitone(.8,CA)*sqrt(1-.8^2);
            bA = 84.9943*Laitone(.8,CA)*sqrt(1-.8^2);
            cA = -38.2779*Laitone(.8,CA)*sqrt(1-.8^2);
    
            CN = aN*Ma^2+bN*Ma+cN;
            CA = aA*Ma^2+bA*Ma+cA;
        else
            CN = Laitone(.8,CN)*sqrt(1-.8^2)/sqrt(Ma^2-1);
            CA = Laitone(.8,CA)*sqrt(1-.8^2)/sqrt(Ma^2-1);
        end
        
    otherwise
        disp('No Compressibility Specified, Defaulting to incompressible');
        
end

end