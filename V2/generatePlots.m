function generatePlots(t,rocketState,drogueState,mainState,a,Ma,R,nDrogueDeploy,Reco,Ldrogue,Tdrogue,Vldrogue,nMainDeploy,Lmain,Tmain,Vlmain,massReport,recoReport)

%% Initialize Inputs:
X = rocketState(1,:); Y = rocketState(2,:); Z = rocketState(3,:);
Vmag = sqrt(sum(rocketState(4:6,:).^2)); Amag = diff(Vmag)*1000;
alpha = rocketState(7,:); beta = rocketState(8,:); gamma = rocketState(9,:);

Zdrogue = drogueState(3,:);
Vmagdrogue = sqrt(sum(drogueState(4:6,:).^2));
Amagdrogue = diff(Vmagdrogue)*1000;

Zmain = mainState(3,:);
Vmagmain = sqrt(sum(mainState(4:6,:).^2));
Amagmain = diff(Vmagmain)*1000;

%% Rocket Data Plots:
plots = figure('WindowState','maximized');
set(gcf,'Color',[1 1 1]);
set(gca, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
sgtitle('Rocket Ascent');

subplot(3,4,[1 5 9]); hold on;
title('Rocket Trajectory'); axis equal; view(45,25);
plot3(X(1:nDrogueDeploy),Y(1:nDrogueDeploy),Z(1:nDrogueDeploy), 'LineWidth', 2,'Color','b'); 
plot3(X(1),Y(1),Z(1),'g*'); plot3(X(nDrogueDeploy),Y(nDrogueDeploy),Z(nDrogueDeploy),'r*');
axis([min(X(1:nDrogueDeploy)) max(X(1:nDrogueDeploy)) min(Y(1:nDrogueDeploy)) max(Y(1:nDrogueDeploy)) min(Z(1:nDrogueDeploy)) max(Z(1:nDrogueDeploy))]);
xlabel('Easting[m]'); ylabel('Northing[m]'); zlabel('Altitude[m]');
hold off;

subplot(3,4,2); hold on;
title(sprintf('Vertical Distance(%.1fft Apogee)',(max(Z)-Z(1))*3.28084));
plot(t(find(Z == max(Z))),Z(find(Z==max(Z)))-Z(1),'r*');
plot(t(1:nDrogueDeploy),Z(1:nDrogueDeploy)-Z(1),'b', 'LineWidth', 2);
xlabel('Time[s]'); ylabel('Altitude[m]');
hold off;

subplot(3,4,3); hold on;
title('Velocity Magnitude');
plot(t(1:nDrogueDeploy),Vmag(1:nDrogueDeploy), 'LineWidth', 2);
xlabel('Time[s]'); ylabel('Velocity[m/s]');
hold off;

subplot(3,4,4); hold on;
title('Acceleration Magnitude');
plot(t(1:nDrogueDeploy),Amag(1:nDrogueDeploy), 'LineWidth', 2);
xlabel('Time[s]'); ylabel('Acceleration[m/s^2]');
hold off;

subplot(3,4,6); hold on;
title('Mach Number');
plot(t(1:nDrogueDeploy),Ma(1:nDrogueDeploy), 'LineWidth', 2);
xlabel('Time[s]'); ylabel('Mach Number');
hold off;

% subplot(3,4,7); hold on;
% title('Stability');
% plot(t(1:nDrogueDeploy),Stability(1:nDrogueDeploy), 'LineWidth', 2);
% plot([t(1) t(nDrogueDeploy)],[R.dn R.dn],'--r');
% plot([t(1) t(nDrogueDeploy)],[2*R.dn 2*R.dn],'--g');
% xlabel('Time[s]'); ylabel('Stability[m]');
% hold off;

subplot(3,4,8); hold on;
title('Angle of Attack');
plot(t(1:nDrogueDeploy),a(1:nDrogueDeploy)*180/pi, 'LineWidth', 2);
xlabel('Time[s]'); ylabel('AoA[deg]');
hold off;

subplot(3,4,10); hold on;
title('Heading');
plot(t(1:nDrogueDeploy),alpha(1:nDrogueDeploy)*180/pi, 'LineWidth', 2);
xlabel('Time[s]'); ylabel('\alpha[deg]');
hold off;

subplot(3,4,11); hold on;
title('Elevation');
plot(t(1:nDrogueDeploy),beta(1:nDrogueDeploy)*180/pi, 'LineWidth', 2);
xlabel('Time[s]'); ylabel('\beta[deg]');
hold off;

subplot(3,4,12); hold on;
title('Tilt');
plot(t(1:nDrogueDeploy),gamma(1:nDrogueDeploy)*180/pi, 'LineWidth', 2);
xlabel('Time[s]'); ylabel('\gamma[deg]');
hold off;

%% Mass Report:
if massReport
    massplots = figure('WindowState','maximized');
    set(gcf,'Color',[1 1 1]);
    set(gca, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');

    subplot(3,3,1); hold on;
    title('Rocket Mass');
    plot([t(1) t(end)],[R.drymass R.drymass],'--g');
    plot([t(1) t(end)],[R.wetmass R.wetmass],'--r');
    plot(t,R.mass(1:end-1));
    legend('Drymass','Wetmass');
    xlabel('Time[s]'); ylabel('Mass[kg]');
    hold off;

    subplot(3,3,2); hold on;
    title('LOX Mass');
    plot(t,R.LOX.mass);
    xlabel('Time[s]'); ylabel('LOX Mass[kg]');
    hold off;

    subplot(3,3,3); hold on;
    title('CH4 Mass');
    plot(t,R.CH4.mass);
    xlabel('Time[s]'); ylabel('CH4 Mass[kg]');
    hold off;

    subplot(3,3,4); hold on;
    title('Rocket Mass Flow');
    plot(t,R.massflow);
    xlabel('Time[s]'); ylabel('Mass Flow Rate[kg/s]');
    hold off;

    subplot(3,3,5); hold on;
    title('LOX Inertia');
    ILOX = cell2mat(R.LOX.I);
    plot(t,ILOX(1,1:3:end));
    plot(t,ILOX(2,2:3:end));
    plot(t,ILOX(3,3:3:end));
    legend('Ix','Iy','Iz');
    xlabel('Time[s]'); ylabel('MoI of LOX[kgm^2]');
    hold off;

    subplot(3,3,6); hold on;
    title('CH4 Inertia');
    ICH4 = cell2mat(R.CH4.I);
    plot(t,ICH4(1,1:3:end));
    plot(t,ICH4(2,2:3:end));
    plot(t,ICH4(3,3:3:end));
    legend('Ix','Iy','Iz');
    xlabel('Time[s]'); ylabel('MoI of CH4[kgm^2]');
    hold off;

    subplot(3,3,7); hold on;
    title('Total Inertia of Rocket');
    I = cell2mat(R.rocketInertia);
    plot(t,I(1,1:3:end));
    plot(t,I(2,2:3:end));
    plot(t,I(3,3:3:end));
    plot([t(1) t(end)],[R.rocketDryInertia(1) R.rocketDryInertia(1)]);
    plot([t(1) t(end)],[R.rocketDryInertia(2) R.rocketDryInertia(2)]);
    plot([t(1) t(end)],[R.rocketDryInertia(3) R.rocketDryInertia(3)]);
    legend('Ix','Iy','Iz','Ixdry','Iydry','Izdry');
    xlabel('Time[s]'); ylabel('MoI of Rocket[kgm^2]');
    hold off;

    subplot(3,3,8); hold on;
    title('Thrust');
    plot(t,R.thrust);
    xlabel('Time[s]'); ylabel('Thrust[N]');
    hold off;

    subplot(3,3,9); hold on;
    title('Center of Masses');
    CM = cell2mat(R.rocketCentroid);
    CMLOX = cell2mat(R.LOX.cent);
    CMCH4 = cell2mat(R.CH4.cent);
    plot(t,CM(3:3:end));
    plot(t,CMLOX(3:3:end));
    plot(t,CMCH4(3:3:end));
    legend('Full Rocket','LOX Tank','CH4 Tank');
    xlabel('Time[s]'); ylabel('Mass[kg]');
    hold off;
end
%% Recovery Plots:
if recoReport
    recoplots = figure('WindowState','maximized');
    set(gcf,'Color',[1 1 1]);
    set(gca, 'FontName', 'TimesNewRoman', 'FontWeight', 'bold');
    sgtitle('Rocket Descent');

    subplot(3,4,[1 5 9]); hold on;
    title('Rocket Trajectory'); axis equal; view(45,25);
    plot3(X(1:nDrogueDeploy),Y(1:nDrogueDeploy),Z(1:nDrogueDeploy), 'LineWidth', 2,'Color',[.25 0 .25]); 
    plot3(X(nDrogueDeploy:nMainDeploy),Y(nDrogueDeploy:nMainDeploy),Z(nDrogueDeploy:nMainDeploy), 'LineWidth', 2,'Color',[.5 0 .5]); 
    plot3(X(nMainDeploy:end),Y(nMainDeploy:end),Z(nMainDeploy:end), 'LineWidth', 2,'Color',[.75 0 .75]);
    legend('Ascent','Drogue Phase','Main Chute Phase');
    plot3(X(nDrogueDeploy),Y(nDrogueDeploy),Z(nDrogueDeploy),'g*'); plot3(X(end),Y(end),Z(end),'r*');
    axis([min(X) max(X) min(Y) max(Y) min(Z) max(Z)]);
    xlabel('Easting[m]'); ylabel('Northing[m]'); zlabel('Altitude[m]');
    hold off;

    subplot(3,4,2); hold on;
    title(sprintf('Vertical Distance(%.1fft Apogee)',(max(Z)-Z(1))*3.28084));
    plot(t(nDrogueDeploy:end),Z(nDrogueDeploy:end)-Z(1), 'LineWidth', 2);
    plot(t(nDrogueDeploy:end),Zdrogue(nDrogueDeploy:end)-Z(1), 'LineWidth', 2);
    plot(t(nMainDeploy:end),Zmain(nMainDeploy:end)-Z(1), 'LineWidth', 2);
    plot(t(find(Z == max(Z))),Z(find(Z==max(Z)))-Z(1),'g*');
    xlabel('Time[s]'); ylabel('Altitude[m]');
    legend('Rocket Height','Drogue Height','Main Chute Height');
    hold off;

    subplot(3,4,3); hold on;
    title('Velocity Magnitude');
    plot(t(nDrogueDeploy:end),Vmag(nDrogueDeploy:end), 'LineWidth', 2);
    plot(t(nDrogueDeploy:end),Vmagdrogue(nDrogueDeploy:end), 'LineWidth', 2);
    plot(t(nMainDeploy:end),Vmagmain(nMainDeploy:end), 'LineWidth', 2);
    xlabel('Time[s]'); ylabel('Velocity[m/s]');
    legend('Rocket','Drogue','Main Chute');
    hold off;

    subplot(3,4,4); hold on;
    title('Acceleration Magnitude');
    plot(t(nDrogueDeploy:end-1),Amag(nDrogueDeploy:end), 'LineWidth', 2);
    plot(t(nDrogueDeploy:end-1),Amagdrogue(nDrogueDeploy:end), 'LineWidth', 2);
    plot(t(nMainDeploy:end-1),Amagmain(nMainDeploy:end), 'LineWidth', 2);
    xlabel('Time[s]'); ylabel('Acceleration[m/s]');
    legend('Rocket','Drogue','Main Chute');
    hold off;

    subplot(3,4,6); hold on;
    title('Line Length');
    plot(t(nDrogueDeploy:end),Ldrogue(nDrogueDeploy:end), 'LineWidth', 2);
    plot(t(nMainDeploy:end),Lmain(nMainDeploy:end), 'LineWidth', 2);
    plot([t(nDrogueDeploy) t(end)],[Reco.Drogue.L Reco.Drogue.L],'--k');
    plot([t(nDrogueDeploy) t(end)],[Reco.Main.L Reco.Main.L],'--k');
    legend('Drogue','Main Chute','Unstretched Rope Lengths');
    xlabel('Time[s]'); ylabel('Distance[m]');
    hold off;

    subplot(3,4,7); hold on;
    title('Rope Tension');
    plot(t(nDrogueDeploy:end),Tdrogue(nDrogueDeploy:end), 'LineWidth', 2);
    plot(t(nMainDeploy:end),Tmain(nMainDeploy:end), 'LineWidth', 2);
    xlabel('Time[s]'); ylabel('Tension[N]');
    hold off;

    subplot(3,4,8); hold on;
    title('Rate of Change of Line Length');
    plot(t(nDrogueDeploy:end),Vldrogue(nDrogueDeploy:end), 'LineWidth', 2);
    plot(t(nMainDeploy:end),Vlmain(nMainDeploy:end), 'LineWidth', 2);
    xlabel('Time[s]'); ylabel('Rate of Change of Line[m/s]');
    legend('Drogue','Main Chute');
    hold off;

    subplot(3,4,10); hold on;
    title('Heading');
    plot(t(nDrogueDeploy:end),alpha(nDrogueDeploy:end)*180/pi, 'LineWidth', 2);
    xlabel('Time[s]'); ylabel('\alpha[deg]');
    hold off;

    subplot(3,4,11); hold on;
    title('Elevation');
    plot(t(nDrogueDeploy:end),beta(nDrogueDeploy:end)*180/pi, 'LineWidth', 2);
    xlabel('Time[s]'); ylabel('\beta[deg]');
    hold off;

    subplot(3,4,12); hold on;
    title('Tilt');
    plot(t(nDrogueDeploy:end),gamma(nDrogueDeploy:end)*180/pi, 'LineWidth', 2);
    xlabel('Time[s]'); ylabel('\gamma[deg]');
    hold off;
end

end