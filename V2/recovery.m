classdef recovery %Recovery System Object
%% Objective: Store dimensions of rocket, calculate variable Inertia tensor, variable CoM. Display functions

%% Constructor Input Variables: 
% Rope:
%   L - rope length
%   k - rope spring constant
%   
% Drogue Chute:
%   Cdc - Coefficient of Drag of Drogue Chute
%   Adc - Area of Drogue Chute
%   mdc - Mass of Drogue Chute
%
% Main Chute:
%   Cmc - Coefficient of Drag of Main Chute
%   Amc - Area of Main Chute
%   mac - Mass of Main Chute
%
% Other:
%   dt - sample rate of simulation(Hz)
%
%% Class Functions:
% 
%% Output Variables: 
%   none
%
%% Functions Called: 
%   
% 
%% Assumptions: 
% massless rope, no angle of attack on parachute, hookes law for rope
% 
%% Model:
% 
%
%% Author: Ethan Foss
% erfoss@ucsd.edu, ethanrfoss@gmail.com
%
%% To Be Implemented:
% 
    properties
        
  
        %Drogue Chute:
        Drogue = struct('L',[],'k',[],'c',[],'Cd',[],'D',[],'A',[],'m',[],'X',[],'Y',[],'Z',[],'hull',[]);
        
        %Main Chute:
        Main = struct('L',[],'k',[],'c',[],'Cd',[],'D',[],'A',[],'m',[],'X',[],'Y',[],'Z',[],'hull',[]);
        
        %Main Deploy Height:
        MainDeployHeight
    end
    
    methods
        function obj = recovery(Ld,kd,cdc,Cddc,Ddc,mdc,Lm,km,cmc,Cdmc,Dmc,mmc,MainDeployHeight)
            
            % Drogue Chute
            obj.Drogue.L = Ld;
            obj.Drogue.k = kd;
            obj.Drogue.c = cdc;
            obj.Drogue.Cd = Cddc;
            obj.Drogue.D = Ddc;
            obj.Drogue.A = pi/4*Ddc^2;
            obj.Drogue.m = mdc;
            %Drogue Coordinates and hull:
            obj.Drogue.X = [Ddc/2*cos(0:pi/16:2*pi)*cos(0) Ddc/2*cos(0:pi/16:2*pi)*cos(pi/8) Ddc/2*cos(0:pi/16:2*pi)*cos(pi/4) Ddc/2*cos(0:pi/16:2*pi)*cos(3*pi/8) 0];
            obj.Drogue.Y = [Ddc/2*sin(0:pi/16:2*pi)*cos(0) Ddc/2*sin(0:pi/16:2*pi)*cos(pi/8) Ddc/2*sin(0:pi/16:2*pi)*cos(pi/4) Ddc/2*sin(0:pi/16:2*pi)*cos(3*pi/8) 0];
            obj.Drogue.Z = [Ddc/2*sin(0)*ones(1,33) Ddc/2*sin(pi/8)*ones(1,33) Ddc/2*sin(pi/4)*ones(1,33) Ddc/2*sin(3*pi/8)*ones(1,33) Ddc/2];
            obj.Drogue.hull = convhull(obj.Drogue.X,obj.Drogue.Y,obj.Drogue.Z);
            
            % Main Chute
            obj.Main.L = Lm;
            obj.Main.k = km;
            obj.Main.c = cmc;
            obj.Main.Cd = Cdmc;
            obj.Main.D = Dmc;
            obj.Main.A = pi/4*Dmc^2;
            obj.Main.m = mmc;
            %Drogue Coordinates and hull:
            obj.Main.X = [Dmc/2*cos(0:pi/16:2*pi)*cos(0) Dmc/2*cos(0:pi/16:2*pi)*cos(pi/8) Dmc/2*cos(0:pi/16:2*pi)*cos(pi/4) Dmc/2*cos(0:pi/16:2*pi)*cos(3*pi/8) 0];
            obj.Main.Y = [Dmc/2*sin(0:pi/16:2*pi)*cos(0) Dmc/2*sin(0:pi/16:2*pi)*cos(pi/8) Dmc/2*sin(0:pi/16:2*pi)*cos(pi/4) Dmc/2*sin(0:pi/16:2*pi)*cos(3*pi/8) 0];
            obj.Main.Z = [Dmc/2*sin(0)*ones(1,33) Dmc/2*sin(pi/8)*ones(1,33) Dmc/2*sin(pi/4)*ones(1,33) Dmc/2*sin(3*pi/8)*ones(1,33) Ddc/2];
            obj.Main.hull = convhull(obj.Main.X,obj.Main.Y,obj.Main.Z);
            
            %Main Deploy Height
            obj.MainDeployHeight = MainDeployHeight;
            
        end
        
        function [Xcat,Ycat,Zcat] = ropeCurve(obj,Xpar,Ypar,Zpar,Xroc,Yroc,Zroc,L)
            
            S = sqrt((Xpar-Xroc)^2+(Ypar-Yroc)^2+(Zpar-Zroc)^2);
            if S>=L
                Xcat = linspace(Xroc,Xpar,100);
                Ycat = linspace(Yroc,Ypar,100);
                Zcat = linspace(Zroc,Zpar,100);
            else
                H = sqrt((Xpar-Xroc)^2+(Ypar-Yroc)^2);
                [X,Y] = catenary([0 Zroc],[H Zpar],50,100);
                Xcat = X*(Xpar-Xroc)/H + Xroc;
                Ycat = X*(Ypar-Yroc)/H +Yroc;
                Zcat = Y;
            end
            
        end
    end
    
end