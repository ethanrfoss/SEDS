classdef rocket %Rocket Object
%% Objective: Store dimensions of rocket, calculate variable Inertia tensor, variable CoM. Display functions

%% Constructor Input Variables: 
% Characteristic lengths:
%   ln - Length of nose cone(m)
%   dn - Diameter of Nose Cone(m)
%   db - Diameter of Rocket body(m) (usually equal to dn)
%   df - Diameter of Rocket at fins(m) (usually equal to dn)
%   du - Diameter of Rocket at end(m) (usually equal to dn)
%   dd - Diameter of Rocket at engine(m)
%   lb - Length of Rocket body(m)
%   lc - Length to engine from bottom of rocket(m)
%   Xf - Length from tip to fins(m)
%   lr - Length of fin at body(m)
%   lt - Length of fin at tip(m)
%   ls - Normal Width of Fin(m)
%   lw - Vertical Length from Start of Fin to Edge of Fin(m)
%   Tf - Fin thickness(m)
%         
% Tanks:
%   lLOXCM - Length to center of mass of LOX tank
%   lLOX - Length of LOX tank
%   dLOX - Diameter of LOX tank
%   m0LOX - Initial Mass of LOX
%         
%   lCH4CM - Length to center of mass of LCH4 tank
%   lCH4 - Length of LCH4 tank
%   dCH4 - Diameter of LCH4 tank
%   m0CH4 - Initial Mass of LCH4
%         
% Burntime:        
%   burntime - burntime(s)(a negative value will prompt a thrust curve)
%         
% Mass:
%   drymass - drymass(kg)
%
% Other:
%   dt - sample rate of simulation(Hz)
%
%% Class Functions:
% displayRocketDimensions(R)
%   Displays Rocket Dimensions
% drawRocket3D
%   Displays 3D Rocket
%% Output Variables: 
%   none
%
%% Functions Called: 
%    drawThrust - Prompts User Input of thrust curve and returns thrust
%    curve.
%    RigidBodyParams - Returns inertia, CoM, and volume of convex hull shape.
% 
%% Assumptions: 
% No fluid sloshing, mass flow rate proportional to thrust, pressurant mass
% invariant, 4 fins, uniform density throughout rocket.
% 
%% Model:
% Basic Inertia Calculations, Center of Mass Calculations
%
%% Author: Ethan Foss
% erfoss@ucsd.edu, ethanrfoss@gmail.com
%
%% To Be Implemented:
% Robust Propulsion Model
%
    properties (Constant)
       
        maximpulse = 40959; %Ns
        defaultthrust = 8896.443; %N
        defaultburntime = 11; %s
        defaultmassflowrate = 2.96783; %kg/s
        n = 4; %Number of fins

    end
    properties
        
        % Characteristic lengths
        ln  %Length of nose cone(m)
        dn  %Diameter of Nose Cone(m)
        db  %Diameter of Rocket body(m) (usually equal to dn)
        df  %Diameter of Rocket at fins(m) (usually equal to dn)
        du  %Diameter of Rocket at end(m) (usually equal to dn)
        dd  %Diameter of Rocket at engine(m) (usually equal to dn)
        lb  %Length of Rocket body(m)
        lc  %Length to engine from bottom of rocket(m)
        Xb  %Length from tip to nose cone(m)
        Xf  %Length from tip to fins(m)
        Xc  %Length from tip to rocket end(m)
        lTR  %Total rocket length(m)
        lr  %Length of fin at body(m)
        lt  %Length of fin at tip(m)
        ls  %Normal Width of Fin(m)
        lw  %Vertical Length from Start of Fin to Edge of Fin(m)
        lm  %Width of Fin(m)
        lTS  %Width of Rocket(m)
        Tf %Fin thickness
        Aref %Reference Area
        
        % Tanks
        lLOXCM %Length to center of mass of LOX tank
        lLOX %Length of LOX tank
        dLOX %Diameter of LOX tank
        m0LOX %Initial Mass of LOX
        LOXvol %Volume of LOX tank
        
        lCH4CM %Length to center of mass of LCH4 tank
        lCH4 %Length of LCH4 tank
        dCH4 %Diameter of LCH4 tank
        m0CH4 %Initial Mass of LCH4
        CH4vol %Volume of LCH4 tank
        
        % Burntime, Thrust
        
        burntime %burntime(s)
        thrust %thrust(N)
        
        % Masses, mass flow
        massflow %massflow(kg/s)
        drymass %drymass(kg)
        wetmass %wetmass(kg)
        mass %(kg)
        
        %Inertias, Centroids, Volumes
        rocketDryCentroid %3d vector of centroid(m) (origin at tip of nose cone)
        rocketCentroid %Struct containg Centroid values at each time step(m)
        rocketDryInertia %3x3 Matrix describing Inertia Tensor With no Fuel(kgm^2)
        rocketInertia %Struct containing Inertia Tensor at each time step(kgm^2)
        rocketVolume %(m^3)
        rocketDryDensity %(kg/m^3)
        
        %Structures of Objects(X dimensions, Y dimensions, Z dimensions,
        %convhull, volume, centroid, inertia tensor)
        
        cone = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[]);
        body = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[]);
        finR = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[]);
        finL = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[]);
        finF = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[]);
        finB = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[]);
        skirt = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[]);
        LOX = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[],'mass',[]);
        CH4 = struct('X',[],'Y',[],'Z',[],'hull',[],'polyshape',[],'vol',[],'cent',[],'I',[],'mass',[]);
        
        % Other
        dt %Sample rate of simulation(Hz)
        Acp %Planform Area
        Xcpb90 %Body Center of pressure at 90 AOA
        Xcp90 %Total Rocket Center of Pressure at 90 AOA
    end
    
    methods
        function obj = rocket(ln,dn,db,df,du,dd,lb,lc,Xf,lr,lt,ls,lw,Tf,lLOXCM,lLOX,dLOX,m0LOX,lCH4CM,lCH4,dCH4,m0CH4,burntime,drymass,dt)
            disp('Calculating Rocket Parameters');
            % Rocket 
            obj.dt = dt;
            obj.ln = ln;
            obj.dn = dn;
            obj.db = db;
            obj.df = df;
            obj.du = du;
            obj.dd = dd;
            obj.lb = lb;
            obj.lc = lc;
            obj.Xb = ln;
            obj.Xf = Xf;
            obj.Xc = ln+lb;
            obj.lTR = ln+lb+lc;
            obj.lr = lr;
            obj.lt = lt;
            obj.ls = ls;
            obj.lw = lw;
            obj.lm = sqrt((lt/2+lw-lr/2)^2+ls^2);
            obj.lTS = 2*ls+du;
            obj.Tf = Tf;
            obj.Aref = pi*dn^2/4;
            
            % Tanks
            obj.lLOXCM = lLOXCM;
            obj.lLOX = lLOX;
            obj.dLOX = dLOX; 
            obj.m0LOX = m0LOX; 
            obj.LOXvol = lLOX*pi*dLOX^2/4;
            
            obj.lCH4CM = lCH4CM;
            obj.lCH4 = lCH4;
            obj.dCH4 = dCH4; 
            obj.m0CH4 = m0CH4; 
            obj.CH4vol = lCH4*pi*dCH4^2/4;
            
            % Mass Flow
            obj.drymass = drymass;
            obj.wetmass = obj.drymass + obj.m0CH4 + obj.m0LOX;
            if burntime>0
                % Positive Burntime Case, Constant Thrust:
                obj.burntime = burntime;
                obj.thrust = obj.maximpulse/burntime*ones(1,length(0:dt:burntime));
                % Mass Flow is Constant:
                obj.massflow = (m0CH4+m0LOX)/(obj.burntime)*ones(1,length(0:dt:burntime));
                
                % Mass:
                obj.mass(1) = obj.wetmass;
                obj.mass = obj.mass(1) - obj.massflow(1)*([0:obj.dt:obj.burntime]);
            else
                % Negative Burntime, Request User Thrust Curve:
                obj.burntime = -1;
                [~,obj.thrust] = drawThrust;
                % Massflow proportional to thrust:
                obj.massflow = obj.thrust*(obj.m0CH4+obj.m0LOX)/(obj.maximpulse);
                
                % Mass
                obj.mass  = zeros(1,length(obj.massflow));
                obj.mass(1) = obj.wetmass;
                % Integrate massflow for mass:
                for i = 1:length(obj.massflow)
                    obj.mass(i+1) = obj.mass(i) - obj.massflow(i)*obj.dt;
                end
            end
            
            
            
            % Construction of Rocket Shapes. The Rocket is split into 7
            % shapes: Cone, Body, Skirt, and Fins. Vertexes of these
            % objects are created based on constructor inputs. A convex
            % hull is created to represent these shapes, along with a
            % polyshape, which represents the 2D form of these shapes.
            % Using the RigidBodyParams method, The volume, Centroid, and
            % Inertia Tensor can be found. Constant density of the shapes
            % is assumed
            
            warning('off'); %Turn off warnings(polyshape gives a lot of warnings)
            
            
            %Cone
            obj.cone.X = [obj.dn/2*cos(0:pi/16:2*pi) 0];
            obj.cone.Y = [obj.dn/2*sin(0:pi/16:2*pi) 0];
            obj.cone.Z = [-obj.Xb*ones(1,length(0:pi/16:2*pi)) 0];
            obj.cone.hull = convhull(obj.cone.X,obj.cone.Y,obj.cone.Z);
            obj.cone.polyshape = polyshape([0 obj.dn/2 -obj.dn/2],[0 -obj.ln -obj.ln]);
            [RBPcone,~] = RigidBodyParams(triangulation(obj.cone.hull,obj.cone.X',obj.cone.Y',obj.cone.Z'));
            obj.cone.vol = RBPcone.volume; obj.cone.cent = RBPcone.centroid; obj.cone.I = RBPcone.inertia_tensor;
            
            %Body
            obj.body.X = [obj.dn/2*cos(0:pi/16:2*pi) obj.df/2*cos(0:pi/16:2*pi) obj.du/2*cos(0:pi/16:2*pi)];
            obj.body.Y = [obj.dn/2*sin(0:pi/16:2*pi) obj.df/2*sin(0:pi/16:2*pi) obj.du/2*sin(0:pi/16:2*pi)];
            obj.body.Z = [-obj.Xb*ones(1,length(0:pi/16:2*pi)) -obj.Xf*ones(1,length(0:pi/16:2*pi)) -obj.Xc*ones(1,length(0:pi/16:2*pi))];
            obj.body.hull = convhull(obj.body.X,obj.body.Y,obj.body.Z);
            obj.body.polyshape = polyshape([obj.dn/2 obj.df/2 obj.du/2 -obj.du/2 -obj.df/2 -obj.dn/2],[-obj.Xb -obj.Xf -obj.Xc -obj.Xc -obj.Xf -obj.Xb]);
            [RBPbody,~] = RigidBodyParams(triangulation(obj.body.hull,obj.body.X',obj.body.Y',obj.body.Z'));
            obj.body.vol = RBPbody.volume; obj.body.cent = RBPbody.centroid; obj.body.I = RBPbody.inertia_tensor;
            
            %Right Fin
            obj.finR.X = [obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2 obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2];
            obj.finR.Y = [obj.Tf/2 obj.Tf/2 obj.Tf/2 obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2];
            obj.finR.Z = [-obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr -obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr];
            obj.finR.hull = convhull(obj.finR.X,obj.finR.Y,obj.finR.Z);
            obj.finR.polyshape = polyshape(obj.finR.X(1:4),obj.finR.Z(1:4));
            [RBPfinR,~] = RigidBodyParams(triangulation(obj.finR.hull,obj.finR.X',obj.finR.Y',obj.finR.Z'));
            obj.finR.vol = RBPfinR.volume; obj.finR.cent = RBPfinR.centroid; obj.finR.I = RBPfinR.inertia_tensor;
            
            %Left Fin
            obj.finL.X = -[obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2 obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2];
            obj.finL.Y = [obj.Tf/2 obj.Tf/2 obj.Tf/2 obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2];
            obj.finL.Z = [-obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr -obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr];
            obj.finL.hull = convhull(obj.finL.X,obj.finL.Y,obj.finL.Z);
            obj.finL.polyshape = polyshape(obj.finL.X(1:4),obj.finL.Z(1:4));
            [RBPfinL,~] = RigidBodyParams(triangulation(obj.finL.hull,obj.finL.X',obj.finL.Y',obj.finL.Z'));
            obj.finL.vol = RBPfinL.volume; obj.finL.cent = RBPfinL.centroid; obj.finL.I = RBPfinL.inertia_tensor;
            
            %Front Fin
            obj.finF.X = [obj.Tf/2 obj.Tf/2 obj.Tf/2 obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2];
            obj.finF.Y = [obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2 obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2];
            obj.finF.Z = [-obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr -obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr];
            obj.finF.hull = convhull(obj.finF.X,obj.finF.Y,obj.finF.Z);
            obj.finF.polyshape = polyshape(obj.finF.X,[obj.finF.Z(1:4) obj.finF.Z(4:-1:1)]);
            [RBPfinF,~] = RigidBodyParams(triangulation(obj.finF.hull,obj.finF.X',obj.finF.Y',obj.finF.Z'));
            obj.finF.vol = RBPfinF.volume; obj.finF.cent = RBPfinF.centroid; obj.finF.I = RBPfinF.inertia_tensor;
            
            %Back Fin
            obj.finB.X = [obj.Tf/2 obj.Tf/2 obj.Tf/2 obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2 -obj.Tf/2];
            obj.finB.Y = -[obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2 obj.df/2 obj.lTS/2 obj.lTS/2 obj.du/2];
            obj.finB.Z = [-obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr -obj.Xf, -obj.Xf-obj.lw, -obj.Xf-obj.lw-obj.lt, -obj.Xf-obj.lr];
            obj.finB.hull = convhull(obj.finB.X,obj.finB.Y,obj.finB.Z);
            obj.finB.polyshape = polyshape(obj.finB.X,[obj.finF.Z(1:4) obj.finB.Z(4:-1:1)]);
            [RBPfinB,~] = RigidBodyParams(triangulation(obj.finB.hull,obj.finB.X',obj.finB.Y',obj.finB.Z'));
            obj.finB.vol = RBPfinB.volume; obj.finB.cent = RBPfinB.centroid; obj.finB.I = RBPfinB.inertia_tensor;
            
            %Skirt
            obj.skirt.X = [obj.du/2*cos(0:pi/16:2*pi) obj.dd/2*cos(0:pi/16:2*pi)];
            obj.skirt.Y = [obj.du/2*sin(0:pi/16:2*pi) obj.dd/2*sin(0:pi/16:2*pi)];
            obj.skirt.Z = [-obj.Xc*ones(1,length(0:pi/16:2*pi)) -obj.lTR*ones(1,length(0:pi/16:2*pi))];
            obj.skirt.hull = convhull(obj.skirt.X,obj.skirt.Y,obj.skirt.Z);
            obj.skirt.polyshape = polyshape([-obj.du/2 obj.du/2 obj.dd/2 -obj.dd/2],[-obj.Xc -obj.Xc -obj.lTR -obj.lTR]);
            [RBPskirt,~] = RigidBodyParams(triangulation(obj.skirt.hull,obj.skirt.X',obj.skirt.Y',obj.skirt.Z'));
            obj.skirt.vol = RBPskirt.volume; obj.skirt.cent = RBPskirt.centroid; obj.skirt.I = RBPskirt.inertia_tensor;
            
            % LOX:
            obj.LOX.X = [obj.dLOX/2*cos(0:pi/16:2*pi) obj.dLOX/2*cos(0:pi/16:2*pi)];
            obj.LOX.Y = [obj.dLOX/2*sin(0:pi/16:2*pi) obj.dLOX/2*sin(0:pi/16:2*pi)];
            obj.LOX.Z = [-(obj.lLOXCM+obj.lLOX/2)*ones(1,length(0:pi/16:2*pi)) -(obj.lLOXCM-obj.lLOX/2)*ones(1,length(0:pi/16:2*pi))];
            obj.LOX.hull = convhull(obj.LOX.X,obj.LOX.Y,obj.LOX.Z);
            obj.LOX.polyshape = polyshape([-obj.dLOX/2 obj.dLOX/2 obj.dLOX/2 -obj.dLOX/2],[-(obj.lLOXCM+obj.lLOX/2) -(obj.lLOXCM+obj.lLOX/2) -(obj.lLOXCM-obj.lLOX/2) -(obj.lLOXCM-obj.lLOX/2)]);
            
            % Mass, Volume, Centroid, Inertia of LOX(changing with time):
            obj.LOX.mass(1) = m0LOX;
            obj.LOX.vol(1) = lLOX*pi/4*dLOX^2;
            obj.LOX.cent{1} = [0,0,-lLOXCM];
            obj.LOX.I{1} = obj.LOX.mass(1)*diag([lLOX^2/12+dLOX^2/16, lLOX^2/12+dLOX^2/16, dLOX^2/16],0);
            for i = 1:length(obj.mass)-1
                obj.LOX.mass(i+1) = obj.LOX.mass(i) - obj.massflow(i)*obj.m0LOX/(obj.m0CH4+obj.m0LOX)*dt;
                obj.LOX.vol(i+1) = obj.LOX.vol(1)/obj.LOX.mass(1)*obj.LOX.mass(i+1);
                l = obj.LOX.vol(i+1)/(pi/4*dLOX^2);
                obj.LOX.cent{i+1} = [0,0,-(lLOXCM+lLOX/2)+l/2];
                obj.LOX.I{i+1} = obj.LOX.mass(i+1)*diag([l^2/12+dLOX^2/16, l^2/12+dLOX^2/16, dLOX^2/16],0);
            end
            
            % CH4:
            obj.CH4.X = [obj.dCH4/2*cos(0:pi/16:2*pi) obj.dCH4/2*cos(0:pi/16:2*pi)];
            obj.CH4.Y = [obj.dCH4/2*sin(0:pi/16:2*pi) obj.dCH4/2*sin(0:pi/16:2*pi)];
            obj.CH4.Z = [-(obj.lCH4CM+obj.lCH4/2)*ones(1,length(0:pi/16:2*pi)) -(obj.lCH4CM-obj.lCH4/2)*ones(1,length(0:pi/16:2*pi))];
            obj.CH4.hull = convhull(obj.CH4.X,obj.CH4.Y,obj.CH4.Z);
            obj.CH4.polyshape = polyshape([-obj.dCH4/2 obj.dCH4/2 obj.dCH4/2 -obj.dCH4/2],[-(obj.lCH4CM+obj.lCH4/2) -(obj.lCH4CM+obj.lCH4/2) -(obj.lCH4CM-obj.lCH4/2) -(obj.lCH4CM-obj.lCH4/2)]);

            % Mass, Volume, Centroid, Inertia of CH4(changing with time):
            obj.CH4.mass(1) = m0CH4;
            obj.CH4.vol(1) = lCH4*pi/4*dCH4^2;
            obj.CH4.cent{1} = [0,0,-lCH4CM];
            obj.CH4.I{1} = obj.CH4.mass(1)*diag([lCH4^2/12+dCH4^2/16, lCH4^2/12+dCH4^2/16, dCH4^2/16],0);
            for i = 1:length(obj.mass)-1
                obj.CH4.mass(i+1) = obj.CH4.mass(i) - obj.massflow(i)*obj.m0CH4/(obj.m0CH4+obj.m0LOX)*dt; 
                obj.CH4.vol(i+1) = obj.CH4.vol(1)/obj.CH4.mass(1)*obj.CH4.mass(i+1);
                l = obj.CH4.vol(i+1)/(pi/4*dCH4^2);
                obj.CH4.cent{i+1} = [0,0,-(lCH4CM+lCH4/2)+l/2];
                obj.CH4.I{i+1} = obj.CH4.mass(i+1)*diag([l^2/12+dCH4^2/16, l^2/12+dCH4^2/16, dCH4^2/16],0);
            end
            
            %Rocket
            obj.rocketVolume = obj.cone.vol+obj.body.vol+obj.finR.vol+obj.finL.vol+obj.finF.vol+obj.finB.vol+obj.skirt.vol;
            obj.rocketDryDensity = obj.drymass/obj.rocketVolume;
            ICM = @(disp,V) V*[disp(2)^2+disp(3)^2,-disp(1)*disp(2),-disp(1)*disp(3);-disp(1)*disp(2),disp(1)^2+disp(3)^2,-disp(2)*disp(3);-disp(1)*disp(3),-disp(2)*disp(3),disp(1)^2+disp(2)^2]; %Translation matrix
            obj.rocketDryCentroid = (obj.cone.cent*obj.cone.vol+obj.body.cent*obj.body.vol+obj.finR.cent*obj.finR.vol+obj.finL.cent*obj.finL.vol+obj.finF.cent*obj.finF.vol+obj.finB.cent*obj.finB.vol+obj.skirt.cent*obj.skirt.vol)/(obj.rocketVolume);
            obj.rocketDryInertia = obj.rocketDryDensity*(obj.cone.I+ICM(obj.rocketDryCentroid-obj.cone.cent,obj.cone.vol) + obj.body.I+ICM(obj.rocketDryCentroid-obj.body.cent,obj.body.vol) + obj.finR.I+ICM(obj.rocketDryCentroid-obj.finR.cent,obj.finR.vol) + obj.finL.I+ICM(obj.rocketDryCentroid-obj.finL.cent,obj.finL.vol) + obj.finB.I+ICM(obj.rocketDryCentroid-obj.finB.cent,obj.finB.vol) + obj.finF.I+ICM(obj.rocketDryCentroid-obj.finF.cent,obj.finF.vol) + obj.skirt.I+ICM(obj.rocketDryCentroid-obj.skirt.cent,obj.skirt.vol));
            for i = 1:length(obj.mass)
                % Centroid and Volume
                obj.rocketCentroid{i} = (obj.drymass*obj.rocketDryCentroid+obj.CH4.cent{i}*obj.CH4.mass(i) + obj.LOX.cent{i}*obj.LOX.mass(i))/(obj.mass(i));
                
                %Inertia of Rocket
                obj.rocketInertia{i} = obj.rocketDryDensity*(obj.cone.I+ICM(obj.rocketCentroid{i}-obj.cone.cent,obj.cone.vol) + obj.body.I+ICM(obj.rocketCentroid{i}-obj.body.cent,obj.body.vol) + obj.finR.I+ICM(obj.rocketCentroid{i}-obj.finR.cent,obj.finR.vol) + obj.finL.I+ICM(obj.rocketCentroid{i}-obj.finL.cent,obj.finL.vol) + obj.finB.I+ICM(obj.rocketCentroid{i}-obj.finB.cent,obj.finB.vol) + obj.finF.I+ICM(obj.rocketCentroid{i}-obj.finF.cent,obj.finF.vol) + obj.skirt.I+ICM(obj.rocketCentroid{i}-obj.skirt.cent,obj.skirt.vol)) + obj.LOX.I{i} + obj.LOX.mass(i)*ICM(obj.rocketCentroid{i}-obj.LOX.cent{i},1) + obj.CH4.I{i} + obj.CH4.mass(i)*ICM(obj.rocketCentroid{i}-obj.CH4.cent{i},1);
            end
            
            % Planform Area:
            [ConeGeom,~,~] = polygeom(obj.cone.polyshape.Vertices(:,1),obj.cone.polyshape.Vertices(:,2));
            [BodyGeom,~,~] = polygeom(obj.body.polyshape.Vertices(:,1),obj.body.polyshape.Vertices(:,2));
            [SkirtGeom,~,~] = polygeom(obj.skirt.polyshape.Vertices(:,1),obj.skirt.polyshape.Vertices(:,2));
            [FinRGeom,~,~] = polygeom(obj.finR.polyshape.Vertices(:,1),obj.finR.polyshape.Vertices(:,2));
            [FinLGeom,~,~] = polygeom(obj.finL.polyshape.Vertices(:,1),obj.finL.polyshape.Vertices(:,2));
            obj.Xcpb90 = -(ConeGeom(1)*ConeGeom(3)+BodyGeom(1)*BodyGeom(3)+SkirtGeom(1)*SkirtGeom(3))/(ConeGeom(1)+BodyGeom(1)+SkirtGeom(1));
            obj.Xcp90 = -(ConeGeom(1)*ConeGeom(3)+BodyGeom(1)*BodyGeom(3)+SkirtGeom(1)*SkirtGeom(3)+FinRGeom(1)*FinRGeom(3)+FinLGeom(1)*FinRGeom(3))/(ConeGeom(1)+BodyGeom(1)+SkirtGeom(1)+FinRGeom(1)+FinLGeom(1));
            obj.Acp = ConeGeom(1)+BodyGeom(1)+SkirtGeom(1)+FinRGeom(1)+FinLGeom(1);
            
            disp('Done Calculating Rocket Parameters');
        end
    end
    
    methods (Static)
        function displayRocketDimensions(obj)
            figure; hold on;
            
            % Poltshapes:
            plot(obj.cone.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(obj.body.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(obj.skirt.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(obj.LOX.polyshape,'FaceColor','g');
            plot(obj.CH4.polyshape,'FaceColor','g');
            plot(obj.finR.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(obj.finL.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(obj.finF.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(obj.finB.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(obj.rocketCentroid{1}(1),obj.rocketCentroid{1}(3),'r*');
            
            % Lines:
            plot([-1.2*obj.lTS/2 1.2*obj.lTS/2],[0 0],'k--');
            plot([-obj.dn/2 3*obj.dn/2],[-obj.Xb -obj.Xb],'k--');
            plot([obj.df/2 obj.lTS/2],[-obj.Xf -obj.Xf],'k--');
            plot([obj.du/2 1.3*obj.lTS/2],[-obj.Xc -obj.Xc],'k--');
            plot([obj.dd/2 1.5*obj.lTS/2],[-obj.lTR -obj.lTR],'k--');
            
            % Horizontal Quivers:
            qdn1 = quiver(0,-obj.ln,obj.dn/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.dn);
            qdn2 = quiver(0,-obj.ln,-obj.dn/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.dn);
            label(qdn1,sprintf('d_{n}=%1.2fm',obj.dn));
            qdf1 = quiver(0,-obj.Xf,obj.df/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.df);
            qdf2 = quiver(0,-obj.Xf,-obj.df/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.df);
            label(qdf1,sprintf('d_{f}=%1.2fm',obj.df));
            qdu1 = quiver(0,-obj.Xc,obj.du/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.du);
            qdu2 = quiver(0,-obj.Xc,-obj.du/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.du);
            label(qdu1,sprintf('d_{u}=%1.2fm',obj.du));
            qdf1 = quiver(0,-obj.lTR,obj.dd/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.dd);
            qdf2 = quiver(0,-obj.lTR,-obj.dd/2,0,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',1/obj.dd);
            label(qdf1,sprintf('d_{d}=%1.2fm',obj.df));
            
            % Vertical Quivers:
            qln = quiver(obj.dn/2,0,0,-obj.ln,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',.5/obj.ln);
            label(qln,sprintf('l_{n}=%1.2fm',obj.ln),'rotation',-90);
            qXb = quiver(2.5*obj.dn/2,0,0,-obj.Xb,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',.5/obj.Xf);
            label(qXb,sprintf('X_{b}=%1.2fm',obj.Xb),'rotation',-90);
            qXf = quiver(4*obj.dn/2,0,0,-obj.Xf,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',.5/obj.Xf);
            label(qXf,sprintf('X_{f}=%1.2fm',obj.Xf),'rotation',-90);
            qXc = quiver(1.2*obj.lTS/2,0,0,-obj.Xc,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',.5/obj.Xc);
            label(qXc,sprintf('X_{c}=%1.2fm',obj.Xc),'rotation',-90);
            qlTR = quiver(1.6*obj.lTS/2,0,0,-obj.lTR,'k','AutoScale','off','LineWidth',3,'MaxHeadSize',.5/obj.lTR);
            label(qlTR,sprintf('l_{TR}=%1.2fm',obj.lTR),'rotation',-90);
            axis equal;
            axis([-1.6*obj.lTS/2 1.6*obj.lTS/2 -1.1*obj.lTR .1]);
            drawnow;
        end
        
        function drawRocket3D(obj)
            figure; hold on;
            trisurf(obj.cone.hull,obj.cone.X,obj.cone.Y,obj.cone.Z,'FaceColor',[0 0.4470 0.7410]);
            trisurf(obj.body.hull,obj.body.X,obj.body.Y,obj.body.Z,'FaceColor',[0 0.4470 0.7410]);
            trisurf(obj.finR.hull,obj.finR.X,obj.finR.Y,obj.finR.Z,'FaceColor',[0 0.4470 0.7410]);
            trisurf(obj.finL.hull,obj.finL.X,obj.finL.Y,obj.finL.Z,'FaceColor',[0 0.4470 0.7410]);
            trisurf(obj.finF.hull,obj.finF.X,obj.finF.Y,obj.finF.Z,'FaceColor',[0 0.4470 0.7410]);
            trisurf(obj.finB.hull,obj.finB.X,obj.finB.Y,obj.finB.Z,'FaceColor',[0 0.4470 0.7410]);
            trisurf(obj.skirt.hull,obj.skirt.X,obj.skirt.Y,obj.skirt.Z,'FaceColor',[0 0.4470 0.7410]);
            axis equal;
            view(135,45);
        end
    end
end