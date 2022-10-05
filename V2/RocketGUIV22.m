classdef RocketGUIV22 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        MassParametersPanel             matlab.ui.container.Panel
        LOXCOMEditFieldLabel            matlab.ui.control.Label
        LOXCOMEditField                 matlab.ui.control.NumericEditField
        LOXTankLengthEditFieldLabel     matlab.ui.control.Label
        LOXTankLengthEditField          matlab.ui.control.NumericEditField
        LOXTankDiamEditFieldLabel       matlab.ui.control.Label
        LOXTankDiamEditField            matlab.ui.control.NumericEditField
        LOXMassEditFieldLabel           matlab.ui.control.Label
        LOXMassEditField                matlab.ui.control.NumericEditField
        CH4COMEditFieldLabel            matlab.ui.control.Label
        CH4COMEditField                 matlab.ui.control.NumericEditField
        CH4TankLengthEditFieldLabel     matlab.ui.control.Label
        CH4TankLengthEditField          matlab.ui.control.NumericEditField
        CH4TankDiamEditFieldLabel       matlab.ui.control.Label
        CH4TankDiamEditField            matlab.ui.control.NumericEditField
        CH4MassEditFieldLabel           matlab.ui.control.Label
        CH4MassEditField                matlab.ui.control.NumericEditField
        DrymassEditFieldLabel           matlab.ui.control.Label
        DrymassEditField                matlab.ui.control.NumericEditField
        GenerateCustomThrustCurveCheckBox  matlab.ui.control.CheckBox
        BurntimeEditFieldLabel          matlab.ui.control.Label
        BurntimeEditField               matlab.ui.control.NumericEditField
        ThrustEditFieldLabel            matlab.ui.control.Label
        ThrustEditField                 matlab.ui.control.NumericEditField
        UpdateMassParametersButton      matlab.ui.control.Button
        LoadRocketPanel                 matlab.ui.container.Panel
        SavedRocketsListBoxLabel        matlab.ui.control.Label
        SavedRocketsListBox             matlab.ui.control.ListBox
        LoadRocketButton                matlab.ui.control.Button
        RocketLengthsPanel              matlab.ui.container.Panel
        NoseLengthEditFieldLabel        matlab.ui.control.Label
        NoseLengthEditField             matlab.ui.control.NumericEditField
        NoseDiamEditFieldLabel          matlab.ui.control.Label
        NoseDiamEditField               matlab.ui.control.NumericEditField
        BodyDiamEditFieldLabel          matlab.ui.control.Label
        BodyDiamEditField               matlab.ui.control.NumericEditField
        DiamatFinsEditFieldLabel        matlab.ui.control.Label
        DiamatFinsEditField             matlab.ui.control.NumericEditField
        EndBodyDiamEditFieldLabel       matlab.ui.control.Label
        EndBodyDiamEditField            matlab.ui.control.NumericEditField
        SkirtDiamEditFieldLabel         matlab.ui.control.Label
        SkirtDiamEditField              matlab.ui.control.NumericEditField
        BodyLengthEditFieldLabel        matlab.ui.control.Label
        BodyLengthEditField             matlab.ui.control.NumericEditField
        SkirtLengthEditFieldLabel       matlab.ui.control.Label
        SkirtLengthEditField            matlab.ui.control.NumericEditField
        LengthtoFinsEditFieldLabel      matlab.ui.control.Label
        LengthtoFinsEditField           matlab.ui.control.NumericEditField
        FinTipLengthEditFieldLabel      matlab.ui.control.Label
        FinTipLengthEditField           matlab.ui.control.NumericEditField
        FinLengthEditFieldLabel         matlab.ui.control.Label
        FinLengthEditField              matlab.ui.control.NumericEditField
        FinHeightEditFieldLabel         matlab.ui.control.Label
        FinHeightEditField              matlab.ui.control.NumericEditField
        FintoTipEditFieldLabel          matlab.ui.control.Label
        FintoTipEditField               matlab.ui.control.NumericEditField
        FinThicknessEditFieldLabel      matlab.ui.control.Label
        FinThicknessEditField           matlab.ui.control.NumericEditField
        UpdateLengthParametersButton    matlab.ui.control.Button
        RecoveryParametersPanel         matlab.ui.container.Panel
        DrogueRopeLengthEditFieldLabel  matlab.ui.control.Label
        DrogueRopeLengthEditField       matlab.ui.control.NumericEditField
        DrogueDragCoefficientEditFieldLabel  matlab.ui.control.Label
        DrogueDragCoefficientEditField  matlab.ui.control.NumericEditField
        DrogueDiamEditFieldLabel        matlab.ui.control.Label
        DrogueDiamEditField             matlab.ui.control.NumericEditField
        DrogueMassEditFieldLabel        matlab.ui.control.Label
        DrogueMassEditField             matlab.ui.control.NumericEditField
        DrogueRopeElasticityEditFieldLabel  matlab.ui.control.Label
        DrogueRopeElasticityEditField    matlab.ui.control.NumericEditField
        UpdateRecoveryParametersButton  matlab.ui.control.Button
        DrogueRopeDampingEditFieldLabel  matlab.ui.control.Label
        DrogueRopeDampingEditField       matlab.ui.control.NumericEditField
        MainRopeLengthEditFieldLabel    matlab.ui.control.Label
        MainRopeLengthEditField         matlab.ui.control.NumericEditField
        MainDragCoefficientEditFieldLabel  matlab.ui.control.Label
        MainDragCoefficientEditField    matlab.ui.control.NumericEditField
        MainDiamEditFieldLabel          matlab.ui.control.Label
        MainDiamEditField               matlab.ui.control.NumericEditField
        MainMassEditFieldLabel          matlab.ui.control.Label
        MainMassEditField               matlab.ui.control.NumericEditField
        MainRopeElasticityEditFieldLabel  matlab.ui.control.Label
        MainRopeElasticityEditField     matlab.ui.control.NumericEditField
        MainRopeDampingEditFieldLabel   matlab.ui.control.Label
        MainRopeDampingEditField        matlab.ui.control.NumericEditField
        Panel                           matlab.ui.container.Panel
        SaveRocketButton                matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
        SimulateandReportPanel          matlab.ui.container.Panel
        GenerateAnimationCheckBox       matlab.ui.control.CheckBox
        GifNameEditFieldLabel           matlab.ui.control.Label
        GifNameEditField                matlab.ui.control.EditField
        SaveGifofAnimationCheckBox      matlab.ui.control.CheckBox
        GenerateMassReportCheckBox      matlab.ui.control.CheckBox
        GenerateRecoveryReportCheckBox  matlab.ui.control.CheckBox
        RunSimulationButton             matlab.ui.control.Button
        LaunchRailHeightEditFieldLabel  matlab.ui.control.Label
        LaunchRailHeightEditField       matlab.ui.control.NumericEditField
        CompressibilityModelDropDownLabel  matlab.ui.control.Label
        CompressibilityModelDropDown    matlab.ui.control.DropDown
        LaunchLatEditFieldLabel         matlab.ui.control.Label
        LaunchLatEditField              matlab.ui.control.NumericEditField
        LaunchLongEditFieldLabel        matlab.ui.control.Label
        LaunchLongEditField             matlab.ui.control.NumericEditField
        LaunchAltitudeEditFieldLabel    matlab.ui.control.Label
        LaunchAltitudeEditField         matlab.ui.control.NumericEditField
        SavedSimulationsPanel           matlab.ui.container.Panel
        SavedSimulationsListBoxLabel    matlab.ui.control.Label
        SavedSimulationsListBox         matlab.ui.control.ListBox
        LoadSimulationDataButton        matlab.ui.control.Button
        
        R % Rocket Object
        Reco % Recovery Object
    end
    % Plotting
    methods (Access = private)
        
        % Plot Rocket
        function plotRocket(app,name)
            cla(app.UIAxes);
            hold(app.UIAxes,'on');
            plot(app.UIAxes,app.R.cone.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(app.UIAxes,app.R.body.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(app.UIAxes,app.R.skirt.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(app.UIAxes,app.R.LOX.polyshape,'FaceColor','g');
            plot(app.UIAxes,app.R.CH4.polyshape,'FaceColor','g');
            plot(app.UIAxes,app.R.finR.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(app.UIAxes,app.R.finL.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(app.UIAxes,app.R.finF.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(app.UIAxes,app.R.finB.polyshape,'FaceColor',[0 0.4470 0.7410]);
            plot(app.UIAxes,app.R.rocketCentroid{1}(1),app.R.rocketCentroid{1}(3),'r*');
            view(app.UIAxes,270,-90);
            axis(app.UIAxes,'equal');
            if nargin == 2
                title(app.UIAxes, name);
            else
                title(app.UIAxes, 'Unnamed Rocket');
            end
            app.UIAxes.PlotBoxAspectRatio = [1 6.32121212121212 1];
        end
    end
    
    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadRocketButton
        function LoadRocketButtonPushed(app, event)
            Rstruc = load([cd '\SavedRockets\' app.SavedRocketsListBox.Value]);
            app.R = Rstruc.R;
            app.Reco = Rstruc.Reco;
            plotRocket(app,app.SavedRocketsListBox.Value);
            
            app.NoseLengthEditField.Value = app.R.ln;
            app.NoseDiamEditField.Value = app.R.dn;
            app.BodyDiamEditField.Value = app.R.db;
            app.DiamatFinsEditField.Value = app.R.df;
            app.EndBodyDiamEditField.Value = app.R.du;
            app.SkirtDiamEditField.Value = app.R.dd;
            app.BodyLengthEditField.Value = app.R.lb;
            app.SkirtLengthEditField.Value = app.R.lc;
            app.LengthtoFinsEditField.Value = app.R.Xf;
            app.FinTipLengthEditField.Value = app.R.lt;
            app.FinLengthEditField.Value = app.R.lr;
            app.FinHeightEditField.Value = app.R.ls;
            app.FintoTipEditField.Value = app.R.lw;
            app.FinThicknessEditField.Value = app.R.Tf;
            
            app.LOXCOMEditField.Value = app.R.lLOXCM;
            app.LOXTankLengthEditField.Value = app.R.lLOX;
            app.LOXTankDiamEditField.Value = app.R.dLOX;
            app.LOXMassEditField.Value = app.R.m0LOX;
            
            app.CH4COMEditField.Value = app.R.lCH4CM;
            app.CH4TankLengthEditField.Value = app.R.lCH4;
            app.CH4TankDiamEditField.Value = app.R.dCH4;
            app.CH4MassEditField.Value = app.R.m0CH4;
            
            app.BurntimeEditField.Value = app.R.burntime;
            app.DrymassEditField.Value = app.R.drymass;
            
            app.DrogueRopeLengthEditField.Value = app.Reco.Drogue.L;
            app.DrogueRopeElasticityEditField.Value = app.Reco.Drogue.k;
            app.DrogueRopeDampingEditField.Value = app.Reco.Drogue.c;
            app.DrogueDragCoefficientEditField.Value = app.Reco.Drogue.Cd;
            app.DrogueDiamEditField.Value = app.Reco.Drogue.A;
            app.DrogueMassEditField.Value = app.Reco.Drogue.m;
            
            app.MainRopeLengthEditField.Value = app.Reco.Main.L;
            app.MainRopeElasticityEditField.Value = app.Reco.Main.k;
            app.MainRopeDampingEditField.Value = app.Reco.Main.c;
            app.MainDragCoefficientEditField.Value = app.Reco.Main.Cd;
            app.MainDiamEditField.Value = app.Reco.Main.A;
            app.MainMassEditField.Value = app.Reco.Main.m;
        end

        % Button pushed function: UpdateLengthParametersButton
        function UpdateLengthParametersButtonPushed(app, event)
            
            app.R = rocket(app.NoseLengthEditField.Value,app.NoseDiamEditField.Value,app.BodyDiamEditField.Value,app.DiamatFinsEditField.Value,app.EndBodyDiamEditField.Value,app.SkirtDiamEditField.Value,app.BodyLengthEditField.Value,app.SkirtLengthEditField.Value,app.LengthtoFinsEditField.Value,app.FinLengthEditField.Value,app.FinTipLengthEditField.Value,app.FinHeightEditField.Value,app.FintoTipEditField.Value,app.FinThicknessEditField.Value,app.R.lLOXCM,app.R.lLOX,app.R.dLOX,app.R.m0LOX,app.R.lCH4CM,app.R.lCH4,app.R.dCH4,app.R.m0CH4,app.R.burntime,app.R.drymass,app.R.dt);
            
            plotRocket(app);
        end

        % Button pushed function: UpdateMassParametersButton
        function UpdateMassParametersButtonPushed(app, event)
            
            if app.GenerateCustomThrustCurveCheckBox.Value
                app.BurntimeEditField.Value = -1;
            end
            app.R = rocket(app.R.ln,app.R.dn,app.R.db,app.R.df,app.R.du,app.R.dd,app.R.lb,app.R.lc,app.R.Xf,app.R.lr,app.R.lt,app.R.ls,app.R.lw,app.R.Tf,app.LOXCOMEditField.Value,app.LOXTankLengthEditField.Value,app.LOXTankDiamEditField.Value,app.LOXMassEditField.Value,app.CH4COMEditField.Value,app.CH4TankLengthEditField.Value,app.CH4TankDiamEditField.Value,app.CH4MassEditField.Value,app.BurntimeEditField.Value,app.DrymassEditField.Value,app.R.dt);
            
            plotRocket(app);
        end

        % Button pushed function: UpdateRecoveryParametersButton
        function UpdateRecoveryParametersButtonPushed(app, event)
            
            app.Reco = recovery(app.DrogueRopeLengthEditField.Value,app.DrogueRopeElasticityEditField.Value,app.DrogueRopeDampingEditField.Value,app.DrogueDragCoefficientEditField.Value,app.DrogueDiamEditField.Value,app.DrogueMassEditField.Value,app.MainRopeLengthEditField.Value,app.MainRopeElasticityEditField.Value,app.MainRopeDampingEditField.Value,app.MainDragCoefficientEditField.Value,app.MainDiamEditField.Value,app.MainMassEditField.Value);
            
        end

        % Button pushed function: SaveRocketButton
        function SaveRocketButtonPushed(app, event)
            saveName = inputdlg('Enter Save Name:');
            R = app.R; Reco = app.Reco;
            save([cd '\SavedRockets\' saveName{1} '.mat'],'R','Reco');
            
            rocketList = dir([cd '\SavedRockets']);
            app.SavedRocketsListBox.Items = {rocketList(3:end).name};
            
            plotRocket(app,[saveName{1} '.mat']);
        end

        % Button pushed function: RunSimulationButton
        function RunSimulationButtonPushed(app, event)
            
            msg = msgbox('Simulation Running...');
            [t,app.R,app.Reco,rocketState,drogueState,mainState,a,Ma,FR,nDrogueDeploy,Ldrogue,Tdrogue,Vldrogue,nMainDeploy,Lmain,Tmain,Vlmain,RopeStart] = rocketTrajectory(app.R,app.Reco,app.LaunchRailHeightEditField.Value,app.LaunchLatEditField.Value,app.LaunchLongEditField.Value,app.LaunchAltitudeEditField.Value,app.CompressibilityModelDropDown.Value,1/1000);
            close(msg);

            msgbox('Simulation Complete');
            
            generatePlots(t,rocketState,drogueState,mainState,a,Ma,app.R,nDrogueDeploy,app.Reco,Ldrogue,Tdrogue,Vldrogue,nMainDeploy,Lmain,Tmain,Vlmain,app.GenerateMassReportCheckBox.Value,app.GenerateRecoveryReportCheckBox.Value);
            
            if app.GenerateAnimationCheckBox.Value
                rocketParachuteAnimation(t,app.R,app.Reco,rocketState,drogueState,mainState,a,FR,nDrogueDeploy,nMainDeploy,RopeStart,H,app.SaveGifofAnimationCheckBox.Value,app.GifNameEditField.Value);
            end
            
            
            saveData = questdlg('Save Simulation Data?','Save','Yes','No','No');
            if isequal(saveData,'Yes')
                saveName = inputdlg('Enter Save Name:');
                R = app.R;
                Reco = app.Reco;
                LaunchLat = app.LaunchLatEditField.Value;
                LaunchLong = app.LaunchLongEditField.Value;
                LaunchAlt = app.LaunchAltitudeEditField.Value;
                H = app.LaunchRailHeightEditField.Value;
                compressibility = app.CompressibilityModelDropDown.Value;
                save([cd '\SavedSims\' saveName{1} '.mat'],'R','Reco','t','X','Y','Z','V','A','F','Fd','Fn','a','alpha','beta','gamma','Re','Ma','Stability','Xdrogue','Ydrogue','Zdrogue','Vdrogue','Adrogue','nDrogueDeploy','RopeStart','Ldrogue','Tdrogue','Vldrogue','Xmain','Ymain','Zmain','Vmain','Amain','nMainDeploy','Lmain','Tmain','Vlmain','LaunchLat','LaunchLong','LaunchAlt','H','compressibility');
                simsList = dir([cd '\SavedSims']);
                app.SavedSimulationsListBox.Items = {simsList(3:end).name};
            end
        end
        % Button pushed function: LoadSimulationDataButton
        function LoadSimulationDataButtonPushed(app, event)
            S = load([cd '\SavedSims\' app.SavedSimulationsListBox.Value]);
            app.R = S.R;
            app.Reco = S.Reco;
            plotRocket(app,app.SavedSimulationsListBox.Value);
            
            app.NoseLengthEditField.Value = app.R.ln;
            app.NoseDiamEditField.Value = app.R.dn;
            app.BodyDiamEditField.Value = app.R.db;
            app.DiamatFinsEditField.Value = app.R.df;
            app.EndBodyDiamEditField.Value = app.R.du;
            app.SkirtDiamEditField.Value = app.R.dd;
            app.BodyLengthEditField.Value = app.R.lb;
            app.SkirtLengthEditField.Value = app.R.lc;
            app.LengthtoFinsEditField.Value = app.R.Xf;
            app.FinTipLengthEditField.Value = app.R.lt;
            app.FinLengthEditField.Value = app.R.lr;
            app.FinHeightEditField.Value = app.R.ls;
            app.FintoTipEditField.Value = app.R.lw;
            app.FinThicknessEditField.Value = app.R.Tf;
            
            app.LOXCOMEditField.Value = app.R.lLOXCM;
            app.LOXTankLengthEditField.Value = app.R.lLOX;
            app.LOXTankDiamEditField.Value = app.R.dLOX;
            app.LOXMassEditField.Value = app.R.m0LOX;
            
            app.CH4COMEditField.Value = app.R.lCH4CM;
            app.CH4TankLengthEditField.Value = app.R.lCH4;
            app.CH4TankDiamEditField.Value = app.R.dCH4;
            app.CH4MassEditField.Value = app.R.m0CH4;
            
            app.BurntimeEditField.Value = app.R.burntime;
            app.DrymassEditField.Value = app.R.drymass;
            
            app.DrogueRopeLengthEditField.Value = app.Reco.Drogue.L;
            app.DrogueRopeElasticityEditField.Value = app.Reco.Drogue.k;
            app.DrogueRopeDampingEditField.Value = app.Reco.Drogue.c;
            app.DrogueDragCoefficientEditField.Value = app.Reco.Drogue.Cd;
            app.DrogueDiamEditField.Value = app.Reco.Drogue.A;
            app.DrogueMassEditField.Value = app.Reco.Drogue.m;
            
            app.MainRopeLengthEditField.Value = app.Reco.Main.L;
            app.MainRopeElasticityEditField.Value = app.Reco.Main.k;
            app.MainRopeDampingEditField.Value = app.Reco.Main.c;
            app.MainDragCoefficientEditField.Value = app.Reco.Main.Cd;
            app.MainDiamEditField.Value = app.Reco.Main.A;
            app.MainMassEditField.Value = app.Reco.Main.m;
            
            app.LaunchLatEditField.Value = S.LaunchLat;
            app.LaunchLongEditField.Value = S.LaunchLong;
            app.LaunchAltitudeEditField.Value = S.LaunchAlt;
            app.LaunchRailHeightEditField.Value = S.H;
            app.CompressibilityModelDropDown.Value = S.compressibility;
            
            generatePlots(S.t,S.X,S.Y,S.Z,S.V,S.A,S.a,S.alpha,S.beta,S.gamma,S.Ma,S.Stability,app.R,S.Zdrogue,S.Vdrogue,S.Adrogue,S.nDrogueDeploy,app.Reco,S.Ldrogue,S.Tdrogue,S.Vldrogue,S.Zmain,S.Vmain,S.Amain,S.nMainDeploy,S.Lmain,S.Tmain,S.Vlmain,app.GenerateMassReportCheckBox.Value,app.GenerateRecoveryReportCheckBox.Value);
            
             if app.GenerateAnimationCheckBox.Value
                rocketParachuteAnimation(S.alpha,S.beta,S.gamma,app.R,S.V,S.F,S.a,S.X,S.Y,S.Z,S.t,S.RopeStart,app.Reco,S.Xdrogue,S.Ydrogue,S.Zdrogue,S.Vdrogue,S.nDrogueDeploy,S.Xmain,S.Ymain,S.Zmain,S.Vmain,S.nMainDeploy,S.H,app.SaveGifofAnimationCheckBox.Value,app.GifNameEditField.Value)
            end
        end
    end
    
    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1712 458];
            app.UIFigure.Name = 'MATLAB App';

            % Create MassParametersPanel
            app.MassParametersPanel = uipanel(app.UIFigure);
            app.MassParametersPanel.TitlePosition = 'centertop';
            app.MassParametersPanel.Title = 'Mass Parameters';
            app.MassParametersPanel.Position = [313 1 523 193];

            % Create LOXCOMEditFieldLabel
            app.LOXCOMEditFieldLabel = uilabel(app.MassParametersPanel);
            app.LOXCOMEditFieldLabel.HorizontalAlignment = 'right';
            app.LOXCOMEditFieldLabel.Position = [35 141 61 22];
            app.LOXCOMEditFieldLabel.Text = 'LOX COM';

            % Create LOXCOMEditField
            app.LOXCOMEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.LOXCOMEditField.Position = [111 141 44 22];

            % Create LOXTankLengthEditFieldLabel
            app.LOXTankLengthEditFieldLabel = uilabel(app.MassParametersPanel);
            app.LOXTankLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.LOXTankLengthEditFieldLabel.Position = [-2 111 98 22];
            app.LOXTankLengthEditFieldLabel.Text = 'LOX Tank Length';

            % Create LOXTankLengthEditField
            app.LOXTankLengthEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.LOXTankLengthEditField.Position = [111 111 44 22];

            % Create LOXTankDiamEditFieldLabel
            app.LOXTankDiamEditFieldLabel = uilabel(app.MassParametersPanel);
            app.LOXTankDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.LOXTankDiamEditFieldLabel.Position = [6 79 90 22];
            app.LOXTankDiamEditFieldLabel.Text = 'LOX Tank Diam';

            % Create LOXTankDiamEditField
            app.LOXTankDiamEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.LOXTankDiamEditField.Position = [111 79 44 22];

            % Create LOXMassEditFieldLabel
            app.LOXMassEditFieldLabel = uilabel(app.MassParametersPanel);
            app.LOXMassEditFieldLabel.HorizontalAlignment = 'right';
            app.LOXMassEditFieldLabel.Position = [34 45 62 22];
            app.LOXMassEditFieldLabel.Text = 'LOX Mass';

            % Create LOXMassEditField
            app.LOXMassEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.LOXMassEditField.Position = [111 45 44 22];

            % Create CH4COMEditFieldLabel
            app.CH4COMEditFieldLabel = uilabel(app.MassParametersPanel);
            app.CH4COMEditFieldLabel.HorizontalAlignment = 'right';
            app.CH4COMEditFieldLabel.Position = [195 141 61 22];
            app.CH4COMEditFieldLabel.Text = 'CH4 COM';

            % Create CH4COMEditField
            app.CH4COMEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.CH4COMEditField.Position = [271 141 44 22];

            % Create CH4TankLengthEditFieldLabel
            app.CH4TankLengthEditFieldLabel = uilabel(app.MassParametersPanel);
            app.CH4TankLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.CH4TankLengthEditFieldLabel.Position = [158 111 98 22];
            app.CH4TankLengthEditFieldLabel.Text = 'CH4 Tank Length';

            % Create CH4TankLengthEditField
            app.CH4TankLengthEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.CH4TankLengthEditField.Position = [271 111 44 22];

            % Create CH4TankDiamEditFieldLabel
            app.CH4TankDiamEditFieldLabel = uilabel(app.MassParametersPanel);
            app.CH4TankDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.CH4TankDiamEditFieldLabel.Position = [166 79 90 22];
            app.CH4TankDiamEditFieldLabel.Text = 'CH4 Tank Diam';

            % Create CH4TankDiamEditField
            app.CH4TankDiamEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.CH4TankDiamEditField.Position = [271 79 44 22];

            % Create CH4MassEditFieldLabel
            app.CH4MassEditFieldLabel = uilabel(app.MassParametersPanel);
            app.CH4MassEditFieldLabel.HorizontalAlignment = 'right';
            app.CH4MassEditFieldLabel.Position = [194 45 62 22];
            app.CH4MassEditFieldLabel.Text = 'CH4 Mass';

            % Create CH4MassEditField
            app.CH4MassEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.CH4MassEditField.Position = [271 45 44 22];

            % Create DrymassEditFieldLabel
            app.DrymassEditFieldLabel = uilabel(app.MassParametersPanel);
            app.DrymassEditFieldLabel.HorizontalAlignment = 'right';
            app.DrymassEditFieldLabel.Position = [342 141 53 22];
            app.DrymassEditFieldLabel.Text = 'Drymass';

            % Create DrymassEditField
            app.DrymassEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.DrymassEditField.Position = [410 141 44 22];

            % Create GenerateCustomThrustCurveCheckBox
            app.GenerateCustomThrustCurveCheckBox = uicheckbox(app.MassParametersPanel);
            app.GenerateCustomThrustCurveCheckBox.Text = 'Generate Custom Thrust Curve';
            app.GenerateCustomThrustCurveCheckBox.Position = [332 111 190 22];

            % Create BurntimeEditFieldLabel
            app.BurntimeEditFieldLabel = uilabel(app.MassParametersPanel);
            app.BurntimeEditFieldLabel.HorizontalAlignment = 'right';
            app.BurntimeEditFieldLabel.Position = [341 79 54 22];
            app.BurntimeEditFieldLabel.Text = 'Burntime';

            % Create BurntimeEditField
            app.BurntimeEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.BurntimeEditField.Position = [410 79 44 22];

            % Create ThrustEditFieldLabel
            app.ThrustEditFieldLabel = uilabel(app.MassParametersPanel);
            app.ThrustEditFieldLabel.HorizontalAlignment = 'right';
            app.ThrustEditFieldLabel.Position = [355 47 40 22];
            app.ThrustEditFieldLabel.Text = 'Thrust';

            % Create ThrustEditField
            app.ThrustEditField = uieditfield(app.MassParametersPanel, 'numeric');
            app.ThrustEditField.Position = [410 47 44 22];

            % Create UpdateMassParametersButton
            app.UpdateMassParametersButton = uibutton(app.MassParametersPanel, 'push');
            app.UpdateMassParametersButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateMassParametersButtonPushed, true);
            app.UpdateMassParametersButton.Position = [186 15 152 22];
            app.UpdateMassParametersButton.Text = 'Update Mass Parameters';

            % Create LoadRocketPanel
            app.LoadRocketPanel = uipanel(app.UIFigure);
            app.LoadRocketPanel.TitlePosition = 'centertop';
            app.LoadRocketPanel.Title = 'Load Rocket';
            app.LoadRocketPanel.Position = [1 293 314 166];

            % Create SavedRocketsListBoxLabel
            app.SavedRocketsListBoxLabel = uilabel(app.LoadRocketPanel);
            app.SavedRocketsListBoxLabel.HorizontalAlignment = 'right';
            app.SavedRocketsListBoxLabel.Position = [17 54 86 61];
            app.SavedRocketsListBoxLabel.Text = 'Saved Rockets';

            % Create SavedRocketsListBox
            app.SavedRocketsListBox = uilistbox(app.LoadRocketPanel);
            rocketList = dir([cd '\SavedRockets']);
            app.SavedRocketsListBox.Items = {rocketList(3:end).name};
            app.SavedRocketsListBox.Position = [118 43 131 74];

            % Create LoadRocketButton
            app.LoadRocketButton = uibutton(app.LoadRocketPanel, 'push');
            app.LoadRocketButton.ButtonPushedFcn = createCallbackFcn(app, @LoadRocketButtonPushed, true);
            app.LoadRocketButton.Position = [134 11 100 22];
            app.LoadRocketButton.Text = 'Load Rocket';

            % Create RocketLengthsPanel
            app.RocketLengthsPanel = uipanel(app.UIFigure);
            app.RocketLengthsPanel.TitlePosition = 'centertop';
            app.RocketLengthsPanel.Title = 'Rocket Lengths';
            app.RocketLengthsPanel.Position = [0 4 314 290];

            % Create NoseLengthEditFieldLabel
            app.NoseLengthEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.NoseLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.NoseLengthEditFieldLabel.Position = [24 237 74 22];
            app.NoseLengthEditFieldLabel.Text = 'Nose Length';

            % Create NoseLengthEditField
            app.NoseLengthEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.NoseLengthEditField.Position = [113 237 47 22];

            % Create NoseDiamEditFieldLabel
            app.NoseDiamEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.NoseDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.NoseDiamEditFieldLabel.Position = [176 237 65 22];
            app.NoseDiamEditFieldLabel.Text = 'Nose Diam';

            % Create NoseDiamEditField
            app.NoseDiamEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.NoseDiamEditField.Position = [256 237 44 22];

            % Create BodyDiamEditFieldLabel
            app.BodyDiamEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.BodyDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.BodyDiamEditFieldLabel.Position = [34 206 64 22];
            app.BodyDiamEditFieldLabel.Text = 'Body Diam';

            % Create BodyDiamEditField
            app.BodyDiamEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.BodyDiamEditField.Position = [113 206 47 22];

            % Create DiamatFinsEditFieldLabel
            app.DiamatFinsEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.DiamatFinsEditFieldLabel.HorizontalAlignment = 'right';
            app.DiamatFinsEditFieldLabel.Position = [168 206 73 22];
            app.DiamatFinsEditFieldLabel.Text = 'Diam at Fins';

            % Create DiamatFinsEditField
            app.DiamatFinsEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.DiamatFinsEditField.Position = [256 206 44 22];

            % Create EndBodyDiamEditFieldLabel
            app.EndBodyDiamEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.EndBodyDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.EndBodyDiamEditFieldLabel.Position = [9 176 89 22];
            app.EndBodyDiamEditFieldLabel.Text = 'End Body Diam';

            % Create EndBodyDiamEditField
            app.EndBodyDiamEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.EndBodyDiamEditField.Position = [113 176 47 22];

            % Create SkirtDiamEditFieldLabel
            app.SkirtDiamEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.SkirtDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.SkirtDiamEditFieldLabel.Position = [180 176 61 22];
            app.SkirtDiamEditFieldLabel.Text = 'Skirt Diam';

            % Create SkirtDiamEditField
            app.SkirtDiamEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.SkirtDiamEditField.Position = [256 176 44 22];

            % Create BodyLengthEditFieldLabel
            app.BodyLengthEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.BodyLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.BodyLengthEditFieldLabel.Position = [25 146 73 22];
            app.BodyLengthEditFieldLabel.Text = 'Body Length';

            % Create BodyLengthEditField
            app.BodyLengthEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.BodyLengthEditField.Position = [113 146 47 22];

            % Create SkirtLengthEditFieldLabel
            app.SkirtLengthEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.SkirtLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.SkirtLengthEditFieldLabel.Position = [171 146 70 22];
            app.SkirtLengthEditFieldLabel.Text = 'Skirt Length';

            % Create SkirtLengthEditField
            app.SkirtLengthEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.SkirtLengthEditField.Position = [256 146 44 22];

            % Create LengthtoFinsEditFieldLabel
            app.LengthtoFinsEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.LengthtoFinsEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthtoFinsEditFieldLabel.Position = [16 116 82 22];
            app.LengthtoFinsEditFieldLabel.Text = 'Length to Fins';

            % Create LengthtoFinsEditField
            app.LengthtoFinsEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.LengthtoFinsEditField.Position = [113 116 47 22];

            % Create FinTipLengthEditFieldLabel
            app.FinTipLengthEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.FinTipLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.FinTipLengthEditFieldLabel.Position = [159 116 82 22];
            app.FinTipLengthEditFieldLabel.Text = 'Fin Tip Length';

            % Create FinTipLengthEditField
            app.FinTipLengthEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.FinTipLengthEditField.Position = [256 116 44 22];

            % Create FinLengthEditFieldLabel
            app.FinLengthEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.FinLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.FinLengthEditFieldLabel.Position = [36 82 62 22];
            app.FinLengthEditFieldLabel.Text = 'Fin Length';

            % Create FinLengthEditField
            app.FinLengthEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.FinLengthEditField.Position = [113 82 47 22];

            % Create FinHeightEditFieldLabel
            app.FinHeightEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.FinHeightEditFieldLabel.HorizontalAlignment = 'right';
            app.FinHeightEditFieldLabel.Position = [181 82 60 22];
            app.FinHeightEditFieldLabel.Text = 'Fin Height';

            % Create FinHeightEditField
            app.FinHeightEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.FinHeightEditField.Position = [256 82 44 22];

            % Create FintoTipEditFieldLabel
            app.FintoTipEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.FintoTipEditFieldLabel.HorizontalAlignment = 'right';
            app.FintoTipEditFieldLabel.Position = [43 50 55 22];
            app.FintoTipEditFieldLabel.Text = 'Fin to Tip';

            % Create FintoTipEditField
            app.FintoTipEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.FintoTipEditField.Position = [113 50 47 22];

            % Create FinThicknessEditFieldLabel
            app.FinThicknessEditFieldLabel = uilabel(app.RocketLengthsPanel);
            app.FinThicknessEditFieldLabel.HorizontalAlignment = 'right';
            app.FinThicknessEditFieldLabel.Position = [161 50 80 22];
            app.FinThicknessEditFieldLabel.Text = 'Fin Thickness';

            % Create FinThicknessEditField
            app.FinThicknessEditField = uieditfield(app.RocketLengthsPanel, 'numeric');
            app.FinThicknessEditField.Position = [256 50 44 22];

            % Create UpdateLengthParametersButton
            app.UpdateLengthParametersButton = uibutton(app.RocketLengthsPanel, 'push');
            app.UpdateLengthParametersButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateLengthParametersButtonPushed, true);
            app.UpdateLengthParametersButton.Position = [77 20 160 22];
            app.UpdateLengthParametersButton.Text = 'Update Length Parameters';

            % Create RecoveryParametersPanel
            app.RecoveryParametersPanel = uipanel(app.UIFigure);
            app.RecoveryParametersPanel.TitlePosition = 'centertop';
            app.RecoveryParametersPanel.Title = 'Recovery Parameters';
            app.RecoveryParametersPanel.Position = [835 1 557 193];

            % Create DrogueRopeLengthEditFieldLabel
            app.DrogueRopeLengthEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.DrogueRopeLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.DrogueRopeLengthEditFieldLabel.Position = [20 141 117 22];
            app.DrogueRopeLengthEditFieldLabel.Text = 'Drogue Rope Length';

            % Create DrogueRopeLengthEditField
            app.DrogueRopeLengthEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.DrogueRopeLengthEditField.Position = [152 141 44 22];

            % Create DrogueDragCoefficientEditFieldLabel
            app.DrogueDragCoefficientEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.DrogueDragCoefficientEditFieldLabel.HorizontalAlignment = 'right';
            app.DrogueDragCoefficientEditFieldLabel.Position = [3 47 134 22];
            app.DrogueDragCoefficientEditFieldLabel.Text = 'Drogue Drag Coefficient';

            % Create DrogueDragCoefficientEditField
            app.DrogueDragCoefficientEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.DrogueDragCoefficientEditField.Position = [152 47 44 22];

            % Create DrogueDiamEditFieldLabel
            app.DrogueDiamEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.DrogueDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.DrogueDiamEditFieldLabel.Position = [412 142 76 22];
            app.DrogueDiamEditFieldLabel.Text = 'Drogue Diam';

            % Create DrogueDiamEditField
            app.DrogueDiamEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.DrogueDiamEditField.Position = [503 142 44 22];

            % Create DrogueMassEditFieldLabel
            app.DrogueMassEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.DrogueMassEditFieldLabel.HorizontalAlignment = 'right';
            app.DrogueMassEditFieldLabel.Position = [412 111 77 22];
            app.DrogueMassEditFieldLabel.Text = 'Drogue Mass';

            % Create DrogueMassEditField
            app.DrogueMassEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.DrogueMassEditField.Position = [503 111 44 22];

            % Create DrogueRopeElasticityEditFieldLabel
            app.DrogueRopeElasticityEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.DrogueRopeElasticityEditFieldLabel.HorizontalAlignment = 'right';
            app.DrogueRopeElasticityEditFieldLabel.Position = [16 112 121 22];
            app.DrogueRopeElasticityEditFieldLabel.Text = 'Drogue Rope Elasticity';

            % Create DrogueRopeElasticityEditField
            app.DrogueRopeElasticityEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.DrogueRopeElasticityEditField.Position = [152 112 44 22];

            % Create UpdateRecoveryParametersButton
            app.UpdateRecoveryParametersButton = uibutton(app.RecoveryParametersPanel, 'push');
            app.UpdateRecoveryParametersButton.ButtonPushedFcn = createCallbackFcn(app, @UpdateRecoveryParametersButtonPushed, true);
            app.UpdateRecoveryParametersButton.Position = [217 13 174 22];
            app.UpdateRecoveryParametersButton.Text = 'Update Recovery Parameters';

            % Create DrogueRopeDampingEditFieldLabel
            app.DrogueRopeDampingEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.DrogueRopeDampingEditFieldLabel.HorizontalAlignment = 'right';
            app.DrogueRopeDampingEditFieldLabel.Position = [15 79 122 22];
            app.DrogueRopeDampingEditFieldLabel.Text = 'Drogue Rope Damping';

            % Create DrogueRopeDampingEditField
            app.DrogueRopeDampingEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.DrogueRopeDampingEditField.Position = [152 79 44 22];

            % Create MainRopeLengthEditFieldLabel
            app.MainRopeLengthEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.MainRopeLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.MainRopeLengthEditFieldLabel.Position = [222 141 104 22];
            app.MainRopeLengthEditFieldLabel.Text = 'Main Rope Length';

            % Create MainRopeLengthEditField
            app.MainRopeLengthEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.MainRopeLengthEditField.Position = [341 141 44 22];

            % Create MainDragCoefficientEditFieldLabel
            app.MainDragCoefficientEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.MainDragCoefficientEditFieldLabel.HorizontalAlignment = 'right';
            app.MainDragCoefficientEditFieldLabel.Position = [205 47 121 22];
            app.MainDragCoefficientEditFieldLabel.Text = 'Main Drag Coefficient';

            % Create MainDragCoefficientEditField
            app.MainDragCoefficientEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.MainDragCoefficientEditField.Position = [341 47 44 22];

            % Create MainDiamEditFieldLabel
            app.MainDiamEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.MainDiamEditFieldLabel.HorizontalAlignment = 'right';
            app.MainDiamEditFieldLabel.Position = [425 79 63 22];
            app.MainDiamEditFieldLabel.Text = 'Main Diam';

            % Create MainDiamEditField
            app.MainDiamEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.MainDiamEditField.Position = [503 79 44 22];

            % Create MainMassEditFieldLabel
            app.MainMassEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.MainMassEditFieldLabel.HorizontalAlignment = 'right';
            app.MainMassEditFieldLabel.Position = [425 47 64 22];
            app.MainMassEditFieldLabel.Text = 'Main Mass';

            % Create MainMassEditField
            app.MainMassEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.MainMassEditField.Position = [503 47 44 22];

            % Create MainRopeElasticityEditFieldLabel
            app.MainRopeElasticityEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.MainRopeElasticityEditFieldLabel.HorizontalAlignment = 'right';
            app.MainRopeElasticityEditFieldLabel.Position = [212 111 114 22];
            app.MainRopeElasticityEditFieldLabel.Text = 'Main Rope Elasticity';

            % Create MainRopeElasticityEditField
            app.MainRopeElasticityEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.MainRopeElasticityEditField.Position = [341 111 44 22];

            % Create MainRopeDampingEditFieldLabel
            app.MainRopeDampingEditFieldLabel = uilabel(app.RecoveryParametersPanel);
            app.MainRopeDampingEditFieldLabel.HorizontalAlignment = 'right';
            app.MainRopeDampingEditFieldLabel.Position = [211 79 115 22];
            app.MainRopeDampingEditFieldLabel.Text = 'Main Rope Damping';

            % Create MainRopeDampingEditField
            app.MainRopeDampingEditField = uieditfield(app.RecoveryParametersPanel, 'numeric');
            app.MainRopeDampingEditField.Position = [341 79 44 22];

            % Create Panel
            app.Panel = uipanel(app.UIFigure);
            app.Panel.Position = [314 193 1078 266];

            % Create SaveRocketButton
            app.SaveRocketButton = uibutton(app.Panel, 'push');
            app.SaveRocketButton.ButtonPushedFcn = createCallbackFcn(app, @SaveRocketButtonPushed, true);
            app.SaveRocketButton.Position = [514 8 100 22];
            app.SaveRocketButton.Text = 'Save Rocket';

            % Create UIAxes
            app.UIAxes = uiaxes(app.Panel);
            title(app.UIAxes, 'Unnamed Rocket');
            app.UIAxes.PlotBoxAspectRatio = [6.32121212121212 1 1];
            app.UIAxes.Position = [11 36 1032 219];

            % Create SimulateandReportPanel
            app.SimulateandReportPanel = uipanel(app.UIFigure);
            app.SimulateandReportPanel.TitlePosition = 'centertop';
            app.SimulateandReportPanel.Title = 'Simulate and Report';
            app.SimulateandReportPanel.Position = [1391 1 316 293];

            % Create GenerateAnimationCheckBox
            app.GenerateAnimationCheckBox = uicheckbox(app.SimulateandReportPanel);
            app.GenerateAnimationCheckBox.Text = 'Generate Animation';
            app.GenerateAnimationCheckBox.Position = [20 244 129 22];

            % Create GifNameEditFieldLabel
            app.GifNameEditFieldLabel = uilabel(app.SimulateandReportPanel);
            app.GifNameEditFieldLabel.HorizontalAlignment = 'right';
            app.GifNameEditFieldLabel.Position = [75 214 56 22];
            app.GifNameEditFieldLabel.Text = 'Gif Name';

            % Create GifNameEditField
            app.GifNameEditField = uieditfield(app.SimulateandReportPanel, 'text');
            app.GifNameEditField.Position = [145 214 100 22];

            % Create SaveGifofAnimationCheckBox
            app.SaveGifofAnimationCheckBox = uicheckbox(app.SimulateandReportPanel);
            app.SaveGifofAnimationCheckBox.Text = 'Save Gif of Animation';
            app.SaveGifofAnimationCheckBox.Position = [165 244 138 22];

            % Create GenerateMassReportCheckBox
            app.GenerateMassReportCheckBox = uicheckbox(app.SimulateandReportPanel);
            app.GenerateMassReportCheckBox.Text = 'Generate Mass Report';
            app.GenerateMassReportCheckBox.Position = [6 179 143 22];

            % Create GenerateRecoveryReportCheckBox
            app.GenerateRecoveryReportCheckBox = uicheckbox(app.SimulateandReportPanel);
            app.GenerateRecoveryReportCheckBox.Text = 'Generate Recovery Report';
            app.GenerateRecoveryReportCheckBox.Position = [151 179 165 22];

            % Create RunSimulationButton
            app.RunSimulationButton = uibutton(app.SimulateandReportPanel, 'push');
            app.RunSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @RunSimulationButtonPushed, true);
            app.RunSimulationButton.Position = [110 26 100 22];
            app.RunSimulationButton.Text = 'Run Simulation';

            % Create LaunchRailHeightEditFieldLabel
            app.LaunchRailHeightEditFieldLabel = uilabel(app.SimulateandReportPanel);
            app.LaunchRailHeightEditFieldLabel.HorizontalAlignment = 'right';
            app.LaunchRailHeightEditFieldLabel.Position = [78 147 107 22];
            app.LaunchRailHeightEditFieldLabel.Text = 'Launch Rail Height';

            % Create LaunchRailHeightEditField
            app.LaunchRailHeightEditField = uieditfield(app.SimulateandReportPanel, 'numeric');
            app.LaunchRailHeightEditField.Position = [200 147 65 22];
            app.LaunchRailHeightEditField.Value = 15;

            % Create CompressibilityModelDropDownLabel
            app.CompressibilityModelDropDownLabel = uilabel(app.SimulateandReportPanel);
            app.CompressibilityModelDropDownLabel.HorizontalAlignment = 'right';
            app.CompressibilityModelDropDownLabel.Position = [39 119 123 22];
            app.CompressibilityModelDropDownLabel.Text = 'Compressibility Model';

            % Create CompressibilityModelDropDown
            app.CompressibilityModelDropDown = uidropdown(app.SimulateandReportPanel);
            app.CompressibilityModelDropDown.Items = {'Laitone-Ackeret', 'Prandtl-Glauert', 'Karman-Tsien'};
            app.CompressibilityModelDropDown.Position = [177 119 119 22];
            app.CompressibilityModelDropDown.Value = 'Laitone-Ackeret';

            % Create LaunchLatEditFieldLabel
            app.LaunchLatEditFieldLabel = uilabel(app.SimulateandReportPanel);
            app.LaunchLatEditFieldLabel.HorizontalAlignment = 'right';
            app.LaunchLatEditFieldLabel.Position = [11 90 65 22];
            app.LaunchLatEditFieldLabel.Text = 'Launch Lat';

            % Create LaunchLatEditField
            app.LaunchLatEditField = uieditfield(app.SimulateandReportPanel, 'numeric');
            app.LaunchLatEditField.Position = [91 90 55 22];
            app.LaunchLatEditField.Value = 35.35;

            % Create LaunchLongEditFieldLabel
            app.LaunchLongEditFieldLabel = uilabel(app.SimulateandReportPanel);
            app.LaunchLongEditFieldLabel.HorizontalAlignment = 'right';
            app.LaunchLongEditFieldLabel.Position = [151 90 75 22];
            app.LaunchLongEditFieldLabel.Text = 'Launch Long';

            % Create LaunchLongEditField
            app.LaunchLongEditField = uieditfield(app.SimulateandReportPanel, 'numeric');
            app.LaunchLongEditField.Position = [241 90 55 22];
            app.LaunchLongEditField.Value = -117.81;

            % Create LaunchAltitudeEditFieldLabel
            app.LaunchAltitudeEditFieldLabel = uilabel(app.SimulateandReportPanel);
            app.LaunchAltitudeEditFieldLabel.HorizontalAlignment = 'right';
            app.LaunchAltitudeEditFieldLabel.Position = [71 58 88 22];
            app.LaunchAltitudeEditFieldLabel.Text = 'Launch Altitude';

            % Create LaunchAltitudeEditField
            app.LaunchAltitudeEditField = uieditfield(app.SimulateandReportPanel, 'numeric');
            app.LaunchAltitudeEditField.Position = [174 58 68 22];
            app.LaunchAltitudeEditField.Value = 622;

            % Create SavedSimulationsPanel
            app.SavedSimulationsPanel = uipanel(app.UIFigure);
            app.SavedSimulationsPanel.Title = 'Saved Simulations';
            app.SavedSimulationsPanel.Position = [1391 293 316 166];

            % Create SavedSimulationsListBoxLabel
            app.SavedSimulationsListBoxLabel = uilabel(app.SavedSimulationsPanel);
            app.SavedSimulationsListBoxLabel.HorizontalAlignment = 'right';
            app.SavedSimulationsListBoxLabel.Position = [26 104 105 22];
            app.SavedSimulationsListBoxLabel.Text = 'Saved Simulations';

            % Create SavedSimulationsListBox
            app.SavedSimulationsListBox = uilistbox(app.SavedSimulationsPanel);
            simsList = dir([cd '\SavedSims']);
            app.SavedSimulationsListBox.Items = {simsList(3:end).name};
            app.SavedSimulationsListBox.Position = [146 54 100 74];

            % Create LoadSimulationDataButton
            app.LoadSimulationDataButton = uibutton(app.SavedSimulationsPanel, 'push');
            app.LoadSimulationDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadSimulationDataButtonPushed, true);
            app.LoadSimulationDataButton.Position = [71 11 130 22];
            app.LoadSimulationDataButton.Text = 'Load Simulation Data';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = RocketGUIV22

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end