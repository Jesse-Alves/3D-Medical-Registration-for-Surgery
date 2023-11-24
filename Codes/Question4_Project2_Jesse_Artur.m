%% ============ PROJECT2 - 3D Medical Registration ==============
%% =================== Jesse and Artur ========================
clear all; close all; clc

%% ==============>   Question 4

% Set the Initial Configuration
InitConfig();
InitPoseEff = Tbase_eff(1:3,4);
flag = 0;

%% Known Parameters
K = [400 0 380;
        0 400 285;
        0 0 1];

%% Obtaining the TCam1_Cam2  WITHOUT the NavSystem

%% First Position
trans = [20 0 0]';
trans = [5 0 -25]';
ang = 0;
MoveCamera(trans, ang);

TCam1_Inst = GetInstrumentPosition(flag);
pixel_target1 = GetTargetPosition(flag);

% DisplayConfig
% return

%% Second Position
trans = [2 0 -5]';
ang = 0;
MoveCamera(trans, ang);

%DisplayConfig

TCam2_Inst = GetInstrumentPosition(flag);
pixel_target2 = GetTargetPosition(flag);

disp('Transformation between the two postions of Camera ')
TCam1_Cam2 = TCam1_Inst*inv(TCam2_Inst)
[RCam1_Cam2, tCam1_Cam2] = GetRt(TCam1_Cam2);

%% Triangulation
m1 = inv(K)*[pixel_target1; 1];
m2 = inv(K)*[pixel_target2; 1];

M = [m1 -RCam1_Cam2*m2];

coef = pinv(M)*tCam1_Cam2;

tCam1_target = coef(1)*m1;
tCam2_target = coef(2)*m2;

%% The position of targert in the Base Frame
Tbase_eff = GetRobotCurrentPosition(flag)

tbase_target = Tbase_eff*(Teff_inst)*inv(TCam1_Inst)*[tCam1_target; 1];
tbase_target = tbase_target(1:3);

%% Moving the instrument to the Target
tbase_trocar = [-350;-100;800]; %Postion obtained in Question 3

%Finding the position of end-effector to do this task
direction = (tbase_target - tbase_trocar)/norm(tbase_target - tbase_trocar);
tbase_effx = tbase_target - 350*direction;

InitConfig
MoveEffPosition(tbase_effx)
DisplayConfig

disp('The Error of the Registration Process')
Error = ComputeTRE()



