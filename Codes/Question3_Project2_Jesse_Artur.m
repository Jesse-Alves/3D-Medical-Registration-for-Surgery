%% ============ PROJECT2 - 3D Medical Registration ==============
%% =================== Jesse and Artur ========================
clear all; close all; clc

%% ==============>   Question 3

% Set the Initial Configuration
InitConfig()
InitPoseEff = Tbase_eff(1:3,4);
TCam_MarkCam = Tcam_mark_cam;
flag = 0;

%% Known Parameters
K = [400 0 380;
        0 400 285;
        0 0 1];

%% Find Teff_markrobot
% Using the Effector-Marker Calibration

%% Getting 3 differents positions
% 1 Position
Tbase_eff1 = GetRobotCurrentPosition();
Tloc_MarkInst1 = GetLocalizerInformation(flag).mark(2).T;

% 2 Position
position2 = InitPoseEff + [10;20;15];
MoveEffPosition(position2)
%DisplayConfig()

Tbase_eff2 = GetRobotCurrentPosition();
Tloc_MarkInst2 = GetLocalizerInformation(flag).mark(2).T;

% 3 Position
position3 = InitPoseEff + [30;-50;5];
MoveEffPosition(position3)
%DisplayConfig()

Tbase_eff3 = GetRobotCurrentPosition();
Tloc_MarkInst3 = GetLocalizerInformation(flag).mark(2).T;

%% Calculating the transforms between Teffi_effj and TMarkInsti_Mrobotj and Slip R and t
%A
Teff2_eff1 = inv(Tbase_eff2)*Tbase_eff1; [Reff2_eff1, teff2_eff1] = GetRt(Teff2_eff1);
Teff3_eff1 = inv(Tbase_eff3)*Tbase_eff1; [Reff3_eff1, teff3_eff1] = GetRt(Teff3_eff1);
Teff3_eff2 = inv(Tbase_eff3)*Tbase_eff2; [Reff3_eff2, teff3_eff2] = GetRt(Teff3_eff2);

%B
TMarkInst2_MarkInst1 = inv(Tloc_MarkInst2)*Tloc_MarkInst1; [RMarkInst2_MarkInst1, tMarkInst2_MarkInst1] = GetRt(TMarkInst2_MarkInst1);
TMarkInst3_MarkInst1 = inv(Tloc_MarkInst3)*Tloc_MarkInst1; [RMarkInst3_MarkInst1, tMarkInst3_MarkInst1] = GetRt(TMarkInst3_MarkInst1);
TMarkInst3_MarkInst2 = inv(Tloc_MarkInst3)*Tloc_MarkInst2; [RMarkInst3_MarkInst2, tMarkInst3_MarkInst2] = GetRt(TMarkInst3_MarkInst2);

%% Obtaning the Eigvectors
[Vector,D] = eig(Reff2_eff1);   veff2_eff1 = Vector(1:3,1);   
[Vector,D] = eig(Reff3_eff1);   veff3_eff1 = Vector(1:3,1);
[Vector,D] = eig(Reff3_eff2);   veff3_eff2 = Vector(1:3,1);

[Vector,D] = eig(RMarkInst2_MarkInst1);   vM2M1 = Vector(1:3,1);
[Vector,D] = eig(RMarkInst3_MarkInst1);   vM3M1 = Vector(1:3,1);
[Vector,D] = eig(RMarkInst3_MarkInst2);   vM3M2 = Vector(1:3,1);

%% Computing Rx

V = [vM2M1' 0 0 0 0 0 0;
       0 0 0 vM2M1' 0 0 0;
       0 0 0 0 0 0 vM2M1';
       vM3M1' 0 0 0 0 0 0;
       0 0 0 vM3M1' 0 0 0;
       0 0 0 0 0 0 vM3M1';
       vM3M2' 0 0 0 0 0 0;
       0 0 0 vM3M2' 0 0 0;
       0 0 0 0 0 0 vM3M2'];

E = [veff2_eff1;veff3_eff1;veff3_eff2];

rx = inv(V)*E;

% The Rotation Matrix
Rx =reshape(rx,3,3)';

%% Computing tx
R_I = [Reff2_eff1-eye(3);
           Reff3_eff1-eye(3);
           Reff3_eff2-eye(3)];

R_t = [Rx*tMarkInst2_MarkInst1 - teff2_eff1;
           Rx*tMarkInst3_MarkInst1 - teff3_eff1;
           Rx*tMarkInst3_MarkInst2 - teff3_eff2];

tx = pinv(R_I)*R_t;
%% Obtaining the Follow Transformations:
clc
disp('The Transformation Effector-MarkRobot')
Teff_markrobot = [Rx tx;0 0 0 1]

disp('The Transformation Instrument-MarkRobot')
Tinst_markrobot = inv(Teff_inst)*Teff_markrobot

disp('The Transformation Base-LocSystem')
Tbase_loc = Tbase_eff3*(Teff_markrobot)*inv(Tloc_MarkInst3)

%% The Instrument Trocar Position
t_eff1 = Tbase_eff1(1:3,4);
t_eff2 = Tbase_eff2(1:3,4);

Z = [Tbase_eff1(1:3,3) -Tbase_eff2(1:3,3)];
f = pinv(Z)*(t_eff2 - t_eff1);
f_eff1 = f(1);
f_eff2 = f(2);

disp('The position of Instrument Trocar')
tbase_trocar = t_eff1 + f_eff1*Tbase_eff1(1:3,3)
    
%% Getting the Target Position

% Position 1
Tloc_MarkCam1 = GetLocalizerInformation(flag).mark(1).T;
pixel_target1 = GetTargetPosition(flag);

% Moving the Camera
trans = [20 0 0]';
ang = 0;
MoveCamera(trans, ang)

% Position 2
Tloc_MarkCam2 = GetLocalizerInformation(flag).mark(1).T;
pixel_target2 = GetTargetPosition(flag);

disp('Transformation between the two postions of Camera ')
Tloc_Cam1 = Tloc_MarkCam1*inv(TCam_MarkCam);
Tloc_Cam2 = Tloc_MarkCam2*inv(TCam_MarkCam);

TCam1_Cam2 = inv(Tloc_Cam1)*Tloc_Cam2
[RCam1_Cam2, tCam1_Cam2] = GetRt(TCam1_Cam2);

%% Triangulation
m1 = inv(K)*[pixel_target1; 1];
m2 = inv(K)*[pixel_target2; 1];

M = [m1 -RCam1_Cam2*m2];

coef = pinv(M)*tCam1_Cam2;

tCam1_target = coef(1)*m1;
tCam2_target = coef(2)*m2;

%% The position of targert in the Base Frame
tbase_target = Tbase_loc*(Tloc_Cam1)*[tCam1_target; 1];
tbase_target = tbase_target(1:3);

%% Moving the instrument to the Target

%Finding the position of end-effector to do this task
direction = (tbase_target - tbase_trocar)/norm(tbase_target - tbase_trocar);
tbase_effx = tbase_target - 350*direction;

InitConfig
MoveEffPosition(tbase_effx)
DisplayConfig

disp('The Error of the Registration Process')
Error = ComputeTRE()







