%% ============ PROJECT2 - 3D Medical Registration ==============
%% =================== Jesse and Artur ========================
clear all; close all; clc

%% ==============>   Question 6

% Set the Initial Configuration
InitConfig()
TCam_MarkCam = Tcam_mark_cam;
flag = 0;

%% Known Parameters
K = [400 0 380;
        0 400 285;
        0 0 1];

%% Steps of the Process
% 1 - Define a W
% 2 - Compute the Jacobian
% LOOP PROCESS
    % 3 - Compute the points with respect to Camera Frame
    %4 - Compute the Error w.r.t. the Camera Frame
    %5 - Find the Error in the Base Frame
    %6 - Integrate
    %7 - Apply the eff position in the robot 

%% Parameters
W = (1/5)*(eye(3));
J = [0 1 0;
       0 0 1;
       1 0 0];
Jinv = inv(J);

%% ========== 1 - Points in the camera frame

%% Obtain target in camera frame (triangulation)

% Position 1
Tnav_MarkCam1 = GetLocalizerInformation(flag).mark(1).T;
pixel_target1 = GetTargetPosition(flag);

% Moving the Camera to a fixed position
trans = [5 0 -25]';
ang = 0;
MoveCamera(trans, ang);

% Position 2
Tnav_MarkCam2 = GetLocalizerInformation(flag).mark(1).T;
pixel_target = GetTargetPosition(flag);

Tnav_Cam1 = Tnav_MarkCam1/TCam_MarkCam;
Tnav_Cam2 = Tnav_MarkCam2/TCam_MarkCam;

TCam1_Cam2 = Tnav_Cam1\Tnav_Cam2;
[RCam1_Cam2, tCam1_Cam2] = GetRt(TCam1_Cam2);

% Triangulation
m1 = inv(K)*[pixel_target1; 1];
m2 = inv(K)*[pixel_target; 1];
M = [m1 -RCam1_Cam2*m2];
coef = pinv(M)*tCam1_Cam2;
tCam1_target = coef(1)*m1;
tCam_target = coef(2)*m2;


%% Transformation Base-Camera

TCam_Inst = GetInstrumentPosition(flag);
tbase_trocar = [-350;-100;800];

Tbase_Cam = Tbase_eff*(Teff_inst)*(inv(TCam_Inst));
[Rbase_Cam, tbase_Cam] = GetRt(Tbase_Cam);

tCam_trocar = (inv(Tbase_Cam))*[tbase_trocar; 1];
tCam_trocar = tCam_trocar(1:3);

%% ====================================================================
%  ======================== CONTROL PROCESS ==========================
%  ====================================================================

% The q0 will be:
q = Tbase_eff(1:3,4);

tbase_target = [-500 -200 600]';

n_iter = 50;
final_n_iter = 0;

e = zeros(1,n_iter);
Iter = zeros(1,n_iter);

for k = 1:n_iter
    %%  Points in Camera frame
    
    % Current Position of the Instrument
    TCam_Inst = GetInstrumentPosition(flag);
    tCam_Inst = TCam_Inst(1:3,4);
    
    direction1 = (tCam_Inst - tCam_trocar)/(norm(tCam_Inst - tCam_trocar));
    tCam_eff = tCam_Inst + 350*direction1;
    
    %% ========== 2 - Error w.r.t. the Camera Frame
    
    direction2 = (tCam_target - tCam_trocar)/(norm(tCam_target - tCam_trocar));
    eCam = (tCam_target + 350*direction2) - tCam_eff;
    
    e(k) = norm(eCam);
    Iter(k) = k;

    if e(k) < 10
        final_n_iter = k;
        break;
    end

    %% ========== 3 - Error w.r.t. the Base Frame
    eBase = Rbase_Cam*eCam;

    dq = Jinv*(W*eBase);
    q = q + dq;

    position = MGD(q, tbase_trocar);
    MoveEffPosition(position(1:3,4));

    if k ~= 1
        DisplayImage
    end
end
% Move Robot

e = e(2:final_n_iter);
Iter = Iter(2:final_n_iter);

figure
plot(Iter,e)
grid on
%DisplayConfig








