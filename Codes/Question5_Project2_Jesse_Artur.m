%% ============ PROJECT2 - 3D Medical Registration ==============
%% =================== Jesse and Artur ========================
clear all; close all; clc

%% ==============>   Question 5

% Set the Initial Configuration
InitConfig();
flag = 0;

%% Variance Propagation

Ex = 3*eye(3);

% Tworld_eff
Tbase_eff = GetRobotCurrentPosition(flag);
tbase_eff  = Tbase_eff(1:3,4);
%tbase_eff = [-305.0199; -70.0133; 859.9735]
tworld_eff = Tworld_base*[tbase_eff; 1];
tworld_eff = tworld_eff(1:3);

%Tworld_trocar
tbase_trocar = [-350; -100; 800];
tworld_trocar = Tworld_base*[tbase_trocar; 1];
tworld_trocar = tworld_trocar(1:3);

%% Computing Jacobian
syms x y z 
X = [x y z]';
Y = tworld_eff + 350*((X - tworld_eff)/(norm(X - tworld_eff)));

J = jacobian(Y);

x = tworld_trocar(1);
y = tworld_trocar(2);
z = tworld_trocar(3);

disp('The Jacobian Matrix')
J = double(subs(J))

%Ey
Ey = J*Ex*J'

%% Many Noisy Configuration
% Let consider tworld_trocar as the right value

Y = tworld_eff + 350*((tworld_trocar - tworld_eff)/(norm(tworld_trocar - tworld_eff)));
Ey_noisy = 0;

n = 100000-1;
for k=1:n
    tworld_trocar_k = tworld_trocar + 3*(-1 -rand([1,3]).*Ex);

    Yk = tworld_eff + 350*((mean(tworld_trocar_k,2) - tworld_eff)/(norm(mean(tworld_trocar_k,2) - tworld_eff)));

    Ey_noisy = (Yk - Y)*((Yk - Y)')  + Ey_noisy;
end

Ey_noisy = Ey_noisy/n

Fro_Norm = norm(Ey_noisy - Ey,'fro')

%% Plotting the Elipsoids

%% Ey Variance Propagation
center = tworld_trocar;
%rot_mat = thetau2r(Ey_noisy*n);
rot_mat = Ex;

[U,S,V] = svd(Ey_noisy);
xaxis = S(1,1)*V(:,1);
yaxis = S(2,2)*V(:,2);
zaxis = S(3,3)*V(:,3);

axe_lengthes = [xaxis  yaxis zaxis];
scale = 1;

plot_ellipsoid(center, rot_mat, axe_lengthes, scale)
xlabel('x')
ylabel('y')


