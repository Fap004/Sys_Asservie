clc;
close all;
clear;

addpath(genpath(pwd));
%%

% sphère
m_s = 0.0030;          
J_s = 7.7028125e-7;        
r_s = 0.019625;           

% gravité
g = 9.8100;             

% partie électrique du moteur
R_m = 1; % à déterminer             
K_m = 0.0076776;         
K_t = 0.0076830;        
n_m = 0.69;           

% partie mécanique du moteur, la charge et engrenage
J_m = 3.9001e-7;         
J_eq = 0.0017728;         
B_eq = 1; % à déterminer           
n_g = 0.9000;       
K_g = 70;              

% poutre
r_arm = 0.0254;          
L_plaque = 0.4254;            

% calculs 
J_c = (K_g^2) * J_m - J_eq;
A_m = (n_m * K_t * n_g * K_g)/R_m;
K_bb = 5 / 7 * (g * r_arm) / L_plaque;

%% 

sim('SM1.slx');

%%

figure();

subplot(2, 1, 1);
plot(t, omega, '--b', 'LineWidth', 2);
hold on;
plot(t, theta, '-r', 'LineWidth', 2);
plot(t, Vm, '-m', 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Amplitude');
legend('Vitesse angulaire', 'Angle','Tension');
title('Moteur-Engrenage-Charge');

subplot(2, 1, 2);
plot(t, v_x, '--b', 'LineWidth', 2);
hold on;
plot(t, x, '-r', 'LineWidth', 2);
plot(t, Vm, '-m', 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Amplitude');
legend({'Vitesse', 'Position','Tension'}, 'Location', 'best');
title('Charge-Sphère');

%% Matrices
A = [ 0, 1, 0, 0;
      0, 0, K_bb, 0;
      0, 0, 0, 1;
      0, 0, 0, -B_eq/J_eq ];

B = [0; 0; 0; A_m/J_eq];

C = [1, 0, 0, 0; 0, 0, 1, 0];

D = [0; 0];

sys = ss(A, B, C, D);

p = pole(sys);

disp('Pôles avec matrices:')
disp(p);

%% Fonction de transfert
tau = J_eq / B_eq;
K_cm = A_m / B_eq;

num_gcm = [K_cm];
den_gcm = [tau, 1, 0]; 
G_cm = tf(num_gcm, den_gcm);

num_gsc = [K_bb];
den_gsc = [1, 0, 0]; 
G_sc = tf(num_gsc, den_gsc);

G = series(G_cm, G_sc);
poles = pole(G);

disp('Pôles avec fonctions de transfert en série:');
disp(poles);

% Pôles individuels de chaque sous-système
poles_gcm = pole(G_cm);
poles_gsc = pole(G_sc);

disp('--- Pôles individuels ---');
disp('Pôles de G_cm(s) :');
disp(poles_gcm);
disp('Pôles de G_sc(s) :');
disp(poles_gsc);
 

