clc;
close all;
clear;

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

sim('SM1_Linear.slx');

%%

figure();

subplot(2, 1, 1);
plot(t, omega, '--b', 'LineWidth', 2);
hold on;
plot(t, theta, '-r', 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Amplitude');
title('Moteur-Engrenage-Charge');

subplot(2, 1, 2);
plot(t, v_x, '--b', 'LineWidth', 2);
hold on;
plot(t, x, '-r', 'LineWidth', 2);
grid on;
xlabel('Temps (s)');
ylabel('Amplitude');
title('Charge-Sphère');