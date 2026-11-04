clc;
close all;
clear;

load('data_1v_4-09_100hz.mat');
t_omega = tsimu(1:end-1);                  
t0 = tsimu(find(Vm > 0.5, 1));

idx = t_omega >= t0;                
K_omega = mean(omega_c(idx))       
omega_63 = 0.632 * K_omega;
idx_63 = find(t_omega >= t0 & omega_c >= omega_63, 1);
tau = t_omega(idx_63) - t0          
K_rad = K_omega * (pi / 180);   

G = tf(K_omega, [tau, 1]);
t_sim_filtre = t_omega(t_omega >= t0) - t0;
omega_c_filtre = lsim(G, ones(size(t_sim_filtre)), t_sim_filtre);

figure();
plot(t_omega, omega_c);
hold on;
plot(t_sim_filtre + t0, omega_c_filtre);
grid on;
xlabel('Temps (s)');
ylabel('Vitesse angulaire (deg/s)');

n_m  = 0.69;           
K_t  = 0.0076830;        
K_m  = 0.0076776;         
J_eq = 0.0017728;         
n_g  = 0.9000;       
K_g  = 70;              
a = n_g * K_g * n_m * K_t;                   
b = n_g * (K_g^2) * n_m * K_t * K_m;         
R_m  = (a * tau) / (K_rad * J_eq)
B_eq = (a/K_rad - b) / R_m

m_s = 0.0030;          
J_s = 7.7028125e-7;        
r_s = 0.019625;  
g = 9.8100;       
J_m = 3.9001e-7;  
r_arm = 0.0254;          
L_plaque = 0.4254;            

J_c = (K_g^2) * J_m - J_eq;
A_m = (n_m * K_t * n_g * K_g) / R_m;
K_bb = 5 / 7 * (g * r_arm) / L_plaque;

sim('SM1_d.slx');

figure();
subplot(2, 1, 1);
plot(t, omega, '--b');
hold on;
plot(t, theta, '-r');
grid on;
xlabel('Temps (s)');
ylabel('Amplitude');
title('Moteur-Engrenage-Charge');

subplot(2, 1, 2);
plot(t, v_x, '--b');
hold on;
plot(t, x, '-r');
grid on;
xlabel('Temps (s)');
ylabel('Amplitude');
title('Charge-Sphère');

% Matrices
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

tau_sys = J_eq / B_eq;
K_cm = A_m / B_eq;

num_gcm = [K_cm];
den_gcm = [tau_sys, 1, 0]; 
G_cm = tf(num_gcm, den_gcm);

num_gsc = [K_bb];
den_gsc = [1, 0, 0]; 
G_sc = tf(num_gsc, den_gsc);

G_total = series(G_cm, G_sc);
poles = pole(G_total)

% Pôles individuels de chaque sous-système
poles_gcm = pole(G_cm);
poles_gsc = pole(G_sc);

disp('--- Pôles individuels ---');
disp('Pôles de G_cm(s) :');
disp(poles_gcm);
disp('Pôles de G_sc(s) :');
disp(poles_gsc);

% Valeurs propres
[num_G, den_G] = tfdata(G_total, 'v'); 
valeurs_propres_tf = roots(den_G);
disp('Valeurs propres avec roots:');
disp(valeurs_propres_tf);

% Conditions d'équilibre
V_m_eq = 0;                 
x_eq = 0;                  
v_x_eq = 0;                 
theta_c_eq = 0;             
omega_c_eq = 0;

X_eq = [x_eq; v_x_eq; theta_c_eq; omega_c_eq];
Y_eq = [x_eq; theta_c_eq];

% Graphique de la réponse
t2 = (0:0.01:20)';
r_echelon  = ones(size(t2));
r_rampe    = t2;
r_parabole = t2.^2 / 2;

% Calcul des réponses
y_echelon  = lsim(G_total, r_echelon, t2);
y_rampe    = lsim(G_total, r_rampe, t2);
y_parabole = lsim(G_total, r_parabole, t2);

figure;

% Échelon
subplot(3, 1, 1);
plot(t2, r_echelon, '--', t2, y_echelon);
title('Échelon');
ylabel('Position [rad]');
legend('Consigne (r)', 'Réponse (y)', 'Location', 'Best');
grid on;

% Rampe
subplot(3, 1, 2);
plot(t2, r_rampe, '--', t2, y_rampe);
title('Rampe');
ylabel('Position [rad]');
legend('Consigne (r)', 'Réponse (y)', 'Location', 'Best');
grid on;

% Parabole
subplot(3, 1, 3);
plot(t2, r_parabole, '--', t2, y_parabole);
title('Parabole');
ylabel('Position [rad]');
xlabel('Temps [s]');
legend('Consigne (r)', 'Réponse (y)', 'Location', 'Best');
grid on;

figure();
plot(t, x, '-r', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Position x (m)');
title('Simulation non linéaire avec saturation de la position');
ylim([-0.25 0.25]);
 