clc;
close all;
clear;

% 1. Chargement et traitement des données expérimentales
load('data_1v_4-09_100hz.mat');
t_omega = tsimu(1:end-1) + diff(tsimu)/2;                  
t0 = tsimu(find(Vm > 0.5, 1));
idx = t_omega >= t0 + 0.5;                
K_omega = mean(omega_c(idx));  
K_rad = K_omega * (pi / 180);   
omega_63 = 0.63 * K_omega;
i2 = find(t_omega >= t0 & omega_c >= omega_63, 1);
i1 = max(1, i2 - 1);
t_63 = t_omega(i1) + (omega_63 - omega_c(i1)) * (t_omega(i2) - t_omega(i1)) / (omega_c(i2) - omega_c(i1));
tau = t_63 - t0;

% 2. Paramètres du système et du moteur
n_m  = 0.69;           
K_t  = 0.0076830;        
K_m  = 0.0076776;         
J_eq = 0.0017728;         
n_g  = 0.9000;       
K_g  = 70;              
a = n_m * K_t * n_g * K_g;                   
B_eq = J_eq / tau; 
R_m  = (n_m * K_t * n_g * K_g) / (K_rad * B_eq); 
m_s = 0.0030;          
J_s = 7.7028125e-7;        
r_s = 0.019625;  
g = 9.8100;       
J_m = 3.9001e-7;  
r_arm = 0.0254;          
L_plaque = 0.2750;            
J_c = -(K_g^2) * J_m + J_eq;
A_m = (n_m * K_t * n_g * K_g) / R_m;
K_bb = 5 / 7 * (g * r_arm) / L_plaque;

% --- Affichage des données expérimentales ---
% figure('Name', 'Données expérimentales', 'NumberTitle', 'off');
% 
% subplot(3, 1, 1);
% plot(tsimu, Vm, 'LineWidth', 1.2);
% grid on;
% xlabel('Temps (s)');
% ylabel('Tension (V)');
% title('Tension d''entrée V_m');
% 
% subplot(3, 1, 2);
% plot(t_omega, omega_c, 'r', 'LineWidth', 1.2);
% grid on;
% xlabel('Temps (s)');
% ylabel('Vitesse (deg/s)');
% title('Vitesse angulaire mesurée \omega_c');
% 
% subplot(3, 1, 3);
% plot(tsimu, servo, 'g', 'LineWidth', 1.2);
% grid on;
% xlabel('Temps (s)');
% ylabel('Signal Servo');
% title('Signal de position du servo');

% --- Signal d'entrée (Pulse) ---
% figure('Name', 'Signal d''entrée', 'NumberTitle', 'off');
% plot(t, pulse, 'b', 'LineWidth', 1.2);
% hold on;
% plot(t1, pulse1, 'r--', 'LineWidth', 1.2);
% grid on;
% xlabel('Temps (s)');
% ylabel('Amplitude (V)');
% title('Signal d''entrée (Pulse)');
% legend('Linéaire', 'Non linéaire', 'Location', 'best');
% 
% % --- Comparaison de la position \theta ---
% figure('Name', 'Simulation de theta', 'NumberTitle', 'off');
% plot(t, theta, 'b', 'LineWidth', 1.2);
% hold on;
% plot(t1, theta1, 'r--', 'LineWidth', 1.2);
% grid on;
% xlabel('Temps (s)');
% ylabel('Position \theta (rad)');
% title('Comparaison de la position angulaire \theta');
% legend('Linéaire', 'Non linéaire', 'Location', 'best');
% 
% % --- Comparaison de la position de la sphère x ---
% figure('Name', 'Simulation de x', 'NumberTitle', 'off');
% plot(t, x, 'b', 'LineWidth', 1.2);
% hold on;
% plot(t1, x1, 'r--', 'LineWidth', 1.2);
% grid on;
% xlabel('Temps (s)');
% ylabel('Position x (m)');
% title('Comparaison de la position de la sphère x');
% legend('Linéaire', 'Non linéaire', 'Location', 'best');

% Matrices
A = [ 0, 1, 0, 0;
      0, 0, K_bb, 0;
      0, 0, 0, 1;
      0, 0, 0, -B_eq/J_eq ];
B = [0; 0; 0; A_m/J_eq];
C = [1, 0, 0, 0; 0, 0, 1, 0];
D = [0; 0];

% % --- Boucle interne ---
% Gcm = tf([A_m],[J_eq B_eq 0]);
% zeta = 0.8;
% K_int = B_eq^2 / (4*zeta^2 * J_eq * A_m)
% 
% FTBF_cm = feedback(K_int*Gcm,1);
% 
% % figure;
% % hold on;
% % rlocus(Gcm);
% % rlocus(FTBF_cm)
% % % [Kcrit,Pcrit] = rlocfind(Gcm)
% % grid on;
% 
% A_int = A - B*K_int*C(2,:);
% B_int = B*K_int;
% C_int = C;
% D_int = D;
% sys_int = ss(A_int,B_int,C_int,D_int);
% 
% G_int   = tf(sys_int);
% Gx_int  = G_int(1);    % de θ_cd à x
% Gth_int = G_int(2);    % de θ_cd à θ_c
% 
% % --- Boucle externe --- compensateur temporel
% 
% Mp= 0.05;                               
% tr = 2;
% tp = 3;
% ts = 4;
% 
% wn_mod = 0.535; % 0.535(0.90) ou 0.6(0.95)
% zeta_modif = 0.21; % +0.26(0.95) ou +0.21(0.90);
% 
% phi   = atan(-pi/log(Mp));
% zeta_ext  = cos(phi) + zeta_modif
% wn_tr   = (1 + 1.1*zeta_ext + 1.4*zeta_ext^2)/tr;
% wn_ts = 4/(zeta_ext*ts);
% wn_tp = pi/(tp*sqrt(1-zeta_ext^2));
% 
% wn = wn_ts + wn_mod;
% 
% 
% sd    = -zeta_ext*wn + 1j*wn*sqrt(1 - zeta_ext^2);
% 
% [num,den] = tfdata(Gx_int,'v');
% 
% Gsd    = polyval(num, sd)/polyval(den, sd);
% dphi   = mod(-180 - rad2deg(angle(Gsd)), 360);
% alpha = rad2deg(angle(sd));                
% 
% za = real(sd) - imag(sd)/tand((alpha + dphi)/2);
% pa = real(sd) - imag(sd)/tand((alpha - dphi)/2);
% ka = 1/abs((sd - za)/(sd - pa)*Gsd);
% Ga = ka*tf([1 -za], [1 -pa]);
% 
% FTBO = Ga*Gx_int;
% FTBF = feedback(FTBO, 1);
% 
% figure;
% hold on;
% rlocus(Gx_int)
% rlocus(FTBO);
% grid on;

% %% Compensateur freq
% PM_des = 45;
% BW_des = 2.3; 
% marge_PM = 0;
% wn2_modif = 0.1474; %2.5 (0.1474) 2.6 (0.2211)
% 
% zeta2 = 1/2 * sqrt(tand(PM_des)*sind(PM_des));
% wn2 = (BW_des /sqrt((1-2*zeta2^2)+sqrt(4*zeta2^4-4*zeta2^2+2))) + wn2_modif
% wg_des = wn2*sqrt(sqrt(1+4*zeta2^4)-2*zeta2^2);
% 
% [mag,pha] = bode(Gx_int,wg_des);
% Kdes = 1/mag;
% PM2   = pha + 180;   
% dphi2 = PM_des - PM2 + marge_PM ;
% 
% alpha = (1 - sind(dphi2))/(1 + sind(dphi2));
% T     = 1/(wg_des*sqrt(alpha));
% Ka2    = Kdes/sqrt(alpha);
% Ga2    = Ka2*tf([1 1/T], [1 1/(alpha*T)]);
% 
% [num_a2, den_a2] = tfdata(Ga2, 'v');
% 
% FTBO2 = Ga2*Gx_int;
% FTBF2 = feedback(FTBO2,1);
% BW = bandwidth(FTBF2)
% [GM2, PM2_m] = margin(FTBO2);
% GM2_dB = 20*log10(GM2)            % plus grand que 12 dB
% PM2_m                             % plus grand que 45 deg
% 
% %% Choix du compensateur
% compensateur = 'temp';
% 
% switch compensateur
%     case 'temp'
%         Gc = Ga;
%     case 'freq'
%         Gc = Ga2;
% end
% [num_a, den_a] = tfdata(Gc, 'v');
% 
% flg_in  = 0;
% retard  = 0;
% bruit_P = 0;
% 
% % Simulation non linéaire (Simulink)
% sim('SM1_d.slx');
% t_nl    = t_nl(:);
% x_nl    = x_nl(:);
% thcd_nl = thcd_nl(:);
% Vm_nl   = Vm_bl(:);
% 
% % Vérification - critères d'acceptation
% x_d = 0.05;
% t = (0:0.001:10)';
% u = x_d*ones(size(t));
% 
% FTBO_c = Gc*Gx_int;
% FTBF_c = feedback(FTBO_c, 1);
% y = lsim(FTBF_c, u, t);
% 
% info = stepinfo(y, t, x_d);
% info_nl = stepinfo(x_nl, t_nl, x_d);
% [GM, PM] = margin(FTBO_c);
% 
% Mp_nl = info_nl.Overshoot           % temp : <  30 %  freq : aucun
% ts_nl = info_nl.SettlingTime        % temp : < 5 s   freq : < 4.5 s
% PM                                  % temp : > 45 deg | freq : > 45 deg
% GM_dB = 20*log10(GM)                % temp : > 10 dB  | freq : > 12 dB
% BW    = bandwidth(FTBF_c)           % temp : aucun     | freq : > 2.5 rad/s
% % Saturation
% thcd_max_nl = max(abs(thcd_nl))*180/pi   % < 56 deg
% Vm_max_nl   = max(abs(Vm_nl))            % < 10 V
% 
% 
% figure;
% plot(t, y, 'b', t_nl, x_nl, 'r--', 'LineWidth', 1.2); hold on;
% yline(1.02*x_d, 'k:'); yline(0.98*x_d, 'k:');
% grid on; xlabel('Temps (s)'); ylabel('Position x (m)');
% legend('Linéaire (lsim)', 'Non linéaire (Simulink)', '\pm2 %', 'Location', 'best');
% title(['Réponse à un échelon de 5 cm - compensateur ', compensateur]);
% %% Robustesse : délai dans la boucle externe
% %[~, PM_t, ~, wg_t] = margin(Ga*Gx_int);
% [~, PM_f, ~, wg_f] = margin(Ga2*Gx_int);
% %DM_t = deg2rad(PM_t)/wg_t                     % délai max théorique (s), temporel
% DM_f = deg2rad(PM_f)/wg_f                     % délai max théorique (s), fréquentiel
% 
% bruit_P = 0;
% for retard = [0 0.02 0.05 0.1]
%     sim('SM1_d.slx');
%     t_nl = t_nl(:); x_nl = x_nl(:); thcd_nl = thcd_nl(:); Vm_nl = Vm_bl(:);
%     info_r = stepinfo(x_nl, t_nl, 0.05);
%     fprintf('retard = %.2f s | Mp = %5.1f %% | ts = %5.2f s | thcd max = %5.1f deg | Vm max = %5.2f V\n', ...
%         retard, info_r.Overshoot, info_r.SettlingTime, max(abs(thcd_nl))*180/pi, max(abs(Vm_nl)));
% end
% retard = 0;
% 
% %% Robustesse : bruit sur la mesure de x
% retard = 0;
% for bruit_P = [0 1e-8 1e-7 1e-6]
%     sim('SM1_d.slx');
%     t_nl = t_nl(:); x_nl = x_nl(:); thcd_nl = thcd_nl(:); Vm_nl = Vm_bl(:);
%     fprintf('sigma = %4.1f mm | thcd max = %5.1f deg | Vm max = %5.2f V | erreur RMS finale = %.2f mm\n', ...
%         sqrt(bruit_P/0.01)*1000, max(abs(thcd_nl))*180/pi, max(abs(Vm_nl)), ...
%         rms(x_nl(t_nl > 6) - 0.05)*1000);
% end
% bruit_P = 0;
% 
% %% sauvegarde
% 
% % save('param_freq_BW26.mat', 'K_int', 'num_a', 'den_a'); %load('param_freq_BW26.mat')
% % save('param_freq_BW25.mat', 'K_int', 'num_a', 'den_a'); %load('param_freq_BW25.mat')
% save('param_temp_0.90.mat','K_int', 'num_a', 'den_a'); %load('param_temp_0.90.mat')
% % save('param_temp_0.95.mat','K_int', 'num_a', 'den_a'); %load('param_temp_0.95.mat')