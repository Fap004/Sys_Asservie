clc
clear
close all

%% Données

g = 9.81;

X = [0; 8; 15; 20; 25];

% Valeurs connues
Y_impose = [30; 19; 20; 16];

%% Matrice pour trouver le polynôme

phi1 = ones(size(X));
phi2 = X;
phi3 = X.^2;
phi4 = X.^3;
phi5 = X.^4;

M = [phi1 phi2 phi3 phi4 phi5];

Psi = M' * M;

%% Recherche de yf

y_test = linspace(10, 15, 1000);

best_y = 0;
min_deriv_abs = inf;

X_E = X(5);

for i = 1:length(y_test)

    % Valeur testée de yf
    Y = [Y_impose; y_test(i)];

    % Calcul des coefficients du polynôme
    Phi = M' * Y;
    A = inv(Psi) * Phi;

    % Dérivée du polynôme au point E
    deriv_E = A(2) ...
            + 2*A(3)*X_E ...
            + 3*A(4)*X_E^2 ...
            + 4*A(5)*X_E^3;

    % Recherche de la pente la plus proche de zéro
    if abs(deriv_E) < min_deriv_abs
        min_deriv_abs = abs(deriv_E);
        best_y = y_test(i);
    end

end

yf= best_y;
fprintf('yf (hauteur au point E) = %.4f m\n', best_y);

%% Interpolation finale

Y = [Y_impose; best_y];

Phi = M' * Y;

A = inv(Psi) * Phi;

a0 = A(1);
a1 = A(2);
a2 = A(3);
a3 = A(4);
a4 = A(5);

fprintf('\nPolynôme d''interpolation de la trajectoire :\n');
fprintf('y(x) = %.6f + %.6f*x + %.6f*x^2 + %.6f*x^3 + %.6f*x^4\n', ...
    a0, a1, a2, a3, a4);

%% Vérification de la dérivée au point E

deriv_E = a1 ...
        + 2*a2*X_E ...
        + 3*a3*X_E^2 ...
        + 4*a4*X_E^3;

fprintf('Pente au point E = %.6f\n', deriv_E);

%% Erreurs sur l'interpolation de la trajectoire

Y_calc = M * A;

% Erreur quadratique
Eq = (Y_calc - Y)' * (Y_calc - Y);

% Erreur RMS
err_rms = sqrt(mean((Y_calc - Y).^2));

% Coefficient R2
y_moy = mean(Y);

R2 = sum((Y_calc - y_moy).^2) / sum((Y - y_moy).^2);

fprintf('\nErreur quadratique E = %.4f\n', Eq);
fprintf('Erreur RMS = %.4f\n', err_rms);
fprintf('R2 = %.4f\n', R2);

%% Graphique

figure;

plot(X, Y, 'bo', 'MarkerSize', 8, 'LineWidth', 2);

hold on;

%% Trajectoire interpolée

XX = linspace(min(X), max(X), 500)';

FF = [ones(size(XX)), ...
      XX, ...
      XX.^2, ...
      XX.^3, ...
      XX.^4];

yy = FF * A;

plot(XX, yy, 'r-', 'LineWidth', 2);

grid on;

xlabel('x');
ylabel('y');

xlim([0 30]);
ylim([0 30]);

title('Trajectoire du participant');

legend('Points imposés', 'Trajectoire interpolée');

hold off;
%%  COEFFICIENT DE FRICTION
ouverture = [0 10 20 30 40 50 60 70 80 90 100];
mu = [0.87 0.78 0.71 0.61 0.62 0.51 0.51 0.49 0.46 0.48 0.46];
% On choisit un polynôme de degré 3
col1 = ones(length(ouverture), 1);
col2 = ouverture';
col3 = ouverture'.^2;
Cmu = [col1 col2 col3];
B = pinv(Cmu)*mu';

% Coefficients
b0 = B(1);
b1 = B(2);
b2 = B(3);

% Fonction du coefficient de friction
mu_calcule = b0 + b1*ouverture + b2*ouverture.^2;

% Erreur RMS
erreur = mu - mu_calcule;
N = length(mu);
erreur_RMS = sqrt(sum(erreur.^2)/N);

% Affichage
fprintf('-Debit eau\n');
fprintf('Polynôme de friction :\n');
fprintf('mu = %.6f + %.6f*u + %.8f*u^2\n',b0,b1,b2);
fprintf('Erreur RMS = %.6f\n',erreur_RMS);

% Graphique du lissage
u_plot = linspace(0,100,1000);
mu_plot = b0 + b1*u_plot + b2*u_plot.^2;

figure
plot(ouverture,mu,'o','MarkerSize',8);
hold on;
plot(u_plot,mu_plot,'LineWidth',2);
grid on;

xlabel('Ouverture de la valve (%)');
ylabel('\mu_f');
title('Coefficient de friction en fonction de l''ouverture');
legend('Données expérimentales','Approximation');
ylim([0 1])

% CHOIX DU COEFFICIENT DE FRICTION
v_min = 20/3.6;
v_max = 25/3.6;
yA = a0;

% mu correspondant à 20 km/h
mu_20 = (yA-yf-v_min^2/(2*g))/25;
mu_25 = (yA-yf-v_max^2/(2*g))/25;
mu_choisi = (mu_20 + mu_25)/2;

fprintf('Coefficient de friction choisi : %.4f\n',mu_choisi);

%  OUVERTURE DE LA VALVE
u = roots([b2 b1 b0-mu_choisi]);
ouverture_valve = u(u >= 0 & u <= 100);
fprintf('Ouverture de la valve = %.2f %%\n',ouverture_valve);

% VITESSE LE LONG DE LA TRAJECTOIRE
x_vitesse = linspace(0,25,1000);
H_x = a0 + a1*x_vitesse + a2*x_vitesse.^2 + a3*x_vitesse.^3 + a4*x_vitesse.^4;
v = sqrt(2*g*(yA - H_x - mu_choisi*x_vitesse));
v_kmh = v*3.6;

figure
plot(x_vitesse,v_kmh,'LineWidth',2);
grid on;

xlabel('x (m)');
ylabel('Vitesse (km/h)');
title('Vitesse du participant le long de la trajectoire');

%  VITESSE FINALE AU POINT E
vE = sqrt(2*g*(yA-yf-mu_choisi*25));

vE_kmh = vE*3.6;
fprintf('\nVitesse finale = %.2f km/h\n',vE_kmh);
% ERREUR SUR LA VITESSE CAUSÉE PAR L''ERREUR RMS

% mu de l'erreur RMS :
mu_plus = mu_choisi + erreur_RMS;
mu_moins = mu_choisi - erreur_RMS;

% Vitesses correspondantes
v_plus = sqrt(2*g*(yA-yf-mu_plus*25));
v_moins = sqrt(2*g*(yA-yf-mu_moins*25));

v_plus_kmh = v_plus*3.6;
v_moins_kmh = v_moins*3.6;

fprintf('Vitesse avec mu + RMS = %.2f km/h\n',v_plus_kmh);
fprintf('Vitesse avec mu - RMS = %.2f km/h\n',v_moins_kmh);

fprintf('Erreur sur la vitesse = +%.2f / -%.2f km/h\n\n', v_moins_kmh-vE_kmh, vE_kmh-v_plus_kmh);
%%  Requis pour design du ballon - mouse
% G1 - Ballon attrapé
mp = 80;
mb = 8;
vb = -1.0;
LT = 3;

vpb_apres = (mp*vE + mb*vb)/(mp+mb);
temps_G1 = LT/vpb_apres;
temps_min_minuterie = temps_G1 + 0.02;
% Affichage des résultats
fprintf('G1 - Ballon attrape\n');
fprintf('Vitesse du participant-ballon apres impact = %.2f km/h ou %.2f m/s\n',vpb_apres*3.6,vpb_apres);
fprintf('Temps pour quitter la trappe = %.3f s\n',temps_G1);
fprintf('La minuterie doit etre superieure a %.3f s\n\n',temps_min_minuterie);
% G2 - Ballon qui rebondit
e = 0.8;

difference = e*(vE-vb);
vp_apres = (mp*vE + mb*vb - mb*difference)/(mp+mb);
vb_apres = vp_apres + difference;
temps_G2 = LT/vp_apres;
temps_max_minuterie = temps_G2 - 0.02;

TempMinuterie = (temps_min_minuterie + temps_max_minuterie) / 2;

% Affichage des résultats
fprintf('G2 - Ballon qui rebondit\n');
fprintf('Vitesse du participant apres impact = %.2f km/h  ou %.2f m/s\n',vp_apres*3.6,vp_apres);
fprintf('Temps pour quitter la trappe = %.3f s\n',temps_G2);
fprintf('La minuterie doit etre inferieure a %.3f s\n\n',temps_max_minuterie);
fprintf('La minuterie est de %.3f s\n\n', TempMinuterie);
%% coussin - trampoline
m = mp + mb;
h0 = 5;
kc = 6000;
%resolution
a = kc/2;
b = -m*g;
c = -m*g*h0;
solutions = roots([a b c]);
dhc = solutions(solutions > 0);
% Affichage
fprintf('Deformation maximale du coussin = %.3f m\n\n',dhc);
%% Design bassin
b_eau = 47;     % coefficient hydrodynamique (kg/m)
kf = 0.95;      % facteur de flottabilité
h = 10;         % hauteur de chute (m)
ffin = 1.10;    % facteur de vitesse sécuritaire

% Vitesse limite dans l'eau (poids apparent = traînée : mp*g*(1-kf) = b*v^2)
v_lim = sqrt(mp*g*(1-kf)/b_eau);
fprintf('Vitesse limite = %.3f m/s\n',v_lim)
fprintf('Vitesse limite = %.3f km/h\n',v_lim*3.6)

% Vitesse à l'entrée dans l'eau
v0 = sqrt(2*g*h);
fprintf('Vitesse a l entree dans l eau = %.3f m/s\n',v0)

% Vitesse sécuritaire
v_securitaire = ffin*v_lim;
fprintf('Vitesse securitaire = %.3f m/s\n',v_securitaire)

% Profondeur sécuritaire (modèle linéarisé)
delta_z = (mp/(2*b_eau))*log((v0 - v_lim)/(v_securitaire - v_lim));
fprintf('Profondeur securitaire = %.3f m\n',delta_z)

% Vitesse le long de la profondeur
z_bassin = linspace(0,delta_z,1000);
v_bassin = v_lim + (v0 - v_lim)*exp(-(2*b_eau/mp)*z_bassin);

% Graphique
figure
plot(z_bassin,v_bassin,'LineWidth',2)
hold on
yline(v_lim,'--','Vitesse limite')
yline(v_securitaire,'--','Vitesse sécuritaire')
grid on12222
xlabel('Profondeur z (m)')
ylabel('Vitesse (m/s)')
title('Vitesse du participant dans le bassin')
legend('Vitesse','Vitesse limite','Vitesse sécuritaire')
%test 12-12<


laisncoawinclakwcxpAMCXPAlmx