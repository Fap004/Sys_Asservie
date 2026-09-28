clc; close all; clear

%%

tfin = 10;
M1 = 0.1;
M2 = 0.75;
g = 9.81;
k = 1;
mu = 0.05;
x1_ini = 0;
x2_ini = -0.15;
F = 0.3;
t_F = 2;

%% 

sim('atelier.slx')

%%

figure();
plot(t, x1, '--b');
hold on;
plot(t, x2, '-r');
title('Déplacement des wagons');
xlabel('Temps (s)');
ylabel('Déplacement (m)');
legend('Wagon1', 'Wagon2', Location='southoutside', Orientation='horizontal');

figure();
subplot(2,1,1);
plot(t, x1, '--b');
title('Déplacement des wagons');
xlabel('Temps (s)');
ylabel('Déplacement (m)');
subplot(2,1,2);
plot(t, x2, '-r');
xlabel('Temps (s)');
ylabel('Déplacement wagon 1 (m)');