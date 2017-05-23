%% ESERCITAZIONE 5
clear all
close all
clc

m1 = 200;
m2 = 150;
m3 = 100;

k1 = 1e9;
k2 = 0.5e9;

M_1 = m1*ones(5,1);
M_2 = m2*ones(5,1);
M_3 = m3*ones(10,1);
M = [0; M_1; M_2; M_3];
M = diag(M);

K_1 = [1 2 2 2 2 2 2 2 1.5 1 1 1 1 1 1 1 1 1 1 1 0.5].*1e9;
K_2 = -k1*ones(8,1);
K_3 = -k2*ones(12,1);
K_conc = [K_2; K_3];

K = diag(K_1) + diag(K_conc,1) + diag(K_conc,-1);

Mbb = M(1,1);
Mbi = M(1,2:end);
Mib = Mbi';
Mii = M([2:1:21],[2:1:21]);
Kbb = K(1,1);
Kbi = K(1,2:end);
Kib = Kbi';
Kii = K([2:1:21],[2:1:21]);

% Dozio check
eigvls=sqrt(eig(Kii,Mii))/(2*pi);
