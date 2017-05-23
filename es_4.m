%% ESERCITAZIONE 4

clear all 
close all
clc

% Define properties of the bodyes
m=60;
k=2.5;
M = [m 0 0; 0 m 0; 0 0 m];
K = [k -k 0; -k 2*k -k; 0 -k k];

Lf = [1 0 0]';
Lu = [0 0 1];

% Frequency domain
omg = linspace(0.01,3,1000);

% Direct recovery of FRF

H = [];

for i = 1:length(omg)    
    % Compute FRF by directly solving 
    % the linear system
    H(i) = Lu*((K-omg(i)^2*M)\Lf);
end

% Recovery of FRF using mode displacements method 

[L,eigv]=eig(K,M);

modal_mass = L'*M*L;
modal_stiffness = L'*K*L;
modal_force = L'*Lf;

H_MDM = zeros(length(omg),1);

eigv = diag(eigv);

% Non è il modo più efficiente di farlo
for i = 1:length(omg)
        H_MDM_1 = L(:,1)*(1/(1*(eigv(1)-omg(i)^2)))*L(:,1)';
        H_MDM_2= L(:,2)*(1/(1*(eigv(2)-omg(i)^2)))*L(:,2)';
        H_MDM_3 = L(:,3)*(1/(1*(eigv(3)-omg(i)^2)))*L(:,3)';
        H_MDM(i)=Lu*(H_MDM_1+H_MDM_2+H_MDM_3)*Lf;
end

% Recovery of FRF using mode acceleration method 



% Direct recovery
figure(1)
loglog(omg,abs(H));
grid on
figure(2)
semilogx(omg,mag2db(abs(H)));
grid on

% Recovery using modes
figure(3)
loglog(omg,abs(H_MDM));
grid on
figure(4)
semilogx(omg,mag2db(abs(H_MDM)));
grid on

