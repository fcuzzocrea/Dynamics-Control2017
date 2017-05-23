%% ESERCITAZIONE 6

% Compute the tip displacemente when the tip force applied on the beam is
% a unit step using Newmark method to integrate directly the equation of 
% motion. Write the equation in chpt 13, and implement the method to have
% the time history. Then compare this solution with a modal solution.
% Select time step and beta and gamma. Integrate from t=0 to t=5s
% U0, V0, = 0. 

close all
clear
clc

N_e = 10;

[M,K] = cf_beam(N_e);
f0 = 1;              
f=zeros(length(M),1);
f(end-1) = f0;
f=f.*1e3;

u0 = zeros(length(M),1);
v0 = zeros(length(M),1);
a0 = (M\(f-K*u0));   %?    

t=(0:0.05:500);
dt=t(2)-t(1);
beta = 1/6;
gamma = 1/2;

%Constants
a_0 = 1/(beta*dt^2);
a_1 = gamma/(beta*dt);
a_2 = 1/(beta*dt);
a_3 = (1/(2*beta))-1;
a_4 = (gamma/beta)-1;
a_5 = dt*((gamma/(2*beta))-1);
a_6 = dt*(1-gamma);
a_7 = dt*gamma;

u = zeros(length(u0),length(t));
v = zeros(length(u0),length(t));
a = zeros(length(u0),length(t));

u(:,1) = u0;
v(:,1) = v0;
a(:,1) = a0;

for i=1:(length(t)-1)
    if i>1
        f(end-1) = 0;
    end
    u(:,i+1) = (a_0*M + K)\(f + M*(a_0*u(:,i) + a_2*v(:,i) + a_3*a(:,i)));
    v(:,i+1) = v(:,i) + a_6*a(:,i) + a_7*a(:,i+1);
    a(:,i+1) = a_0*(u(:,i+1)-u(:,i)) - a_2*v(:,i) - a_3*a(:,i);  
end

figure(1)
plot(t,u(end-1,:))
figure(2)
plot(t,v(end-1,:))
figure(3)
plot(t,a(end-1,:))
