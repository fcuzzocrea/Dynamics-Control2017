clear all
close all
clc

%% RITZ MODEL FOR ROD

% Number of terms to be used for the approximation

N = 4;

% Exact solution for comparison
lambda_e = zeros(N,1);

for n = 1:N
    lambda_e(n,1) = ((2*n-1)^2*pi^2)/4;
end

% Compute generalized mass and stiffness matrices

M = zeros(N,N);
K = zeros(N,N);

for i = 1:N
    for j = 1:N
        M(i,j) = 1/(i+j+1);
        K(i,j) = (i*j)/(i+j-1);
    end
end

% Solve the eigenvalue problem

[U,Lambda] = eig(K,M);

% Vector of approximated frequency parameters

lambda = diag(Lambda);

x = linspace(0,1);
mode=zeros(N,length(x));

% Exact mode shapes 
for n = 1:N
    alphan= ((2*n-1)*pi)/(2);
    mode_shape_ex = sqrt(2)*sin(alphan*x);
end

% Approximated mode shapes
for n = 1:N
    for i = 1:length(x)
        for j = 1:N
            mode(n,i) =  mode(n,i) + (x(i)^j)*U(j,n);
        end
    end
end

% Plot shape functions 
figure(1)
grid on
hold on
for i=1:N
    plot(x,x.^i)
end

figure(2)
grid on;
hold on;
plot(x,mode_shape_ex,'k');
plot(x,mode(n,:),'k--');

%Comportamento strano dei polinomi ? Se N è pari il grafico viene
%specchiato rispetto ad x, se N è dispari viene corretto.
