clear all,close all, clc

sym ('a');
sym ('am');
sym ('ak');

A = [cos(a) sin(a)  cosh(a) sin(a) am/ak*a^4 -1;
    1       0           1       0             0;
    0       0           0       1             0;
    cos(a)  sin(a)  -cosh(a) -sinh(a)         0;
    sin(a) -cos(a)   sinh(a)  cosh(a)     -am*a]

chareq = simplify(det(A));

% Analisi parametrica fai una mesh di vari valori di ak e am

% Function
fun = @(a, ak, am) (ak-am*a.^4).*(1+cos(a))+ak*am*a.*(cos(a).*sinh(a)-sin(a).*cosh(a));
a = linsplace(1,10,1000);
y = fun(a,alphak,alpham);

%Plot
plot (a,y)

% Solve with FSolve
options = optimoption('fsolve');
options.Algorithm = 'levenber-marquadt';
option.Displa='iter';
option.FiniteDifferenceType = 'central';
option.StepTolerance = 1e-12;
sol = fsolve(@(a) fun(a,alphak,alpham), 1.4:7,options);