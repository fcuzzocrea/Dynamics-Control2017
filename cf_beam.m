
function [K,M] = cf_beam(N_e)

% N_e : number of elements

% beam properties

E = 70e9;           % Young modulus (Pa)
rho = 2700;         % mass density (kg/m^3)
A = 0.02;           % cross section (m^2)
l = 10;             % length (m)
I_yy = 1.6667e-5;   % moment of inertia about y axis (m^4)
m = rho*A;          % mass per unit length (kg/m)
EJ = E*I_yy;        % bending stiffness (xz plane)

% number of nodes

N_n = N_e + 1;

 % number of FE degrees of freedom

N_dof = 2*N_n;     


% elemental stiffness and mass matrices

h_e = l/N_e;        % length of the beam element

K_e = (EJ/h_e^3)*[  12          6*h_e       -12         6*h_e   ;
                    6*h_e       4*h_e^2     -6*h_e      2*h_e^2 ;
                    -12         -6*h_e      12          -6*h_e  ;
                    6*h_e       2*h_e^2     -6*h_e      4*h_e^2 ];
                
M_e = (m*h_e/420)*[ 156         22*h_e      54          -13*h_e  ;
                    22*h_e      4*h_e^2     13*h_e      -3*h_e^2 ;
                    54          13*h_e      156         -22*h_e  ;
                    -13*h_e     -3*h_e^2    -22*h_e     4*h_e^2  ];

 
% stiffness and mass matrices (before applications of boundary conditions)
                
KK = zeros(N_dof, N_dof);
MM = zeros(N_dof, N_dof);

% assembly

for i = 1 : 2: 2*N_e-1
    
    KK(i : i + 3, i : i + 3) = KK(i : i + 3, i : i + 3) + K_e;
    MM(i : i + 3, i : i + 3) = MM(i : i + 3, i : i + 3) + M_e;

end

% stiffness and mass matrices for CF bcs

KK(1:2,:) = [];
KK(:,1:2) = [];
MM(1:2,:) = [];
MM(:,1:2) = [];

K = KK;
M = MM;


