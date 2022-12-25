%% Define all the constants 
clear all; close all; clc;
Ls = 2.47*1e-6;
Rs = 5*1e-3;
Lc = 100*1e-9;
c = 95*1e-6;
Rc = 2*1e-3;
Lcab = 0.0;
Rcab = 0.0;
Lp = 8.21e-6;
Rp = 5*1e-3;
L = 0.9*1e-6; % in vacuum: 1.85e-6
R = 100*1e-3; % in vacuum: 20e-3 
Vspa = 400;
a = 1/(Ls+Lc);
b = 1/(Lcab+Lp+Lc);
d= a*(Rp+Rc); % Rdelta should be here?
f = a*Rc;
g = a*Lc;
p = Rp/(L+Lp);
h = b*Rc;
k = b*(Rcab+Rp+Rc);
l = b*Rp;
m = b*Lc;
n = b*Lp;
q = (R+Rp)/(L+Lp);
v = Lp/(L+Lp);
%% Define the matrics for xdot = Ax + Bu, y = Cx
A = 1/(m*g-n*v)*[g*h-d*(1-n*v), f*(1-n*v)-g*(k-n*p), g*(l-n*q), g*b-a*(1-n*v);
    h-m*d, m*f-k + n*p, l-n*q, b-m*a;
    v*(h-m*d), p*(1-m*g)+v*(m*f-k), v*l-q*(1-m*g), v*(b-m*a);
    (1-m*g-n*v)/c, -(1-m*g-n*v)/c, 0, 0];
B = 1/(m*g-n*v)*[a*(1-n*v),m*a,v*m*a,0]';
C = [0,0,1,0];
D = zeros(size(C,1),size(B,2));
%% Check observability and controllability
obsv(A,C)
ctrb(A,B)
det(obsv(A,C))
sys = ss(A,B,C,D);
%det(gram(sys,'o'))
sysFullOutput = ss(A,B,eye(4),zeros(4,size(B,2)));  % system with full state output, disturbance, no noise
%% Build LQR controller
Q = eye(4);
R = 0.0001;
K = lqr(A,B,Q,R);
%% Add in sensor noise and disturbances
Vd = 0.1*eye(4);
Vn = 1;
Bf = [B Vd 0*B];
sysC = ss(A,Bf,C,[0 0 0 0 0 Vn]);
%% Build Kalman filter
[Kf,P,E] = lqe(A,Vd,C,Vd,Vn);
sysKf = ss(A-Kf*C,[B Kf],eye(4),0*[B Kf]);
%% Solve the system with the Kalman filter
%t = 0:dt:50;
%ud = randn(4,size(t,2));
%un = rand(size(t));
%u = 0*t;
%u(200:400) = Vspa;
%u(600:800) = -Vspa;
%u_augmented = [u; Vd*Vd*ud; un];
%[y,t] = lsim(sysC,u_augmented,t);
%plot(t,y)