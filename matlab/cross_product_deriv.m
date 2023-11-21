%% I did use this file to compute derivatives symbolically
%  Which was required in the formulation of Kalman Filter 
%  jacobian matrix

%%
syms p q r dp dq dr real
syms I11 I22 I33 real

% Define vectors A and B
omega = [p; q; r];
delta_omega = [dp; dq; dr];

% Define the inertia matrix I
Inertia__ = [I11, 0, 0; 0, I22, 0; 0, 0, I33];

% Calculate the cross product A x C
term1 = cross(omega, Inertia__*delta_omega);
term2 = cross(delta_omega, Inertia__*omega);
term3 = cross(delta_omega, Inertia__*delta_omega);


% Calculate the derivative
derivative1 = jacobian(term1, [dp, dq, dr]);
derivative2 = jacobian(term2, [dp, dq, dr]);
derivative3 = jacobian(term3, [dp, dq, dr]);

delta_omegaDot_delta_omega = -derivative1 - derivative2 - derivative3;

%%
syms p q r dp dq dr real
syms I11 I22 I33 real
syms dq1 dq2 dq3 real

omega = [p; q; r];
delta_omega = [dp; dq; dr];
dq = [dq1; dq2; dq3];

term1 = cross(0.5*(-delta_omega -2*omega), dq);
term2 = cross(0.5*(-delta_omega -2*omega), dq) + delta_omega;

derivative1 = jacobian(term1, [dq1, dq2, dq3]);
derivative2 = jacobian(term2, [p, q, r]);





