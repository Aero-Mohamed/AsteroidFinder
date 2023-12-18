global Inertia_3rw_cg

reactionWheelParams;

m = 2.6 + 3*mass_rw; %kilograms
% Inertia_3rw_cg can be found in reactionWheel.m
Inertia_satellite = [0.9 0 0; 0 0.9 0; 0 0 0.3];
I = Inertia_satellite  + Inertia_3rw_cg;

invI = inv(I);