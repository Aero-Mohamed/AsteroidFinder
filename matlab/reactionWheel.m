
reactionWheelParams; % All 3 Reaction Wheels Parameters

global Inertia_3rw_cg;
global Inertia_rw1_cg Inertia_rw2_cg Inertia_rw3_cg;
global Inertia_wr1_B Inertia_wr2_B Inertia_wr3_B;
global acc_rw_max maxSpeed;

Iyy_rw_by_mass = (1/12)*(3*radius_rw^2 + height_rw^2);
% Inertia on reaction wheel frame
Inertia_rw = mass_rw*[...
    radius_rw^2/2, 0, 0;
    0, Iyy_rw_by_mass, 0
    0, 0, Iyy_rw_by_mass];

% Max Reaction Wheel Acceleration & velocity
maxSpeed = 4000*2*pi/60; %rad/s % max speed 4000 rpm
acc_rw_max = torque_rw_max/Inertia_rw(1, 1);

% Rotation matrices
R1 = Rscrew(n1);
R2 = Rscrew(n2);
R3 = Rscrew(n3);

% Inertia of each reaction wheel on satellite body frame
Inertia_wr1_B = R1*Inertia_rw*R1';
Inertia_wr2_B = R2*Inertia_rw*R2';
Inertia_wr3_B = R3*Inertia_rw*R3';

% shift inertia to the center of mass of the satellite
skewr1 = skew(r1); skewr2 = skew(r2); skewr3 = skew(r3);
Inertia_rw1_cg = Inertia_wr1_B + mass_rw*(skewr1')*skewr1;
Inertia_rw2_cg = Inertia_wr2_B + mass_rw*(skewr2')*skewr2;
Inertia_rw3_cg = Inertia_wr3_B + mass_rw*(skewr3')*skewr3;

% Total 3 reaction wheels inertia on satellite body frame
% with respect to the c.g.
Inertia_3rw_cg = Inertia_rw1_cg + Inertia_rw2_cg + Inertia_rw3_cg;




