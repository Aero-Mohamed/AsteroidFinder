% Reaction Wheel Model RWP015
% MaxTorque: 0.004Nm
% Mass: 0.13 kg for each reaction wheel
% volume: 42x42x19 mm

% Reaction Wheel Model RWP050
% MaxTorque: 0.007Nm
% Mass: 0.24 kg for each reaction wheel
% volume: 58x58x25 mm




% Mass and dimensions
mass_rw = 0.24;
height_rw = 25/1000;
radius_rw = 58/1000;
torque_rw_max = 0.007;


% orientation
n1 = [1; 0; 0];
n2 = [0; 1; 0];
n3 = [0; 0; 1];

% offset from c.g.
r1 = [40; 0; 0]/1000;
r2 = [0; 40; 0]/1000;
r3 = [0; 0; 40]/1000;


