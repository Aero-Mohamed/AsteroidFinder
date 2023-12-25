%%% Planet_Radius, M, G, Mu
planet;  %%% Radius, Mass, Gravity Const., Mu

%% Orbital Elements Init.
a_orbit = 600*1000+Planet_Radius;     % Meters
%%% caution to exceed the max eccentricity (1-Planet_radius/a_orbit)
e_orbit = 0.05;          % Eccentricity (0 < e < 1)
anomaly_orbit = 0;      % True anomaly (rads)
RAN_orbit = 50*pi/180;        % Right Ascension of aecending node (rads)
i_orbit = 30*pi/180;            % Inclination (rads) 0=> equatorial orbit
omega_orbit = 20*pi/180;% Argument of perigee (rads)

%% Calculations
T_orbit = 2*pi*sqrt(a_orbit^3/Mu); %% One Orbital Preiod 
anomaly_orbit_dot = 2*pi/T_orbit;  %% mean orbital motion
%%% Eccentric Anomaly
E_anomaly = 2* atan(sqrt((1-e_orbit)/(1+e_orbit)));
%%% position vector magnitude
r_mag = a_orbit*(1-e_orbit^2)/(1+e_orbit*cos(anomaly_orbit));
%%% velocity vector magnitude
r_mag_dot = sqrt(Mu*(2/r_mag - 1/a_orbit));
%%% r_p_frame: poistion vector in Perifoal Ref. Frame.
r_p_frame = [r_mag*cos(anomaly_orbit);
             r_mag*sin(anomaly_orbit);
             0];
        
%%% v_p_frame: velocity vector in Perifoal Ref. Frame.
%%% P_orbit: P = h^2/Mu (h: spesific angular momentum vector Mag.,)
P_orbit = a_orbit*(1-e_orbit^2);
v_p_frame = sqrt(Mu/P_orbit).*[-sin(anomaly_orbit); (cos(anomaly_orbit) + e_orbit); 0];

%% Transform from Perifocal Ref. Frame to ECI
s_omega = sin(omega_orbit); c_omega = cos(omega_orbit);
s_RAN = sin(RAN_orbit); c_RAN = cos(RAN_orbit);
s_i = sin(i_orbit); c_i = cos(i_orbit);
%%% Direct Cosine Matrix From Perifocal Frame to ECI;
Q_IP = [(-s_RAN*c_i*s_omega + c_RAN*c_omega), (-s_RAN*c_i*c_omega - c_RAN*s_omega), (s_RAN*s_i);
        (c_RAN*c_i*s_omega + s_RAN*c_omega),  (c_RAN*c_i*c_omega - s_RAN*s_omega),  (-c_RAN*s_i);
        (s_i*s_omega),                        (s_i*c_omega),                        (c_i)];
xyz_0 = Q_IP * r_p_frame;
xyz_dot_0 = Q_IP * v_p_frame;




