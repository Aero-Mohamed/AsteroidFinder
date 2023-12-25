%% Notes
% -  I used the International Geomagnetic Reference Field (IGRF) Model
% is an internationally agreed upon mathematical model of 
% the Earth's magnetic field. File can be found on igrf/ folder
% and can be downloaded from the following link
% https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model?s_tid=mwa_osa_a
addpath 'igrf/' % International Geomagnetic Reference Field Model

%% Initialization & Planet Parameters
% profile on
clear; clc; close all;  tic;
global BB; %% Mag Field Body Frame
global magFieldTimeStep currentTime lastMagFieldTime;
global BB_meas pqr_meas;
global w123_dot_rw_current; %% Reaction Wheels Command Signal
orbitElement;       %% Determine Orbital inital Conditions from Elements
reactionWheel;      %% Reaction Wheels Inertia Calculations
satellite_inertia;  %% Mass, Inertia
filter_data;        %% gyroCovariance, Q, R, H_jac, P_0 for kalman filter

%% Initial Conditions 

% State Vector Initialization
phi_theta_psi_0 = [0, 0, 0*pi/180];
q0123_0 = EulerAngles2Quaternions(phi_theta_psi_0)';
pqr_0 = [3*pi/180, 5*pi/180, -0*pi/180];
w_reactionWheel_0 = [0, 0, 0]; % reaction wheel anglure velocity
initail_state = [xyz_0', xyz_dot_0', q0123_0, pqr_0, w_reactionWheel_0]';

%% Time Frame of an full orbit
timeStep = 1;
t_span = 0:timeStep:T_orbit;
% Calculate magField from IGRF model every n step
magFieldTimeStep = 100; lastMagFieldTime = 0;
% Get Sensor Date every n step
sensorTimeStep = 1; lastSensorTime = 0;

%% Locating matrices sizes
state = NaN(length(t_span), length(initail_state));
state(1, :) = initail_state';
% Store Magnitic Field readings from IGRF model & Sensor
BB_out = NaN(length(t_span), 3);
BB_meas_out = NaN(length(t_span), 3);
% Store pqr reading from gyro sensor
pqr_meas_out = NaN(length(t_span), 3);

% Kalman Filter Covariance Matrix
P_cov = NaN(9, 9, length(t_span)); % multiple diminsions 
P_cov(:, :, 1) = P_0;

% For External Forces and Moments
dF = [0, 0, 0]';
dM = [0, 0, 0]';

% initiale reaction wheel angular acceleration commands
w123_dot_rw = NaN(length(t_span), 3);
w123_dot_rw(1, :) = [0, 0, 0]; 
prevErr = 0; integral = 0;


% Simulation Start
for i = 1:length(t_span)-1

    currentTime = t_span(i);
    w123_dot_rw_current = w123_dot_rw(i, :)';
    
    % RK4 for non linear equations of motion
    state(i+1, :) = SatelliteRK4(timeStep, state(i, :), dF, dM);
    BB_out(i, :) = BB; % Get Magnetic Field from IGRF model
    
    %%% Compute Kalman Filter F jacobian matrix
    currentState = state(i+1, :);
    prevState = state(i, :);
    F = kalmanSystemJacobian(timeStep, currentState, prevState, I, invI);
    P_cov(:, :, i+1) = F*P_cov(:, :, i)*F' + Q;
    
    % Sensor Measurements | Assume measurements are available every second.
    if (currentTime >= lastSensorTime)
        pqr = currentState(11:13);
        [BB_meas, pqr_meas] = Sensors(BB, pqr);
        BB_meas_out(i, :) = BB_meas;
        pqr_meas_out(i, :) = pqr_meas;
        lastSensorTime = currentTime + sensorTimeStep;
        
        % let the state be [ q234 pqr TxTyTz] 
        % T is disturbance Torque, 
        % as suggested by the given Diploma Thesis
        state_estimate = [state(i+1, 8:10), state(i+1, 11:13), [0 0 0]];
        [corrected_state, p_cov_corrected] = ES_EKF(varGyro, P_cov(:, :, i+1), state_estimate, pqr_meas_out(i, :));
        P_cov(:, :, i+1) = p_cov_corrected;
        
        state(i+1, 11:13) = corrected_state(5:7);
        state(i+1, 7:10) = corrected_state(1:4);
    end 
    
    % Reaction Wheel PID Controller 
    [rw_command, prevErr, integral] = reactionWheelController(timeStep, state(i+1, :), prevErr, integral);
    w123_dot_rw(i+1, :) = rw_command;
    
    % Satellite Frame to Calculate moments and forces
    % e.g. moment due to reaction wheel acceleration
    [dF, dM] = satelliteFrame(w123_dot_rw(i+1, :)'); % state(14-16) == w123_rw
    
end


%% Extract Data for Plotting
xyz_out_rk4 = state(:, 1:3);
q_out = state(:,7:10);
pqr_out = state(:, 11:13);
phi_theta_psi_out = Quaternions2EulerAngles(q_out);
w123_rw = state(:, 14:16); % reaction wheel anglure velocity


%% Plot Results
%%% Orbit
fig = figure();
set(fig, 'color', 'white');
plot3(xyz_out_rk4(:, 1), xyz_out_rk4(:, 2), xyz_out_rk4(:, 3), 'b-', 'LineWidth', 4);
grid on 
hold on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Orbit Locus Simulation');
%%% Earth Sphere
[X, Y, Z] = sphere(100);
X = X*Planet_Radius;
Y = Y*Planet_Radius;
Z = Z*Planet_Radius;
surf(X, Y, Z, 'edgeColor', 'none'); 
colormap([0 0 0]);
axis equal

%% Magnetic Field
% fig2 = figure();
% set(fig2, 'color', 'white')
% subplot(3, 1, 1);
% plot(t_span, BB_out(:, 1), 'b-', 'LineWidth', 2); hold on;
% plot(t_span, BB_meas_out(:, 1), 'k--', 'LineWidth', 2);hold on;
% legend('Bx (IGRF)', 'Bx Sensor'); ylabel('(Tesla)');
% grid on;
% title("Magnitic Field (IGRF model) vs (Sensor Measures)");
% 
% subplot(3, 1, 2);
% plot(t_span, BB_out(:, 2), 'g-', 'LineWidth', 2); hold on;
% plot(t_span, BB_meas_out(:, 2), 'k--', 'LineWidth', 2); hold on;
% legend('By (IGRF)', 'By Sensor'); ylabel('(Tesla)'); grid on;
% 
% subplot(3, 1, 3);
% plot(t_span, BB_out(:, 3), 'r-', 'LineWidth', 2); hold on;
% plot(t_span, BB_meas_out(:, 3), 'k--', 'LineWidth', 2); hold on;
% grid on; legend('Bz (IGRF)', 'Bz Sensor');
% xlabel('Time(sec)')
% ylabel('(Tesla)');
% 
% 
% %%% Magnatic Field Magnitute
% fig3 = figure();
% set(fig3, 'color', 'white');
% plot(t_span, sqrt(sum(BB_out.^2, 2)), 'LineWidth', 2);
% xlabel('time (sec)'); ylabel('Tesla');
% title('Magnitute of Magnitic Field (T)');
% grid on;
%% Velocities
figAngles = figure();
set(figAngles, 'color', 'white');
subplot(2, 1, 1);
plot(t_span, state(:, 4)*18/5, 'r-', 'LineWidth', 2); hold on;
plot(t_span, state(:, 5)*18/5, 'b-','LineWidth', 2); hold on;
plot(t_span, state(:, 6)*18/5, 'g-','LineWidth', 2); hold on;
grid on;
ylabel('Km/hr');
legend('v_x', 'v_y', 'v_z');
title('Velicties in Inertial Earth Centerd Frame');

subplot(2, 1, 2);
velocity_mag = sqrt(sum(state(:, 4:6).^2, 2));
plot(t_span, velocity_mag*18/5, 'k-','LineWidth', 2); hold on;
grid on;
title('Velocity Mangintude In ECIF');
xlabel('time (sec)');
ylabel('Km/hr');
%% Angles
figAngles = figure();
set(figAngles, 'color', 'white');
subplot(2, 1, 1);
plot(t_span, phi_theta_psi_out(:, 1)*180/pi, 'r-', 'LineWidth', 2); hold on;
plot(t_span, phi_theta_psi_out(:, 2)*180/pi, 'b-','LineWidth', 2); hold on;
plot(t_span, phi_theta_psi_out(:, 3)*180/pi, 'g-','LineWidth', 2); hold on;
grid on;
ylabel('deg');
legend('phi', 'theta', 'psi');
title('Eulre Angles');

subplot(2, 1, 2);
plot(t_span, q_out(:, 1), 'k-','LineWidth', 2); hold on;
plot(t_span, q_out(:, 2), 'r-','LineWidth', 2); hold on;
plot(t_span, q_out(:, 3), 'b-','LineWidth', 2); hold on;
plot(t_span, q_out(:, 4), 'g-','LineWidth', 2); hold on;
grid on;
legend('q0', 'q1', 'q2', 'q3');
title('Quaternions');
xlabel('time (sec)');
%% Angle Rates
figAngleRates = figure();
subplot(3, 1, 1);
set(figAngleRates, 'color', 'white');
plot(t_span, pqr_meas_out(:, 1)*180/pi, 'b--', 'LineWidth', 2);hold on;
plot(t_span, pqr_out(:, 1)*180/pi, 'r-','LineWidth', 2); hold on;
yline(0, 'LineWidth', 2); grid on;
legend('P Gyro', 'P ES-EKF', 'Command'); title('P (deg/s)');ylabel('deg/s');

subplot(3, 1, 2);
plot(t_span, pqr_meas_out(:, 2)*180/pi, 'b--', 'LineWidth', 2);hold on;
plot(t_span, pqr_out(:, 2)*180/pi, 'r-','LineWidth', 2); hold on;
yline(0, 'LineWidth', 2); grid on;
legend('Q Gyro', 'Q ES-EKF', 'Command'); title('Q (deg/s)');ylabel('deg/s');

subplot(3, 1, 3);
plot(t_span, pqr_meas_out(:, 3)*180/pi, 'b--', 'LineWidth', 2);hold on;
plot(t_span, pqr_out(:, 3)*180/pi, 'r-','LineWidth', 2); hold on;
yline(0, 'LineWidth', 2);
legend('R Gyro', 'R ES-EKF', 'Command'); title('R (deg/s)');ylabel('deg/s');
grid on;
xlabel('time (sec)');
ylabel('deg/sec');

%% Reaction Wheel anguler velocities
figAngleRatesRW = figure();
subplot(2, 1, 1);
set(figAngleRatesRW, 'color', 'white');
plot(t_span, w123_rw*180/pi, 'LineWidth', 2);
grid on;
xlabel('Time (sec)'); ylabel('deg/s');
title('Reaction Wheel Angular Velocities');
legend('Reaction Wheel 1', 'Reaction Wheel 2', 'Reaction Wheel 3');
%%% Reaction Wheel anguler acceleration
subplot(2, 1, 2);
plot(t_span, w123_dot_rw*180/pi, 'LineWidth', 2);
grid on;
xlabel('Time (sec)'); ylabel('deg/s^2');
title('Reaction Wheel Angular Acceleration');
legend('Reaction Wheel 1', 'Reaction Wheel 2', 'Reaction Wheel 3');


%% Total Execution Time
toc;


