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
planet              %% Radius, Mass, Gravity Const., Mu
satellite_inertia;  %% Mass, Inertia
filter_data;        %% gyroCovariance, Q, R, H_jac, P_0 for kalman filter

%% Initial Conditions 

altitute = 600*1000; %% Meters
xyz_0 = [Planet_Radius + altitute, 0, 0];

% Assume Circular Orbit for simplicity.
semi_major = norm(xyz_0);
v_circ = sqrt(Mu/semi_major);  % Planet Orbital Velocity
inclination = 56*pi/180;       % Orbit Inclination

xyz_dot_0 = [0, v_circ*cos(inclination), v_circ*sin(inclination)];
phi_theta_psi_0 = [0, 0, 0*pi/180];
q0123_0 = EulerAngles2Quaternions(phi_theta_psi_0)';
pqr_0 = [0.0*pi/180, 0*pi/180, -0*pi/180];

initail_state = [xyz_0, xyz_dot_0, q0123_0, pqr_0]';

%% Time Frame
period = 2*pi/sqrt(Mu)*semi_major^(3/2);  % one full orbit time
timeStep = 1;
t_span = 0:timeStep:period;
magFieldTimeStep = 100; lastMagFieldTime = 0;
sensorTimeStep = 1; lastSensorTime = 0;

%% Simulation
state = NaN(length(t_span), length(initail_state));
state(1, :) = initail_state';
% Store Magnitic Field readings from IGRF model & Sensor
BB_out = NaN(length(t_span), 3);
BB_meas_out = NaN(length(t_span), 3);
% Store pqr reading from sensor
pqr_meas_out = NaN(length(t_span), 3);
% Kalman Filter Covariance Matrix
P_cov = NaN(9, 9, length(t_span)); % multiple diminsions 
P_cov(:, :, 1) = P_0;

for i = 1:length(t_span)-1

    currentTime = t_span(i);
    % For External Forces and Moments
    forces = [0, 0, 0]';
    moments = [0, 0, 0]';
    % 
    state(i+1, :) = SatelliteRK4(timeStep, state(i, :), forces, moments);
    BB_out(i, :) = BB; % Get Magnetic Field from IGRF model
    
    %%% Kalman Filter F jacobian matrix
    currentState = state(i+1, :);
    prevState = state(i, :);
    
    F = kalmanSystemJacobian(currentState, prevState, I, invI);
    P_cov(:, :, i+1) = F*P_cov(:, :, i)*F' + Q;
    
    % Sensor Measurements | Assume measurements are available every second.
    if (currentTime >= lastSensorTime)
        pqr = currentState(11:13);
        [BB_meas, pqr_meas] = Sensors(BB, pqr);
        BB_meas_out(i, :) = BB_meas;
        pqr_meas_out(i, :) = pqr_meas;
        lastSensorTime = currentTime + sensorTimeStep;
        
        % let the state be [ q234 pqr TxTyTz] 
        % T is disturbance Torque
        state_estimate = [state(i+1, 8:10), state(i+1, 11:13), [0 0 0]];
        [corrected_state, p_cov_corrected] = ES_EKF(varGyro, P_cov(:, :, i+1), state_estimate, pqr_meas_out(i, :));
        P_cov(:, :, i+1) = p_cov_corrected;
        
        state(i+1, 11:13) = corrected_state(5:7);
        state(i+1, 7:10) = corrected_state(1:4);
        
    end 
    
end
%% Extract Data for Plotting
xyz_out_rk4 = state(:, 1:3);
q_out = state(:,7:10);
pqr_out = state(:, 11:13);
phi_theta_psi_out = Quaternions2EulerAngles(q_out);


%% Plot
%%% Orbit
fig = figure();
set(fig, 'color', 'white');
plot3(xyz_out_rk4(:, 1), xyz_out_rk4(:, 2), xyz_out_rk4(:, 3), 'b-', 'LineWidth', 4);
grid on 
hold on;
xlabel('X'); ylabel('Y'); zlabel('Z');
%%% Earth Sphere
[X, Y, Z] = sphere(100);
X = X*Planet_Radius;
Y = Y*Planet_Radius;
Z = Z*Planet_Radius;
surf(X, Y, Z, 'edgeColor', 'none');
axis equal

%% Magnetic Field
fig2 = figure();
set(fig2, 'color', 'white')
subplot(3, 1, 1);
plot(t_span, BB_out(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(t_span, BB_meas_out(:, 1), 'k--', 'LineWidth', 2);hold on;
legend('Bx (IGRF)', 'Bx Sensor'); ylabel('(Tesla)');

subplot(3, 1, 2);
plot(t_span, BB_out(:, 2), 'g-', 'LineWidth', 2); hold on;
plot(t_span, BB_meas_out(:, 2), 'k--', 'LineWidth', 2); hold on;
legend('By (IGRF)', 'By Sensor'); ylabel('(Tesla)');

subplot(3, 1, 3);
plot(t_span, BB_out(:, 3), 'r-', 'LineWidth', 2); hold on;
plot(t_span, BB_meas_out(:, 3), 'k--', 'LineWidth', 2); hold on;
grid on; legend('Bz (IGRF)', 'Bz Sensor');
xlabel('Time(sec)')
ylabel('(Tesla)');
title("Magnitic Field (IGRF model) vs (Sensor Measures)");

%%% Magnatic Field Magnitute
fig3 = figure();
set(fig3, 'color', 'white');
plot(t_span, sqrt(sum(BB_out.^2, 2)), 'LineWidth', 2);
xlabel('time (sec)'); ylabel('Tesla');
title('Magnitute of Magnitic Field (nT)');
grid on;

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
legend('p sensor', 'p model'); title('P (deg/s)');ylabel('deg/s');

subplot(3, 1, 2);
plot(t_span, pqr_meas_out(:, 2)*180/pi, 'b--', 'LineWidth', 2);hold on;
plot(t_span, pqr_out(:, 2)*180/pi, 'r-','LineWidth', 2); hold on;
legend('q sensor', 'q model'); title('Q (deg/s)');ylabel('deg/s');

subplot(3, 1, 3);
plot(t_span, pqr_meas_out(:, 3)*180/pi, 'b--', 'LineWidth', 2);hold on;
plot(t_span, pqr_out(:, 3)*180/pi, 'r-','LineWidth', 2); hold on;
legend('r sensor', 'r model'); title('R (deg/s)');ylabel('deg/s');
grid on;
xlabel('time (sec)');
ylabel('deg/sec');




%% Total Execution Time
toc;


