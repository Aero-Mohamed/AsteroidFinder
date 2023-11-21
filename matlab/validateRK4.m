
%% I used this code to validate RK4 Results with Ode45

[tout, state_out] = ode45(@Satellite, t_span, initail_state);
xyz_out_ode45 = state_out(:, 1:3);

%% Plot
%%% Orbit
figValidateOrbit = figure();
set(figValidateOrbit, 'color', 'white');
title('Orbit Validation ode45 vs rk4');
grid on 
plot3(xyz_out_ode45(:, 1), xyz_out_ode45(:, 2), xyz_out_ode45(:, 3), 'b-', 'LineWidth', 4);
hold on;
plot3(xyz_out(:, 1), xyz_out(:, 2), xyz_out(:, 3), 'r-', 'LineWidth', 4);
xlabel('X'); ylabel('Y'); zlabel('Z');
%%% Earth Sphere
[X, Y, Z] = sphere(100);
X = X*R;
Y = Y*R;
Z = Z*R;
surf(X, Y, Z, 'edgeColor', 'none');
legend({'ODE45', 'RK4', 'Earth'});
axis equal
