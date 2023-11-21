% RBD Solver is function that implements Runge-Kutta-4
% Algorithm in order to integrate 6 DOF set of equations
% with a give Inital Conditions ICs, for single dt time step.
function state = SatelliteRK4(dt, ICs, forces, moments)
    K = zeros(length(ICs), 4);
    
    K(:, 1) = dt*Satellite(0, ICs', forces, moments);
    K(:, 2) = dt*Satellite(0, ICs'+0.5*K(:, 1), forces, moments );
    K(:, 3) = dt*Satellite(0, ICs'+0.5*K(:, 2), forces, moments );
    K(:, 4) = dt*Satellite(0, ICs'+K(:, 3), forces, moments );
    
    state = ICs' + (...
        K(:, 1)+...
        2*K(:, 2)+...
        2*K(:, 3)+...
        K(:, 4))/6;
    state = state';
end