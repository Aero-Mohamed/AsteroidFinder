
global H_jac

varGyro = 1e-2;
H_jac = [zeros(3, 3), eye(3, 3), zeros(3, 3)];


Q = 1e-5*eye(9, 9);
%Q(4:6, 4:6) = 1e-5*eye(3, 3);

P_0 = eye(9, 9)*10;


