function F = kalmanSystemJacobian(currentState, prevState, I, invI)
    %%% values needed in F Jacobian Computation
    pqr = currentState(11:13);
    p = pqr(1); q = pqr(2); r = pqr(3);
    I11 = I(1, 1); I22 = I(2, 2); I33 = I(3, 3);
    delta_pqr = currentState(11:13) - prevState(11:13);
    dp = delta_pqr(1); dq = delta_pqr(2); dr = delta_pqr(3);
    
    delta_quat = currentState(7:9) - prevState(7:9);
    dq1 = delta_quat(1); dq2 = delta_quat(2); dq3 = delta_quat(3);

    
    % the following Jacobian matrix is as suggested by the given 
    % Diploma Thesis (Development of ADC for AsteroidFinder Satellite)
    F = zeros(9, 9);    
    F(1:3, 1:3) = [          0,   dr/2 + r, - dq/2 - q;
                    - dr/2 - r,          0,   dp/2 + p;
                      dq/2 + q, - dp/2 - p,          0];
    F(1:3, 4:6) = [    0, -dq3,  dq2;
                     dq3,    0, -dq1;
                    -dq2,  dq1,    0];
    
    F(4:6, 4:6) = invI*[                              0, I22*dr - I33*dr + I22*r - I33*r, I22*dq - I33*dq + I22*q - I33*q;
                        I33*dr - I11*dr - I11*r + I33*r,                               0, I33*dp - I11*dp - I11*p + I33*p;
                        I11*dq - I22*dq + I11*q - I22*q, I11*dp - I22*dp + I11*p - I22*p,                               0];
    F(4:6, 7:9) = invI;
    F(7:9, 7:9) = eye(3, 3);
end