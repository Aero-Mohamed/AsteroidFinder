function R = RotationIB(phi, theta, psi)

[S, C, ~] = SCT([phi, theta, psi]);

R = [...
    C.theta*C.epsi, (S.phi*S.theta*C.epsi - C.phi*S.epsi), (C.phi*S.theta*C.epsi + S.phi*S.epsi);
    C.theta*S.epsi, (S.phi*S.theta*S.epsi + C.phi*C.epsi), (C.phi*S.theta*S.epsi - S.phi*C.epsi);
    -S.theta, S.phi*C.theta, C.phi*C.theta];

end