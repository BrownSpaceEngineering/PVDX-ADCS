function m = torque2moment3axis(tau_req, B)
    % torque2moment3axis
    % Converts desired torque into magnetic moments for 3 orthogonal magnetorquers
    %
    % Inputs:
    %   tau_req : 3x1 desired torque [N*m]
    %   B       : 3x1 magnetic field [Tesla]
    %
    % Output:
    %   m       : 3x1 magnetic moment [mx my mz] [A m^2]
    
    tau_req = tau_req(:);
    B = B(:);

    % Skew-symmetric matrix of B
    Bx = [  0     -B(3)   B(2);
           B(3)     0    -B(1);
          -B(2)   B(1)     0 ];
    
    % Solve minimum-norm magnetic dipole
    m = pinv(Bx) * tau_req;

end
