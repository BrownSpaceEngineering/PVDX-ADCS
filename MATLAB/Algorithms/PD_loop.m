% INPUTS:
% omega - (3 x 1), deg/s: w_eci2b in body coordinates, that is angular velocity wr/ to eci
% written in body coordinates
% q_error - (4 x 1): error quaternion from pointing_error, scalar first

% RETURN:
% tau - torque in N-m, body coordinates


function tau = PD_loop(omega, q_error)

    arguments
        omega (3,1) double
        q_error (4, 1) double
    end
    %find shortest possible rotation
    if (q_error(1) < 0)
        q_error = -q_error;
    end

   
    angle = 2 * acos(q_error(1));

    if angle > 1e-12
        axis = q_error(2:4) / sin(angle/2);
    else
        axis = [0; 0; 0];
    end
    r_e = [axis; angle];
    % Convert 3x1 error omega vector back to axis-angle
    omega_rad = omega * (pi/180); % omega is in deg/s. converting to rad/s
    omega_mag = norm(omega_rad);
   
    if omega_mag > 1e-12
        omega_axis = omega_rad / omega_mag;
    else
        omega_axis = [0; 0; 0];
    end
    
    r_omega = [omega_axis; omega_mag];

    % PD controller
    Kp0 = 0.034;
    Kd0 = 0.42;
    max_tau = 5;
    min_tau = -5;
    I_body = [3.0054115e-02, -5.1674000e-05, 2.5075000e-05 ;-5.1674000e-05, 1.1430611e-02, -3.6481690e-03 ;2.5075000e-05, -3.6481690e-03, 2.3052295e-02];
    
    Kp = [Kp0, Kp0, Kp0];
    Kd = [Kd0, Kd0, Kd0];
    tau = zeros(3,1);

    for i = 1:3
        Pi = Kp(i) * r_e(4) * r_e(i);
        Di = Kd(i) * r_omega(4) * r_omega(i);
        tau(i) = Pi - Di;
    end
    tau = I_body * tau;
    tau = max(min(tau, max_tau), min_tau);
end