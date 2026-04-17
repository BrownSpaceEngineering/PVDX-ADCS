
function simulate_1vec()

    %setting constants
    dt = 0.5;
    true_body_to_ref = Quaternion([0.2812, -0.6998, 0.5497, -0.3592]).quaternion_normalize();%the true body --> reference rotation
    current_guess = [0.2812, -0.6998, 0.5497, -0.3592] + normrnd(0.05, 0.1, [1,4]);%the estimated body --> reference rotation
    current_guess = Quaternion(current_guess / norm(current_guess));

    true_omega = deg2rad([2, 0.5, -1]);%true angular velocity in the reference frame

    %in Simulink, you need to take the omega_icrf2b and do
    %omega_b = true_ref_to_body.apply_rotation(-omega_icrf2b);
    %true_omega = reshape(deg2rad(omega_b), [1,3]);
    
    gyro_bias = [0.001, 0.001, 0.001];
    gyro_noise = 0.001;
    measurement_noise = 0.03;

    state = zeros([1,6]);%starting state for the filter
    cov = eye(6) * 0.1;%starting covariance for the filter
    Q = eye(6) * 0.01;%process noise matrix
    R = diag([0.03, 0.03, 0.03, 1, 1, 1]);%measurement noise matrix for 1 vector

    reference = [40, -40, 30];
    last_reference = reference;
    reference_to_body = true_body_to_ref.quaternion_inverse();
    last_body = reference_to_body.apply_rotation(last_reference);
    err = 0;
    total_abs_err = [0 ,0, 0];
    for i = 1:1000 %iterate
        simulated_gyro_measurement = true_omega + normrnd(gyro_bias, gyro_noise);

        delta_q = Quaternion.rotation_vec2quaternion(true_omega * dt);
        true_body_to_ref = true_body_to_ref.quaternion_multiply(delta_q);
        true_body_to_ref = true_body_to_ref.quaternion_normalize();

        reference_to_body = true_body_to_ref.quaternion_inverse();%true reference to body rotation

        %get one dummy measurement in the body frame
        
        body = reference_to_body.apply_rotation(reference) + normrnd(0, measurement_noise, [1,3]);
        
        %Taking the time derivative. Works here but not in Simulink.(in
        %Simulink, something completely different works, which is 
        % 
        % finite_difference_reference = (reference - last_reference) / dt + cross(simulated_gyro_measurement, reference);
        % finite_difference_body = (body - last_body) / dt;

        %currently debugging reasons with profs
        finite_difference_reference = (reference - last_reference) / dt;
        finite_difference_body = (body - last_body) / dt + cross(simulated_gyro_measurement, body);

        ref_readings = [reference / norm(reference), finite_difference_reference];
        body_msmts = [body / norm(body), finite_difference_body];

        [state, current_guess, cov] = iterate(state, current_guess, cov, ref_readings, body_msmts,...
            simulated_gyro_measurement, Q, R, dt);
        state(1:6) = zeros([1,6]);%for 1 vector --> no bias estimation
        %state(1:3) = zeros([1,3]);
        last_body = body;
        last_reference = reference;
        reference = reference + normrnd(0, 0.001, [1, 3]);
        if mod(i, 50) == 0
            disp("Iteration");
            %disp(current_guess)
            %disp(true_rotation)
            disp(rad2deg(Quaternion.quat_diff(true_body_to_ref, current_guess)))
            
        end
        err = err + norm(rad2deg(Quaternion.quat_diff(true_body_to_ref, current_guess)));
        total_abs_err = total_abs_err + abs(rad2deg(Quaternion.quat_diff(true_body_to_ref, current_guess)));
        
        
    end
    disp(err / 1000);
    disp(total_abs_err / 1000);
    
end


function [new_error_state, new_guess, new_cov] = iterate(current_error_state, current_guess, current_cov,...
    ref_readings, body_msmts, gyro, Q, R, dt)
    n = 6;%size of the state
    current_cov = ensure_psd(current_cov);
    msmt_size = length(ref_readings);
    alpha = 0.01;
    beta = 2;
    lambda = calculate_lambda(n, alpha);
    error_sigmas = get_sigma_points(n, lambda, current_error_state, current_cov + Q);
    num_sigmas = size(error_sigmas, 1);
    quat_sigmas = error_sigmas_to_quat_sigmas(error_sigmas, current_guess);
    propagated_sigmas = propagate_sigmas(quat_sigmas, gyro, dt);

    sigma_measurements = get_sigma_measurements(propagated_sigmas, ref_readings);
    [cov_weights, mean_weights] = get_weights(lambda, n, alpha, beta);

    quaternions_of_sigmas = propagated_sigmas(:, 1:4);
    [average_quaternion, propagated_error_vectors] = gradient_descent(quaternions_of_sigmas, current_guess.to_array());
    propagated_errors = zeros([num_sigmas, 6]);
    propagated_errors(:, 1:3) = propagated_error_vectors(:, 1:3);
    propagated_errors(:, 4:6) = propagated_sigmas(:, 5:7);

    
    mean_error = zeros([1, n]);
    for col = 1:n
        var_mean = 0;
        for row = 1:num_sigmas
            var_mean = var_mean + mean_weights(row) * propagated_errors(row, col);
        end
        mean_error(col) = var_mean;
    end
    mean_msmt = zeros([1, msmt_size]);
    for col = 1:msmt_size
        var_mean = 0;
        for row = 1:num_sigmas
            var_mean = var_mean + mean_weights(row) * sigma_measurements(row, col);
        end
        
        mean_msmt(col) = var_mean;
    end
    P_hat = zeros([n, n]);

    for i = 1:num_sigmas
        err = propagated_errors(i, :) - mean_error;
        P_hat = P_hat + cov_weights(i) * (err' * err);
    end
    P_xz = zeros([n, msmt_size]);
    for i = 1:num_sigmas
        err = propagated_errors(i, :) - mean_error;
        msmt_err = sigma_measurements(i, :) - mean_msmt;
        P_xz = P_xz + cov_weights(i) * (err' * msmt_err);
    end

    P_zz = zeros([msmt_size, msmt_size]);
    for i = 1:num_sigmas
        msmt_err = sigma_measurements(i, :) - mean_msmt;
        P_zz = P_zz + cov_weights(i) * (msmt_err' * msmt_err);
    end
    P_vv = P_zz + R;
    k = P_xz * inv(P_vv);
    x_hat = mean_error + (k*(body_msmts - mean_msmt)')';
    

    P = P_hat - k * P_vv * k';
    P = ensure_psd(P);
    
    x_hat_rot = Quaternion.rotation_vec2quaternion(x_hat(1:3));
    new_error_state = x_hat;
    new_guess = average_quaternion.quaternion_multiply(x_hat_rot);
    new_cov = P;
end

function mat = ensure_psd(mat)
    eps = 0.0000001;
    mat = (mat + mat')/2;
    eigvals = real(eig(mat));
    if(min(eigvals) < eps)
        mat = mat + eye(size(mat, 1)) * (eps - min(eigvals) + eps);
    end
end
function lam = calculate_lambda(n, alpha)
    kappa = 3-n;
    lam = alpha^2 * (n + kappa) - n;
end

function sigma_points = get_sigma_points(n, lam, state, cov)
    cov = ensure_psd(cov * (lam + n));
    U = chol(cov)';
    sigmas = zeros([2 * n + 1, n]);
    sigmas(1, :) = state;
    for k = 1:n
        idx = k;
        
        sigmas(k + 1, :) = state + U(idx, :);
    end
    for k = n+1 : 2*n
        idx = k - n;
        sigmas(k+1, :) = state - U(idx, :);
    end
    sigma_points = sigmas;
end

function [cov_weights, mean_weights] = get_weights(lambda, n, alpha, beta)
    c = 0.5 / (n + lambda);
    cov_weights = ones([1, 2 * n + 1]) * c;
    mean_weights = ones([1, 2 * n + 1]) * c;
    cov_weights(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);
    mean_weights(1) = lambda / (n + lambda);
end

function quat_sigmas = error_sigmas_to_quat_sigmas(error_sigmas, current_guess)
    n = size(error_sigmas, 1);
    quat_sigmas = zeros([n, 7]);
    for i = 1:n
        error_rot = Quaternion.rotation_vec2quaternion(error_sigmas(i, 1:3));
        new_rotation = current_guess.quaternion_multiply(error_rot);
        new_rotation = new_rotation.quaternion_normalize();
        quat_sigmas(i, 1:4) = new_rotation.to_array();
        quat_sigmas(i, 5:7) = error_sigmas(i, 4:6);
    end
end

function propagated_sigmas = propagate_sigmas(quat_sigmas, gyro_measurement, dt)
    n = size(quat_sigmas, 1);
    propagated_sigmas = zeros([n, 7]);
    for i = 1:n
        rotation = quat_sigmas(i, 1:4);
        bias = quat_sigmas(i, 5:7);
        omega = gyro_measurement - bias;%active rotation from body old -> body new
        delta_q = Quaternion.rotation_vec2quaternion(omega * dt);
        %active rotation from body old -> body new
        %which is the same as a passive rotation from body new -> body old
        propagated_rotation = Quaternion(rotation).quaternion_multiply(delta_q);%body new -> body old -> reference
        propagated_rotation = propagated_rotation.quaternion_normalize();
        propagated_sigmas(i, 1:4) = propagated_rotation.to_array();

        propagated_sigmas(i, 5:7) = quat_sigmas(i, 5:7);
    end
end

function msmts = get_sigma_measurements(quat_sigmas, ref_readings)
    n = size(quat_sigmas, 1);
    msmts = zeros([n, length(ref_readings)]);
    for i = 1:n
        sigma_rotation = Quaternion(quat_sigmas(i, 1:4));
        sigma_ref_to_body = sigma_rotation.quaternion_inverse();
        expected_reading_one = sigma_ref_to_body.apply_rotation(ref_readings(1:3));
        msmts(i, 1:3) = expected_reading_one;
        if(length(ref_readings) == 6)
            expected_reading_two = sigma_ref_to_body.apply_rotation(ref_readings(4:6));
            msmts(i, 4:6) = expected_reading_two;
        end
        
        
    end
end
function error_vecs = get_error_vectors(other_quats, x)
    n = length(other_quats);
    x = Quaternion(x);
    x_inv = x.quaternion_inverse();
    error_vecs = zeros([n, 3]);
    for i = 1:n
        other_quat = Quaternion(other_quats(i, :));
        error_quaternion = x_inv.quaternion_multiply(other_quat);
        error_vecs(i, :) = error_quaternion.quaternion2rotation_vec();
    end
end
function [avg_quaternion, errors] = gradient_descent(propagated_quats, current_guess)
    for i = 1:50
        error_vecs = get_error_vectors(propagated_quats, current_guess);
        average = mean(error_vecs, 1);
        if(norm(average) < 0.001)
            break;
        end
        average_error_quat = Quaternion.rotation_vec2quaternion(average);
        current_guess = Quaternion(current_guess).quaternion_multiply(average_error_quat);
        current_guess = current_guess.to_array();
    end
   
    errors = error_vecs;
    avg_quaternion = Quaternion(current_guess);
end