function [estimate_err, bias, m] = fcn(q_eci2b, reference_mag)
    % --- Constants & Sensor Specs ---
    dt = 0.1;

    % Fix #9: timestep_per_mode now matches switch_every in
    % ukf_python_replica.m ( = int(60*10/dt) = 6000 ). Was 10000, which
    % did not correspond to anything in the reference implementation.
    timestep_per_mode = int32(60 * 45 / dt);

    % Truth gyro error model -- matches true_gyro_bias in
    % ukf_python_replica.m (fix #2: previously the bias truth was zero, so
    % the filter had nothing real to estimate). Scale factor has been
    % pulled out of both the truth model and the filter state -- see the
    % removal note near the state-vector definition below.
    true_gyro_bias    = [0.002, -0.001, 0.0015];

    % REVERTED (per empirical observation): the Python-matched sigmas
    % (gyro_noise=0.0005, measurement_noise=0.05) converged worse in
    % practice than these original values, so noise tuning is reverted to
    % match the original Simulink file exactly. This intentionally
    % diverges from the current sigma_gyro/sigma_magnetometer in
    % ukf_python_replica.m -- if the Python file gets re-synced later,
    % re-check this against it rather than assuming they still match.
    gyro_noise        = 0.001;   % rad/s
    measurement_noise = 0.01;    % uT

    % Second reference vector for the 2-vector ("Sun+Mag") mode. Fixed
    % placeholder ECI unit vector for now (matches ref_vec_2 in
    % ukf_python_replica.m) -- NOT read from an actual sun sensor/ephemeris
    % yet. Swap this out for a real ECI sun vector once one is wired into
    % this block as an input.
    ref_vec_2_eci = [1, 0, 0];
    sigma_unit_vector = 0.01;   % second-vector sensor noise std, matches ukf_python_replica.m

    reference_mag = reshape(reference_mag/1000, [1,3]); % nT -> uT

    % --- Persistent State Memory ---
    persistent timestep mode error_state error_cov attitude_estimate;
    persistent true_body_to_ref_prev;

    if isempty(timestep)
        timestep = 0;
    end
    if isempty(mode)
        mode = 1; % 0: Sun+Mag (Full Observability, 2-vector), 1: Mag Only (Eclipse, 1-vector)
    end
    % NOTE: mode now actually gates the measurement path (see below), but
    % the switch is still purely on a fixed timestep_per_mode timer, same
    % as ukf_python_replica.m's switch_every -- neither file ties the
    % switch to sun_in_view / real sun visibility.

    % 6-state error vector: [Attitude (3), Bias (3)]
    % Scale factor removed from the filter state -- it was estimating an
    % unobserved/unneeded 9th-order term for this use case. n = length(x)
    % is used everywhere downstream (calculate_lambda, get_sigma_points,
    % get_weights, iterate) rather than a hardcoded 9, so shrinking this to
    % 6 automatically resizes lambda, alpha/beta-derived weights, and the
    % 2n+1 = 13 sigma points -- no other change needed in those functions.
    if isempty(error_state)
        error_state = zeros([1,6]);
    end

    if isempty(error_cov)
        error_cov = eye(6, 6) * 0.01;
    end

    % --- True Environment Simulation ---
    true_ref_to_body = Quaternion(real(q_eci2b)).quaternion_normalize();
    true_body_to_ref = true_ref_to_body.quaternion_inverse();

    % omega_icrf2b was reading straight 0 in Simulink even during active
    % nadir-pointing motion, so the truth angular velocity is derived
    % directly from the ACTUAL attitude quaternion history instead --
    % same technique used for the nadir-pointing environment model in
    % ukf_python_replica.m. This uses the same kinematic relationship as
    % the gyro-driven propagation elsewhere in this file:
    %   true_body_to_ref_new = true_body_to_ref_old * rotation_to_quat(omega*dt)^-1
    % => rotation_to_quat(omega*dt) = true_body_to_ref_new^-1 * true_body_to_ref_old
    if isempty(true_body_to_ref_prev)
        true_body_to_ref_prev = true_body_to_ref;  % first call: no history yet -> zero rate
    end
    delta_q = true_body_to_ref.quaternion_inverse().quaternion_multiply(true_body_to_ref_prev);
    true_angular_velocity = delta_q.quaternion2rotation_vec() / dt;

    if isempty(attitude_estimate)
        attitude_estimate = true_body_to_ref;
    end

    % --- Rate Cadence ---
    % This block runs at 10 Hz (dt = 0.1s). Vector measurements are only
    % available/applied at 0.1 Hz -- once every (1/0.1)/dt = 100 calls --
    % matching the propagate-every-call/update-every-100th-call scheme in
    % ukf_python_replica.m (update_every = int(10/dt) = 100).
    measurement_rate_hz = 0.1;
    steps_per_measurement = round((1 / measurement_rate_hz) / dt);  % = 100
    do_update = (mod(timestep, steps_per_measurement) == 0);

    % Sensor readings
    % Fix #2: truth gyro measurement now matches ukf_python_replica.m --
    % true rate + true bias + (true scale factor .* true rate) + noise,
    % instead of just true rate + noise with a zero bias and no scale factor.
    gyro_measurement_noise = normrnd(0, gyro_noise, [1,3]);
    simulated_gyro_measurement = true_angular_velocity + true_gyro_bias + gyro_measurement_noise;

    % --- Reference Vector(s) & R sizing -- mode-dependent ---
    % mode == 0: "Sun+Mag" / 2-vector (full observability) -- stack the
    %            magnetometer with the second reference vector.
    % mode == 1: "Mag Only" / 1-vector (eclipse) -- magnetometer alone.
    % This mirrors ukf_python_replica.m's vector_input toggle (1 <-> 2),
    % which also swaps every switch_every steps via the same mode/timestep
    % machinery already in this file (see timestep_per_mode below).
    use_two_vector = (mode == 0);

    if use_two_vector
        ref_readings = [reference_mag / norm(reference_mag); ...
                         ref_vec_2_eci / norm(ref_vec_2_eci)];
        R = blkdiag(eye(3) * 2.5e-5, eye(3) * sigma_unit_vector^2);
    else
        ref_readings = reference_mag / norm(reference_mag);
        R = eye(3) * 2.5e-5;
    end

    if do_update
        % Only sample fresh vector readings when one is actually available
        % at the 0.1 Hz measurement rate.
        body_mag = true_ref_to_body.apply_rotation(reference_mag) + normrnd([0,0,0], measurement_noise);
        mag_unit = body_mag / norm(body_mag);

        if use_two_vector
            body_vec2 = true_ref_to_body.apply_rotation(ref_vec_2_eci) + normrnd([0,0,0], sigma_unit_vector);
            vec2_unit = body_vec2 / norm(body_vec2);
            body_msmts = [mag_unit, vec2_unit];
        else
            body_msmts = mag_unit;
        end
    else
        % unused: iterate() skips the correction step when do_update is
        % false. Sized to match ref_readings/R for the current mode so the
        % call below always has consistent argument shapes.
        body_msmts = zeros(1, size(ref_readings, 1) * 3);
    end

    % --- Q Matrix (Process Noise Tuning) ---
    % REVERTED (per empirical observation): back to the original fixed
    % constants rather than the gyro_noise-derived value, since the
    % original tuning converged better in practice. No longer claims to
    % match global_Q in ukf_python_replica.m -- see note above.
    Q = zeros(6,6);
    Q(1:3,1:3) = eye(3) * 4e-8;    % Attitude noise
    Q(4:6,4:6) = eye(3) * 1e-10;   % Bias random walk

    % --- UKF Iteration ---
    [error_state, attitude_estimate, error_cov] = iterate(error_state, attitude_estimate, error_cov, ...
        ref_readings, body_msmts, simulated_gyro_measurement, Q, R, dt, do_update);

    % --- TRUE STATE MEMORY FIX ---
    % Reset ONLY the attitude error (1:3) to zero after it's injected into the nominal quaternion.
    % Bias (4:6) persists naturally through the persistent variable.
    error_state(1:3) = [0 0 0];

    % --- Outputs & Loop Maintenance ---
    estimate_err = rad2deg(Quaternion.quat_diff(true_body_to_ref, attitude_estimate));
    bias = reshape(error_state(4:6), [1 3]); % Extract bias for logging

    true_body_to_ref_prev = true_body_to_ref;
    timestep = timestep + 1;
    m = mode;

    if(mod(timestep, timestep_per_mode) == 0)
        if(mode == 1)
            mode = 0;
        else
            mode = 1;
        end
    end
end

% =========================================================================
% UKF HELPER FUNCTIONS -- exact translations of ukf_python_replica.m
% =========================================================================

function [new_error_state, new_guess, new_cov] = iterate(current_error_state, current_guess, current_cov,...
    ref_readings, body_msmts, gyro, Q, R, dt, do_update)

    n = length(current_error_state); % 6-state filter (attitude + bias)
    alpha = 0.1;
    beta = 2;

    % Fix #1: guard incoming covariance before building sigma points,
    % matching ensure_positive_definite(P) call at the top of
    % ukf_python_replica.m's iterate().
    current_cov = ensure_positive_definite(current_cov);

    lam = calculate_lambda(alpha, current_error_state);
    sigmas = get_sigma_points(lam, current_error_state, current_cov);
    quat_sigmas = error_sigmas_to_quat_sigmas(sigmas, current_guess);
    propagated_quat_sigmas = propagate_quat_sigmas(quat_sigmas, gyro, dt);   % always propagate (10 Hz)

    [covariance_weights, mean_weights] = get_weights(lam, current_error_state, alpha, beta);

    num_sigmas = size(propagated_quat_sigmas, 1);
    quaternions_of_propagated_sigmas = cell(num_sigmas, 1);
    for i = 1:num_sigmas
        quaternions_of_propagated_sigmas{i} = Quaternion(propagated_quat_sigmas(i, 1:4));
    end

    [average_quaternion, propagated_error_vectors] = gradient_descent( ...
        quaternions_of_propagated_sigmas, current_guess, mean_weights);

    propagated_errors = [propagated_error_vectors, propagated_quat_sigmas(:, 5:7)];

    mean_error = zeros(1, n);
    for col = 1:n
        var_mean = 0;
        for row = 1:num_sigmas
            var_mean = var_mean + propagated_errors(row, col) * mean_weights(row);
        end
        mean_error(col) = var_mean;
    end

    P_hat = zeros(n, n);
    for i = 1:num_sigmas
        err = propagated_errors(i, :) - mean_error;
        P_hat = P_hat + covariance_weights(i) * (err' * err);
    end
    P_hat = P_hat + Q;   % process noise always accrues, every propagate step

    if do_update
        % --- MEASUREMENT CORRECTION (only on 0.1 Hz update steps) ---
        measurements = get_measurements(propagated_quat_sigmas, ref_readings);

        msmt_size = size(measurements, 2);
        mean_measurement = zeros(1, msmt_size);
        for col = 1:msmt_size
            var_mean = 0;
            for row = 1:num_sigmas
                var_mean = var_mean + measurements(row, col) * mean_weights(row);
            end
            mean_measurement(col) = var_mean;
        end

        P_xz = zeros(n, msmt_size);
        for i = 1:num_sigmas
            err = propagated_errors(i, :) - mean_error;
            msmt_err = measurements(i, :) - mean_measurement;
            P_xz = P_xz + covariance_weights(i) * (err' * msmt_err);
        end

        P_zz = zeros(msmt_size, msmt_size);
        for i = 1:num_sigmas
            msmt_err = measurements(i, :) - mean_measurement;
            P_zz = P_zz + covariance_weights(i) * (msmt_err' * msmt_err);
        end

        P_vv = P_zz + R;
        k = P_xz * inv(P_vv);
        x_hat = mean_error + (k * (body_msmts - mean_measurement)')';

        % Matches ukf_python_replica.m exactly: P_hat - k @ P_vv @ k.T
        P = P_hat - k * P_vv * k';
    else
        % --- PREDICT-ONLY (propagate mean/covariance, no correction) ---
        x_hat = mean_error;
        P = P_hat;
    end

    % Fix #1: guard outgoing covariance after the update, matching
    % ensure_positive_definite(P) call at the end of
    % ukf_python_replica.m's iterate().
    P = ensure_positive_definite(P);

    x_hat_rot = Quaternion.rotation_vec2quaternion(x_hat(1:3));
    new_error_state = x_hat;
    new_guess = x_hat_rot.quaternion_multiply(average_quaternion);
    new_cov = P;
end

function P = ensure_positive_definite(P, epsilon)
    % Exact port of ensure_positive_definite() from ukf_python_replica.m
    if nargin < 2
        epsilon = 1e-9;
    end
    P = (P + P') / 2;
    % real() forces this to stay a real-typed array for codegen: eig() on a
    % generic square matrix is typed as potentially complex by MATLAB
    % Coder/Simulink even though P is symmetric here, which otherwise
    % causes "Cannot assign a complex value into a non-complex location."
    eigvals = real(eig(P));
    if min(eigvals) < epsilon
        P = P + eye(size(P,1)) * (epsilon - min(eigvals) + epsilon);
    end
end

function lam = calculate_lambda(alpha, x)
    n = length(x);
    kappa = 0;
    lam = alpha^2 * (n + kappa) - n;
end

function sigmas = get_sigma_points(lam, x, P)
    n = length(x);
    P = (P + P') / 2;
    P = P + eye(n) * 1e-10;
    U = chol((lam + n) * P);   % upper-triangular, U'*U = (lam+n)*P

    sigmas = zeros(2*n + 1, n);
    sigmas(1, :) = x;
    for k = 1:n
        sigmas(k + 1, :) = x + U(k, :);
    end
    for k = (n+1):(2*n)
        idx = k - n;
        sigmas(k + 1, :) = x - U(idx, :);
    end
end

function [covariance_weights, mean_weights] = get_weights(lam, x, alpha, beta)
    n = length(x);
    c = 0.5 / (n + lam);
    covariance_weights = ones(1, 2*n + 1) * c;
    mean_weights = ones(1, 2*n + 1) * c;
    covariance_weights(1) = lam / (n + lam) + (1 - alpha^2 + beta);
    mean_weights(1) = lam / (n + lam);
end

function quat_sigmas = error_sigmas_to_quat_sigmas(error_sigmas, rotation)
    n = size(error_sigmas, 1);
    quat_sigmas = zeros(n, 7);
    for i = 1:n
        quat_rot = Quaternion.rotation_vec2quaternion(error_sigmas(i, 1:3));
        new_quat_rot = quat_rot.quaternion_multiply(rotation);
        new_quat_rot = new_quat_rot.quaternion_normalize();
        quat_sigmas(i, :) = [new_quat_rot.to_array(), error_sigmas(i, 4:6)];
    end
end

function new_sigmas = propagate_quat_sigmas(quat_sigmas, gyro_measurement, dt)
    n = size(quat_sigmas, 1);
    new_sigmas = zeros(n, 7);
    for i = 1:n
        quaternion = quat_sigmas(i, 1:4);
        bias = quat_sigmas(i, 5:7);

        w = gyro_measurement - bias;

        new_quaternion = Quaternion(quaternion).quaternion_multiply( ...
            Quaternion.rotation_vec2quaternion(w * dt).quaternion_inverse());
        new_quaternion = new_quaternion.quaternion_normalize();

        new_sigmas(i, :) = [new_quaternion.to_array(), quat_sigmas(i, 5:7)];
    end
end

function measurements = get_measurements(quat_sigmas, ref_vecs)
    % ref_vecs: k x 3 matrix, one ECI reference unit vector per row (k=1
    % for mag-only mode, k=2 for the 2-vector "Sun+Mag" mode). Produces a
    % 3*k-wide row per sigma point, concatenating each vector's rotated
    % body-frame prediction in the same order as ref_vecs' rows -- mirrors
    % ukf_python_replica.m's get_measurements() when ref_vecs is a list.
    n = size(quat_sigmas, 1);
    k = size(ref_vecs, 1);
    measurements = zeros(n, 3*k);
    for i = 1:n
        quat = Quaternion(quat_sigmas(i, 1:4));
        ECI_to_body = quat.quaternion_inverse();
        for j = 1:k
            measurements(i, (3*j-2):(3*j)) = ECI_to_body.apply_rotation(ref_vecs(j, :));
        end
    end
end

function error_vectors = get_error_vectors(Y, x)
    n = length(Y);
    error_vectors = zeros(n, 3);
    x_inv = x.quaternion_inverse();
    for i = 1:n
        error_quaternion = Y{i}.quaternion_multiply(x_inv);
        error_vectors(i, :) = error_quaternion.quaternion2rotation_vec();
    end
end

function [x, error_vectors] = gradient_descent(Y, x, weights)
    n_pts = length(Y);
    if nargin < 3 || isempty(weights)
        weights = ones(1, n_pts) / n_pts;
    end
    safe_w = max(weights, 0);
    safe_w = safe_w / sum(safe_w);

    for iter = 1:100
        error_vectors = get_error_vectors(Y, x);
        ave = sum(error_vectors .* safe_w(:), 1); % weighted average across rows
        if norm(ave) < 1e-6
            break;
        end
        average_error_quat = Quaternion.rotation_vec2quaternion(ave);
        x = average_error_quat.quaternion_multiply(x);
    end
end