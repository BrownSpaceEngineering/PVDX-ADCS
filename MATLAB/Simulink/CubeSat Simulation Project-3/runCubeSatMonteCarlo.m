function results = runCubeSatMonteCarlo(opts)
%runCubeSatMonteCarlo Monte Carlo settling study for the CubeSat PD controller.
%
%   results = runCubeSatMonteCarlo() runs 50 trials with random initial
%   inertial attitude and random initial body rates, and returns a table with
%   one row per trial containing the sampled initial conditions and the
%   measured settling time of the attitude error.
%
%   Initial conditions are injected with Simulink.SimulationInput, so the
%   model on disk is never modified.  The two ICs that matter are dialog
%   parameters of
%       asbCubeSat/Vehicle Model/Vehicle Dynamics/Spacecraft Dynamics
%   namely 'attitude'     (4x1 quaternion, ICRF->body, scalar first) and
%          'attitudeRate' (3x1 body rates, deg/s -- block angleUnits is Degrees).
%   Note the model currently holds a hard-coded literal in 'attitude' rather
%   than initCond.euler, and takes 'attitudeRate' from initCond.pqr in
%   asbCubeSatModelData.sldd.  Overriding the block parameter bypasses both.
%
%   Name-value arguments:
%     NumTrials      Number of trials                            (default 50)
%     StopTime       Simulation stop time, s                     (default 600)
%     AttitudeMode   "uniform" for uniformly random attitudes,
%                    "axisangle" for a random axis with tilt angle
%                    uniform in [MinTiltDeg MaxTiltDeg]          (default "uniform")
%                    Note the tilt is measured from ICRF identity, not from the
%                    commanded attitude -- the pointing error the controller
%                    actually sees depends on the selected pointing mode, and is
%                    reported per trial as initErr_deg.
%     MinTiltDeg     Lower tilt bound, axisangle mode only       (default 0)
%     MaxTiltDeg     Upper tilt bound, axisangle mode only       (default 180)
%     RateRangeDeg   [min max] magnitude of initial body rate, deg/s
%                                                                (default [0 5])
%     SettleTolDeg   Attitude error must stay below this, deg    (default 2)
%     SettleRateTol  Body rate magnitude must stay below this, deg/s.
%                    Inf disables the rate criterion.  Note the Earth-pointing
%                    modes hold a nonzero orbit rate (~0.07 deg/s), so do not
%                    set this below about 0.2.                   (default Inf)
%     SteadyWindow   Fraction of the run, measured back from StopTime, used
%                    for the steady-state statistics                (default 0.2)
%     Seed           RNG seed for repeatability                  (default 0)
%     Plot           Draw summary figures                        (default true)
%
%   Steady-state columns are computed over that trailing window rather than
%   from the single last sample (which finalErr_deg still reports).  Read
%   ssErrP2P_deg and ssErrSlope_degps together before trusting ssErrMean_deg:
%   a large p2p means the error is rippling or limit-cycling, and a nonzero
%   slope means it has not converged and the "steady-state" value is just a
%   reading at one point on a trend.  This orbit's period is ~5560 s, so short
%   runs sample only a few percent of an orbit and a nonzero slope there is
%   expected.
%
%   Example:
%     r = runCubeSatMonteCarlo(NumTrials=100, StopTime=900, SettleTolDeg=1.5);
%     summary(r)

    arguments
        opts.NumTrials     (1,1) double {mustBePositive, mustBeInteger} = 50
        opts.StopTime      (1,1) double {mustBePositive} = 600
        opts.AttitudeMode  (1,1) string {mustBeMember(opts.AttitudeMode, ["uniform","axisangle"])} = "uniform"
        opts.MinTiltDeg    (1,1) double = 0
        opts.MaxTiltDeg    (1,1) double = 180
        opts.RateRangeDeg  (1,2) double = [0 5]
        opts.SettleTolDeg  (1,1) double {mustBePositive} = 2
        opts.SettleRateTol (1,1) double {mustBePositive} = Inf
        opts.SteadyWindow  (1,1) double {mustBePositive} = 0.2
        opts.Seed          (1,1) double = 0
        opts.Plot          (1,1) logical = true
    end

    mdl    = 'asbCubeSat';
    dynBlk = [mdl '/Vehicle Model/Vehicle Dynamics/Spacecraft Dynamics'];

    load_system(mdl);
    cleanup = prepareModel(mdl); %#ok<NASGU> % restores model state on exit

    % This harness must never alter the model's timing. CubeSatTimeStep sets the
    % sensor Unit Delay, the q_est/w_est delays, and every rate that inherits
    % from them, so a silent change here would invalidate the whole batch.
    % Captured now, verified untouched at the end.
    dict          = Simulink.data.connect('asbCubeSatModelData.sldd');
    timeStepAtRun = dict.get('CubeSatTimeStep');
    fprintf('CubeSatTimeStep = %g s (%.4g Hz) -- harness will not modify this\n', ...
        timeStepAtRun, 1/timeStepAtRun);

    %% Sample initial conditions
    rng(opts.Seed);
    n  = opts.NumTrials;
    q0 = sampleAttitude(n, opts);
    w0 = sampleRate(n, opts.RateRangeDeg);

    %% Build the simulation batch
    in(1:n) = Simulink.SimulationInput(mdl);
    for k = 1:n
        in(k) = setBlockParameter(in(k), dynBlk, 'attitude',     mat2str(q0(k,:), 12));
        in(k) = setBlockParameter(in(k), dynBlk, 'attitudeRate', mat2str(w0(k,:), 12));
        in(k) = setModelParameter(in(k), 'StopTime', num2str(opts.StopTime));
        in(k) = in(k).setUserString(sprintf('trial%03d', k));
    end

    %% Run
    if license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        out = parsim(in, 'ShowProgress', 'on', 'StopOnError', 'off');
    else
        out = sim(in, 'ShowProgress', 'on', 'StopOnError', 'off');
    end

    %% Post-process
    results = table('Size', [n 13], ...
        'VariableTypes', [{'double','cell','cell'}, repmat({'double'}, 1, 9), {'logical'}], ...
        'VariableNames', {'Trial','q0','w0_degps','w0Mag_degps','initErr_deg', ...
                          'peakErr_deg','finalErr_deg','ssErrMean_deg','ssErrRMS_deg', ...
                          'ssErrP2P_deg','ssErrSlope_degps','settleTime_s','settled'});
    ssCols = {'ssErrMean_deg','ssErrRMS_deg','ssErrP2P_deg','ssErrSlope_degps'};

    traces = cell(n,1);
    for k = 1:n
        results.Trial(k)       = k;
        results.q0{k}          = q0(k,:);
        results.w0_degps{k}    = w0(k,:);
        results.w0Mag_degps(k) = norm(w0(k,:));

        if ~isempty(out(k).ErrorMessage)
            warning('Trial %d failed: %s', k, out(k).ErrorMessage);
            results{k, [{'initErr_deg','peakErr_deg','finalErr_deg'}, ssCols, {'settleTime_s'}]} = NaN;
            results.settled(k) = false;
            continue
        end

        [t, errDeg, rateDeg] = extractSignals(out(k));
        traces{k} = [t errDeg];

        % The first couple of AttitudeError samples are identically zero while
        % the FSW pipeline fills, so they are excluded from every metric below.
        first = find(errDeg > 1e-9, 1, 'first');
        if isempty(first), first = 2; end
        s = first:numel(t);
        results.initErr_deg(k)  = errDeg(s(1));
        results.peakErr_deg(k)  = max(errDeg(s));
        results.finalErr_deg(k) = errDeg(end);

        % Steady-state statistics over the trailing window, never reaching back
        % into the leading zero samples.
        win = t >= (1 - opts.SteadyWindow) * t(end);
        win(1:s(1)-1) = false;
        tailT = t(win);
        tailE = errDeg(win);
        if numel(tailE) >= 2
            results.ssErrMean_deg(k) = mean(tailE);
            results.ssErrRMS_deg(k)  = rms(tailE);
            results.ssErrP2P_deg(k)  = max(tailE) - min(tailE);
            fit = polyfit(tailT - tailT(1), tailE, 1);
            results.ssErrSlope_degps(k) = fit(1);
        else
            results{k, ssCols} = NaN;
        end

        bad = errDeg(s) > opts.SettleTolDeg;
        if isfinite(opts.SettleRateTol)
            bad = bad | rateDeg(s) > opts.SettleRateTol;
        end
        last = find(bad, 1, 'last');
        if isempty(last)
            results.settleTime_s(k) = t(s(1));   % already inside tolerance
            results.settled(k)      = true;
        elseif last == numel(s)
            results.settleTime_s(k) = NaN;       % still outside at stop time
            results.settled(k)      = false;
        else
            results.settleTime_s(k) = t(s(last) + 1);
            results.settled(k)      = true;
        end
    end

    % Timing guard: confirm nothing in this run moved the model's step size.
    timeStepAfter = dict.get('CubeSatTimeStep');
    if ~isequal(timeStepAfter, timeStepAtRun)
        warning(['CubeSatTimeStep changed from %g to %g during this run. The harness ' ...
                 'does not write it, so something else did. Restoring, but treat these ' ...
                 'results as suspect.'], timeStepAtRun, timeStepAfter);
        dict.set('CubeSatTimeStep', timeStepAtRun);
    end

    results.Properties.UserData = struct('options', opts, 'traces', {traces}, ...
                                         'CubeSatTimeStep', timeStepAtRun);

    %% Report
    ok = results.settled;
    fprintf('\n--- Monte Carlo summary (%d trials, %g s each) ---\n', n, opts.StopTime);
    fprintf('settled within tolerance (%.2f deg): %d/%d (%.0f%%)\n', ...
        opts.SettleTolDeg, sum(ok), n, 100*sum(ok)/n);
    if any(ok)
        st = results.settleTime_s(ok);
        fprintf('settling time  mean %.1f s  median %.1f s  95th pct %.1f s  max %.1f s\n', ...
            mean(st), median(st), prctile(st, 95), max(st));
    end
    if any(~ok)
        fprintf('non-settling trials: %s\n', mat2str(results.Trial(~ok)'));
    end

    win = 100 * opts.SteadyWindow;
    fprintf('steady-state over last %.0f%% of run: mean %.4f deg  RMS %.4f deg  p2p %.4f deg\n', ...
        win, mean(results.ssErrMean_deg, 'omitnan'), mean(results.ssErrRMS_deg, 'omitnan'), ...
        mean(results.ssErrP2P_deg, 'omitnan'));
    slope = mean(results.ssErrSlope_degps, 'omitnan');
    fprintf('steady-state drift: %.2e deg/s\n', slope);
    if abs(slope) * opts.StopTime > 0.05 * opts.SettleTolDeg
        fprintf(['  NOTE: the error is still trending, so the values above are a ' ...
                 'reading on a\n        trend, not a converged steady state. ' ...
                 'This orbit''s period is ~5560 s;\n        %g s covers %.1f%% of ' ...
                 'one orbit. Re-run longer before quoting them.\n'], ...
                 opts.StopTime, 100 * opts.StopTime / 5563.6);
    end

    if opts.Plot
        plotResults(results, traces, opts);
    end
end

% ------------------------------------------------------------------------
function cleanup = prepareModel(mdl)
%prepareModel Turn on the signal logging the study needs and disable the
%Satellite Scenario Playback block, whose StopFcn callback errors out on the
%second and later programmatic sim() calls ("Invalid or deleted object").
%Everything is put back by the returned onCleanup object.

    % A model left paused or compiled by an earlier run rejects set_param.
    for attempt = 1:50
        if strcmp(get_param(mdl, 'SimulationStatus'), 'stopped'), break, end
        set_param(mdl, 'SimulationCommand', 'stop');
        pause(0.2);
    end

    % An open block dialog with unapplied edits aborts sim() with
    % "Apply or cancel unapplied changes ...", so close the one we override.
    try
        close_system([mdl '/Vehicle Model/Vehicle Dynamics/Spacecraft Dynamics']);
    catch
    end

    restore = struct('handle', {}, 'params', {});

    acsPort = outPortHandle(mdl, [mdl '/Vehicle Model/Vehicle Flight Software/Attitude Control System'], 1);
    restore(end+1) = saveAndSet(acsPort, ...
        {'DataLogging','DataLoggingNameMode','DataLoggingName'}, {'on','Custom','ACSOut'});

    vd  = [mdl '/Vehicle Model/Vehicle Dynamics'];
    wPort = outPortHandle(mdl, vd, portNumber(vd, 'angular_velocity'));
    restore(end+1) = saveAndSet(wPort, ...
        {'DataLogging','DataLoggingNameMode','DataLoggingName'}, {'on','Custom','omega_body'});

    satScen = [mdl '/Satellite Scenario Playback'];
    if getSimulinkBlockHandle(satScen) > 0
        restore(end+1) = saveAndSet(getSimulinkBlockHandle(satScen), {'Commented'}, {'on'});
    end

    % Model-level params: pass the model name, not a block handle
    % (getSimulinkBlockHandle returns -1 for a block diagram root).
    restore(end+1) = saveAndSet(mdl, ...
        {'SignalLogging','SignalLoggingName'}, {'on','logsout'});

    cleanup = onCleanup(@() undo(restore));

    function undo(r)
        for j = 1:numel(r)
            try
                set_param(r(j).handle, r(j).params{:});
            catch
            end
        end
    end
end

function s = saveAndSet(h, names, values)
    old = cell(1, 2*numel(names));
    for k = 1:numel(names)
        old{2*k-1} = names{k};
        old{2*k}   = get_param(h, names{k});
        set_param(h, names{k}, values{k});
    end
    s = struct('handle', h, 'params', {old});
end

function h = outPortHandle(~, blk, idx)
    ph = get_param(blk, 'PortHandles');
    h  = ph.Outport(idx);
end

function n = portNumber(sys, name)
    b = find_system(sys, 'SearchDepth', 1, 'LookUnderMasks', 'all', ...
                    'BlockType', 'Outport', 'Name', name);
    n = str2double(get_param(b{1}, 'Port'));
end

% ------------------------------------------------------------------------
function q = sampleAttitude(n, opts)
%sampleAttitude Random ICRF->body quaternions, scalar first, q0 >= 0.

    switch opts.AttitudeMode
        case "uniform"
            % Normalised 4-D Gaussian is Haar-uniform on the rotation group.
            q = randn(n, 4);
            q = q ./ vecnorm(q, 2, 2);
        case "axisangle"
            axis  = randn(n, 3);
            axis  = axis ./ vecnorm(axis, 2, 2);
            theta = deg2rad(opts.MinTiltDeg + (opts.MaxTiltDeg - opts.MinTiltDeg) * rand(n,1));
            q     = [cos(theta/2), sin(theta/2) .* axis];
    end
    flip = q(:,1) < 0;
    q(flip,:) = -q(flip,:);
end

function w = sampleRate(n, range)
%sampleRate Random body rates, deg/s: uniform direction, uniform magnitude.

    dir = randn(n, 3);
    dir = dir ./ vecnorm(dir, 2, 2);
    mag = range(1) + (range(2) - range(1)) * rand(n, 1);
    w   = mag .* dir;
end

% ------------------------------------------------------------------------
function [t, errDeg, rateDeg] = extractSignals(simOut)
%extractSignals Pointing error angle and body rate magnitude from one run.
%
%   ACSOut.AttitudeError is a bus of Roll/Pitch/Yaw in RADIANS (the Aerospace
%   "Quaternions to Rotation Angles" block, XYZ order).  The three angles are
%   folded back into a quaternion so the reported error is the true single-axis
%   rotation angle rather than a norm of Euler angles.

    ae = simOut.logsout.get('ACSOut').Values.AttitudeError;
    t  = ae.Roll.Time(:);
    qe = angle2quat(ae.Roll.Data(:), ae.Pitch.Data(:), ae.Yaw.Data(:), 'XYZ');
    errDeg = rad2deg(2 * acos(min(1, abs(qe(:,1)))));

    w = simOut.logsout.get('omega_body').Values;         % deg/s
    rateDeg = interp1(w.Time(:), vecnorm(w.Data, 2, 2), t, 'linear', 'extrap');
end

% ------------------------------------------------------------------------
function plotResults(results, traces, opts)
    figure('Name', 'CubeSat PD Monte Carlo', 'Color', 'w');
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(tl, [1 2]); hold on; grid on
    for k = 1:numel(traces)
        if isempty(traces{k}), continue, end
        c = [0.5 0.5 0.5]; if ~results.settled(k), c = [0.85 0.2 0.2]; end
        plot(traces{k}(:,1), traces{k}(:,2), 'Color', [c 0.45]);
    end
    yline(opts.SettleTolDeg, 'b--', 'LineWidth', 1.2);
    set(gca, 'YScale', 'log');
    xlabel('time (s)'); ylabel('pointing error (deg)');
    title(sprintf('%d trials (red = did not settle)', height(results)));

    nexttile; grid on
    histogram(results.settleTime_s(results.settled));
    xlabel('settling time (s)'); ylabel('count'); title('Settling time');

    nexttile; hold on; grid on
    ok = results.settled;
    scatter(results.w0Mag_degps(ok),  results.settleTime_s(ok), 24, results.initErr_deg(ok), 'filled');
    scatter(results.w0Mag_degps(~ok), repmat(opts.StopTime, sum(~ok), 1), 36, 'rx', 'LineWidth', 1.2);
    cb = colorbar; cb.Label.String = 'initial error (deg)';
    xlabel('|\omega_0| (deg/s)'); ylabel('settling time (s)');
    title('Settling vs initial tumble rate');
end
