function out = readmeStats(file, runIndex, signalName, opts)
%READMESTATS  Print README-ready numbers from a saved SDI run, and plot them.
%
%   readmeStats(file)                    list the runs in the file, then stop
%   readmeStats(file, runIndex)          analyse that run
%   readmeStats(file, runIndex, signal)  override the attitude-error signal
%
%   The figure is two stacked panels sharing a time axis: attitude error on
%   top, gyro bias error below, each autoscaled to itself so that deg and
%   deg/hr keep their own range.  The 1-vector (magnetometer only) windows
%   are shaded.
%
%   NAME-VALUE
%     Plot      Draw the figure.                                  (true)
%     SegmentSec  Segment length handed to attitudeStats.         (2700)
%
%   OUTPUT
%     out.attitude   attitudeStats struct for the error signal
%     out.bias       attitudeStats struct for b_error, if it was logged
%     out.signals    names available in the run
%
%   Example:
%       readmeStats("True_anom_285.mldatx")        % see what is in there
%       readmeStats("True_anom_285.mldatx", 25)
%
%   See also attitudeStats, loadRunData.

    arguments
        file            {mustBeTextScalar}
        runIndex        (1,1) double = -1
        signalName      {mustBeTextScalar} = ""
        opts.Plot       (1,1) logical = true
        opts.SegmentSec (1,1) double = 2700
    end

    Simulink.sdi.load(file);
    ids = Simulink.sdi.getAllRunIDs;

    % ---- no run chosen: list what is on offer and stop --------------------
    if runIndex < 0
        fprintf("\n%d runs in %s\n\n", numel(ids), file);
        fprintf("  %3s  %-34s  %7s  %9s\n", "idx", "name", "signals", "t_end (s)");
        for k = 1:numel(ids)
            r = Simulink.sdi.getRun(ids(k));
            tEnd = NaN;
            if r.SignalCount > 0
                tv = r.getSignalByIndex(1).Values.Time;
                if ~isempty(tv), tEnd = tv(end); end
            end
            fprintf("  %3d  %-34s  %7d  %9.0f\n", k, r.Name, r.SignalCount, tEnd);
        end
        fprintf("\nPick one:  readmeStats(""%s"", <idx>)\n\n", file);
        out = struct();
        return
    end

    % ---- what signals does this run actually have? -----------------------
    r     = Simulink.sdi.getRun(ids(runIndex));
    names = strings(r.SignalCount, 1);
    for k = 1:r.SignalCount
        names(k) = string(r.getSignalByIndex(k).Name);
    end
    out.signals = names;

    if strlength(signalName) == 0
        signalName = firstPresent(["deg_error" "MATLAB Function:2" "AttitudeError" "error"], names);
        if signalName == ""
            fprintf("Signals in run %d:\n", runIndex);  fprintf("   %s\n", names);
            error("readmeStats:NoSignal", ...
                  "No recognised attitude-error signal. Pass one explicitly.");
        end
    end
    fprintf("\nRun %d of %d  (%s)   attitude error = ""%s""\n", ...
            runIndex, numel(ids), r.Name, signalName);

    % ---- attitude error, split by measurement mode -----------------------
    out.attitude = attitudeStats(file, SignalName=signalName, Units="deg", ...
                                 SegmentSec=opts.SegmentSec, RunIndex=runIndex);

    fprintf("\n  ---- README: Attitude knowledge - MUKF ----\n");
    printByMode(out.attitude, "deg");
    if ~isnan(out.attitude.peak.value)
        fprintf("  1-vector peak  %.2f deg at t = %.0f s (segment %d)\n", ...
                out.attitude.peak.value, out.attitude.peak.time, out.attitude.peak.segment);
    end

    % ---- gyro bias error, if it was logged -------------------------------
    biasName = firstPresent(["b_error" "MATLAB Function:5"], names);
    out.bias = [];
    if biasName ~= ""
        out.bias = attitudeStats(file, SignalName=biasName, Units="deg/hr", ...
                                 SegmentSec=opts.SegmentSec, RunIndex=runIndex);
        [tb, yb] = loadRunData(file, biasName, RunIndex=runIndex);
        fprintf("\n  ---- README: Gyro bias ----\n");
        fprintf("  final value        %.4f  (at t = %.0f s)\n", yb(end), tb(end));
        tol = 0.1 * abs(yb(end));
        bad = find(abs(yb - yb(end)) > tol, 1, "last");
        if isempty(bad)
            fprintf("  converged in       < first sample\n");
        elseif bad < numel(tb)
            fprintf("  converged in       %.0f s  (within 10%% of final)\n", tb(bad+1));
        else
            fprintf("  converged in       never settles inside 10%% of final\n");
        end
        printByMode(out.bias, "deg/hr");
    else
        fprintf("\n  (no bias-error signal in this run - skipping the bias table)\n");
    end
    fprintf("\n");

    if opts.Plot
        drawSummary(file, runIndex, names, signalName, biasName, opts.SegmentSec, r.Name);
    end
end

% =========================================================================

function name = firstPresent(candidates, names)
%FIRSTPRESENT  First candidate that exists in NAMES, or "" if none do.
    hit  = candidates(ismember(candidates, names));
    name = "";
    if ~isempty(hit), name = hit(1); end
end

function printByMode(s, units)
    modes = fieldnames(s.byMode);
    for k = 1:numel(modes)
        m = s.byMode.(modes{k});
        if ~isfield(m, "rms")
            fprintf("  %-12s  (no samples)\n", modes{k});  continue
        end
        fprintf("  %-12s  RMS %8.3f %-7s 3-sigma %8.3f   max %8.3f   n=%d\n", ...
                modes{k}, m.rms, units, 3*m.std, m.max, m.n);
    end
end

% -------------------------------------------------------------------------

function drawSummary(file, runIndex, names, signalName, biasName, segSec, runName)
%DRAWSUMMARY  deg_error and b_error against time, one stacked panel each so
%   that deg and deg/hr keep their own scale.  The 1-vector (magnetometer
%   only) windows are shaded.

    [t, e] = loadRunData(file, signalName, RunIndex=runIndex);

    hasBias = biasName ~= "";
    if hasBias
        [tb, yb] = loadRunData(file, biasName, RunIndex=runIndex);
    end

    % Real mode flag if the run logged one, otherwise fall back to the same
    % alternating-segment assumption attitudeStats makes.
    modeName = firstPresent(["Mode:Value" "Mode" "MATLAB Function:6"], names);
    if modeName ~= ""
        [tm, ym] = loadRunData(file, modeName, RunIndex=runIndex);
    else
        tm = t;
        ym = double(mod(floor((t - t(1)) / segSec), 2) == 1);
    end

    f  = figure(Name=sprintf("readmeStats - %s", runName), Color="w");
    tl = tiledlayout(f, 1 + hasBias, 1, TileSpacing="compact", Padding="compact");
    title(tl, sprintf("%s  (run %d)   shaded = 1-vector", runName, runIndex), ...
          Interpreter="none");

    ax1 = plotTrace(nexttile(tl), t, e, tm, ym, signalName, "attitude error (deg)");
    xlabel(ax1, "time (s)");

    if hasBias
        ax2 = plotTrace(nexttile(tl), tb, yb, tm, ym, biasName, "bias error (deg/hr)");
        xlabel(ax2, "time (s)");
        xlabel(ax1, "");
        linkaxes([ax1 ax2], "x");
    end
end

function ax = plotTrace(ax, t, y, tm, ym, sigName, ylab)
%PLOTTRACE  One signal, autoscaled to itself, with the 1-vector windows shaded.
    hold(ax, "on");  grid(ax, "on");

    % Pin the y range from the data first, so the shading patches cannot
    % drag the axis limits out to their own corner coordinates.
    yl  = [min(y) max(y)];
    pad = 0.05 * max(diff(yl), eps);
    yl  = yl + [-pad pad];

    on = ym(:) > 0.5;
    d  = diff([false; on; false]);
    lo = find(d ==  1);
    hi = find(d == -1) - 1;
    for k = 1:numel(lo)
        x = [tm(lo(k)) tm(hi(k)) tm(hi(k)) tm(lo(k))];
        patch(ax, x, [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
              EdgeColor="none", FaceAlpha=0.45, HandleVisibility="off");
    end

    plot(ax, t, y, LineWidth=1.2);
    ylim(ax, yl);
    ylabel(ax, ylab);
    title(ax, sigName, Interpreter="none", FontWeight="normal");
end
