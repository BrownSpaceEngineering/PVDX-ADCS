function out = readmeStats(file, runIndex, signalName, opts)
%READMESTATS  Print README-ready numbers from a saved SDI run, and plot them.
%
%   readmeStats(file)                    list the runs in the file, then stop
%   readmeStats(file, runIndex)          analyse that run
%   readmeStats(file, -1)                analyse the last run in the file
%   readmeStats(file, runIndex, signal)  override the attitude-error signal
%
%   The figure plots the chosen signal against time, with the 1-vector
%   (magnetometer only) windows shaded.
%
%   NAME-VALUE
%     Plot      Draw the figure.                                  (true)
%     SegmentSec  Segment length handed to attitudeStats.         (2700)
%     TStart    Ignore samples before this time (s), in both the
%               stats and the plot.                                    (0)
%     Title     Plot title.  Defaults to a name for the signal, e.g.
%               "Gyro bias error" for b_error.
%     YMax      Top of the y axis; Inf autoscales to the data.       (Inf)
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
        opts.TStart     (1,1) double = 0
        opts.Title      {mustBeTextScalar} = ""
        opts.YMax       (1,1) double {mustBePositive} = Inf
    end

    Simulink.sdi.load(file);
    ids = Simulink.sdi.getAllRunIDs;

    % ---- no run chosen: list what is on offer and stop --------------------
    % Only when runIndex was omitted.  An explicit negative index means
    % "the last run", matching loadRunData's convention.
    if nargin < 2
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

    % ---- resolve a negative index to the last run ------------------------
    if runIndex < 0
        runIndex = numel(ids);
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
    [defaultTitle, units] = describeSignal(signalName);
    if strlength(opts.Title) == 0
        opts.Title = defaultTitle;
    end
    fprintf("\nRun %d of %d  (%s)   signal = ""%s""\n", ...
            runIndex, numel(ids), r.Name, signalName);

    % ---- chosen signal, split by measurement mode ------------------------
    out.attitude = attitudeStats(file, SignalName=signalName, Units=units, ...
                                 SegmentSec=opts.SegmentSec, TStart=opts.TStart, ...
                                 RunIndex=runIndex);

    fprintf("\n  ---- README: %s ----\n", opts.Title);
    printByMode(out.attitude, units);
    if ~isnan(out.attitude.peak.value)
        fprintf("  1-vector peak  %.2f %s at t = %.0f s (segment %d)\n", ...
                out.attitude.peak.value, units, out.attitude.peak.time, out.attitude.peak.segment);
    end

    % ---- gyro bias error, if it was logged -------------------------------
    biasName = firstPresent(["b_error" "MATLAB Function:5"], names);
    out.bias = [];
    if biasName ~= "" && biasName ~= signalName
        out.bias = attitudeStats(file, SignalName=biasName, Units="deg/s", ...
                                 SegmentSec=opts.SegmentSec, TStart=opts.TStart, ...
                                 RunIndex=runIndex);
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
        printByMode(out.bias, "deg/s");
    else
        fprintf("\n  (no bias-error signal in this run - skipping the bias table)\n");
    end
    fprintf("\n");

    if opts.Plot
        drawSummary(file, runIndex, names, signalName, units, opts);
    end
end

% =========================================================================

function name = firstPresent(candidates, names)
%FIRSTPRESENT  First candidate that exists in NAMES, or "" if none do.
    hit  = candidates(ismember(candidates, names));
    name = "";
    if ~isempty(hit), name = hit(1); end
end

function [label, units] = describeSignal(name)
%DESCRIBESIGNAL  Readable title and units for a logged signal name.
    switch char(name)
        case {'b_error', 'MATLAB Function:5'}
            label = "Gyro bias error";  units = "deg/s";
        case {'deg_error', 'MATLAB Function:2', 'AttitudeError', 'error'}
            label = "Attitude error";   units = "deg";
        otherwise
            label = string(name);       units = "deg";
    end
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

function drawSummary(file, runIndex, names, signalName, units, opts)
%DRAWSUMMARY  The chosen signal against time, with the 1-vector
%   (magnetometer only) windows shaded.

    [t, e] = loadRunData(file, signalName, RunIndex=runIndex);

    % Real mode flag if the run logged one, otherwise fall back to the same
    % alternating-segment assumption attitudeStats makes.
    modeName = firstPresent(["Mode:Value" "Mode" "MATLAB Function:6"], names);
    if modeName ~= ""
        [tm, ym] = loadRunData(file, modeName, RunIndex=runIndex);
    else
        tm = t;
        ym = double(mod(floor((t - t(1)) / opts.SegmentSec), 2) == 1);
    end

    % Drop samples before TStart only now, so the fallback segments above
    % stay anchored at the start of the run, as they are in attitudeStats.
    k  = t  >= opts.TStart;   t  = t(k);   e  = e(k, :);
    k  = tm >= opts.TStart;   tm = tm(k);  ym = ym(k);

    f  = figure(Name=sprintf("readmeStats - %s", opts.Title), Color="w");
    ax = plotTrace(axes(f), t, e, tm, ym, sprintf("error (%s)", units), opts.YMax);
    title(ax, sprintf("%s", opts.Title), ...
          Interpreter="none");
    xlabel(ax, "time (s)");
end

function ax = plotTrace(ax, t, y, tm, ym, ylab, yMax)
%PLOTTRACE  One signal with the 1-vector windows shaded.  The y axis
%   autoscales to the data unless YMAX is finite, which fixes the top.
    hold(ax, "on");  grid(ax, "on");

    % Pin the y range first, so the shading patches cannot drag the axis
    % limits out to their own corner coordinates.
    top = max(y);
    if isfinite(yMax), top = yMax; end
    yl  = [min(y) top];
    pad = 0.05 * max(diff(yl), eps);
    yl(1) = yl(1) - pad;
    if ~isfinite(yMax), yl(2) = yl(2) + pad; end

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
end
