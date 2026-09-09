function stats = attitudeStats(source, opts)
%ATTITUDESTATS  Segment-wise error statistics for a logged run.
%
%   Splits a run into fixed-length segments (default 2700 s), labels them
%   alternately by measurement mode starting with the first, and reports
%   statistics per segment plus pooled per mode.
%
%   stats = attitudeStats()                     simulate, then analyse
%   stats = attitudeStats(simOut)               analyse an existing run
%   stats = attitudeStats("run.mldatx")         analyse a saved SDI file
%
%   Attitude error:
%       s = attitudeStats(simOut);
%
%   Gyro bias error (same maths, different signal and units):
%       b = attitudeStats(simOut, SignalName="b_error", Units="deg/hr");
%
%   NAME-VALUE
%     SignalName   Logged signal to read.                    ("deg_error")
%     SegmentSec   Segment length in seconds.                      (2700)
%     ModeNames    Labels for the alternating modes, first one
%                  applying to the first segment.  (["2-vector" "1-vector"])
%     TStart       Ignore samples before this time. Segment boundaries
%                  still anchor at the start of the run.               (0)
%     Percentiles  Percentiles to report.                        ([95 99])
%     Units        Units string used in the printout.             ("deg")
%     Plot         Draw the trace with segment boundaries.        (false)
%     Model        Model to simulate when SOURCE is omitted. ("asbCubeSat")
%
%   OUTPUT
%     stats.segments   table, one row per segment
%     stats.byMode     struct keyed by mode name, pooled over all its
%                      segments: rms, mean, std, max, percentiles, n
%     stats.peak       worst sample in the second mode (the eclipse peak):
%                      .value, .time, .segment
%     stats.signal, stats.units, stats.segmentSec
%
%   Percentiles are computed by sorting, so no Statistics Toolbox needed.
%
%   See also loadRunData, runCubeSatMonteCarlo.

    arguments
        source                          = []
        opts.SignalName  {mustBeTextScalar} = "MATLAB Function:2"
        opts.SegmentSec  (1,1) double   = 2700
        opts.ModeNames   (1,2) string   = ["2-vector" "1-vector"]
        opts.TStart      (1,1) double   = 0
        opts.Percentiles (1,:) double   = [95 99]
        opts.Units       {mustBeTextScalar} = "deg"
        opts.Plot        (1,1) logical  = false
        opts.Model       {mustBeTextScalar} = "asbCubeSat"
        opts.RunIndex    (1,1) double   = -1
    end

    % ---- project, so this runs from a cold MATLAB ------------------------
    try
        currentProject;
    catch
        here = fileparts(mfilename("fullpath"));
        openProject(fullfile(here, "CubeSatSimulationProject.prj"));
    end

    % ---- simulate if no source given -------------------------------------
    if isempty(source)
        if ~bdIsLoaded(opts.Model)
            load_system(opts.Model);
        end
        fprintf("Simulating %s ...\n", opts.Model);
        source = sim(Simulink.SimulationInput(opts.Model));
    end

    % ---- pull the signal --------------------------------------------------
    [t, e] = loadRunData(source, opts.SignalName, RunIndex=opts.RunIndex);
    if size(e,2) ~= 1
        error("attitudeStats:NotScalar", ...
              "Signal ""%s"" has %d columns; expected a scalar error.", ...
              opts.SignalName, size(e,2));
    end
    e = double(e);
    t = double(t);

    % ---- segment boundaries, anchored at the start of the run -------------
    % TStart trims samples but does NOT shift the boundaries, so segment k
    % always covers the same wall-clock interval regardless of trimming.
    t0   = t(1);
    span = t(end) - t0;
    nSeg = max(1, ceil(span / opts.SegmentSec));

    segIdx  = zeros(nSeg,1);
    segMode = strings(nSeg,1);
    tLo     = zeros(nSeg,1);
    tHi     = zeros(nSeg,1);
    nPts    = zeros(nSeg,1);
    rmsV    = nan(nSeg,1);
    meanV   = nan(nSeg,1);
    stdV    = nan(nSeg,1);
    maxV    = nan(nSeg,1);
    tPeak   = nan(nSeg,1);
    pctV    = nan(nSeg, numel(opts.Percentiles));

    modeOfSeg = strings(nSeg,1);
    keepMask  = false(numel(t), nSeg);

    for k = 1:nSeg
        lo = t0 + (k-1)*opts.SegmentSec;
        hi = lo + opts.SegmentSec;

        if k == nSeg
            m = t >= lo & t <= hi;      % last segment closes inclusively
        else
            m = t >= lo & t <  hi;
        end
        m = m & (t >= opts.TStart);

        % alternate modes, first segment gets ModeNames(1)
        modeOfSeg(k) = opts.ModeNames(2 - mod(k,2));

        segIdx(k)  = k;
        segMode(k) = modeOfSeg(k);
        tLo(k)     = lo;
        tHi(k)     = hi;
        nPts(k)    = nnz(m);
        keepMask(:,k) = m;

        if nPts(k) < 2
            continue                     % leave stats as NaN
        end

        ek = e(m);
        tk = t(m);
        s  = computeStats(ek, opts.Percentiles);

        rmsV(k)   = s.rms;
        meanV(k)  = s.mean;
        stdV(k)   = s.std;
        maxV(k)   = s.max;
        pctV(k,:) = s.pct;

        [~, iMax] = max(ek);
        tPeak(k)  = tk(iMax);
    end

    T = table(segIdx, segMode, tLo, tHi, nPts, rmsV, meanV, stdV, maxV, tPeak, ...
        'VariableNames', {'segment','mode','tStart','tEnd','n', ...
                          'rms','mean','std','max','tPeak'});
    for j = 1:numel(opts.Percentiles)
        T.(sprintf('p%g', opts.Percentiles(j))) = pctV(:,j);
    end
    stats.segments = T;

    % ---- pool across all segments of each mode ---------------------------
    % Pool the raw samples rather than averaging the per-segment RMS values,
    % which would be wrong whenever the segments hold different sample counts.
    stats.byMode = struct();
    for mn = opts.ModeNames
        sel = false(numel(t),1);
        for k = find(modeOfSeg == mn)'
            sel = sel | keepMask(:,k);
        end
        field = matlab.lang.makeValidName(mn);
        if nnz(sel) < 2
            stats.byMode.(field) = struct('n', nnz(sel));
            continue
        end
        s = computeStats(e(sel), opts.Percentiles);
        s.n    = nnz(sel);
        s.mode = mn;
        stats.byMode.(field) = s;
    end

    % ---- peak in the second mode (eclipse) --------------------------------
    stats.peak = struct('value', NaN, 'time', NaN, 'segment', NaN);
    eclipseSegs = find(modeOfSeg == opts.ModeNames(2))';
    bestVal = -inf;
    for k = eclipseSegs
        if nPts(k) < 1, continue, end
        if maxV(k) > bestVal
            bestVal = maxV(k);
            stats.peak.value   = maxV(k);
            stats.peak.time    = tPeak(k);
            stats.peak.segment = k;
        end
    end

    stats.signal     = string(opts.SignalName);
    stats.units      = string(opts.Units);
    stats.segmentSec = opts.SegmentSec;

    % ---- report -----------------------------------------------------------
    u = opts.Units;
    fprintf("\n  %s  —  %g s segments  (%s)\n", ...
            opts.SignalName, opts.SegmentSec, u);
    fprintf("  %s\n", repmat('-', 1, 76));
    fprintf("  %3s  %-10s  %8s  %8s  %8s  %8s  %8s  %9s\n", ...
            "seg", "mode", "t0 (s)", "n", "RMS", "p95", "max", "t_peak (s)");
    for k = 1:nSeg
        if nPts(k) < 2
            fprintf("  %3d  %-10s  %8.0f  %8d   (too few samples)\n", ...
                    k, segMode(k), tLo(k), nPts(k));
            continue
        end
        p95col = pctV(k, find(opts.Percentiles == 95, 1));
        if isempty(p95col), p95col = NaN; end
        fprintf("  %3d  %-10s  %8.0f  %8d  %8.3f  %8.3f  %8.3f  %9.0f\n", ...
                k, segMode(k), tLo(k), nPts(k), rmsV(k), p95col, ...
                maxV(k), tPeak(k));
    end

    fprintf("\n  pooled by mode\n");
    for mn = opts.ModeNames
        field = matlab.lang.makeValidName(mn);
        s = stats.byMode.(field);
        if ~isfield(s, 'rms')
            fprintf("  %-10s  (no samples)\n", mn);
            continue
        end
        fprintf("  %-10s  n=%-7d  RMS %7.3f %s   mean %7.3f   max %7.3f\n", ...
                mn, s.n, s.rms, u, s.mean, s.max);
        for j = 1:numel(opts.Percentiles)
            fprintf("  %-10s  %-9s  p%-4g %7.3f %s\n", "", "", ...
                    opts.Percentiles(j), s.pct(j), u);
        end
    end

    if ~isnan(stats.peak.value)
        fprintf("\n  peak in %s: %.3f %s at t = %.0f s (segment %d)\n\n", ...
                opts.ModeNames(2), stats.peak.value, u, ...
                stats.peak.time, stats.peak.segment);
    else
        fprintf("\n  no %s segments in this run\n\n", opts.ModeNames(2));
    end

    % ---- optional plot ----------------------------------------------------
    if opts.Plot
        figure("Name", opts.SignalName);
        plot(t, e, "LineWidth", 1); hold on
        for k = 1:nSeg
            xline(tLo(k), "-", sprintf("%s", segMode(k)), ...
                  "Color", [0.6 0.6 0.6], "LabelOrientation", "horizontal");
        end
        if ~isnan(stats.peak.value)
            plot(stats.peak.time, stats.peak.value, "rv", ...
                 "MarkerFaceColor", "r");
            text(stats.peak.time, stats.peak.value, ...
                 sprintf("  peak %.2f %s", stats.peak.value, u), ...
                 "VerticalAlignment", "bottom");
        end
        grid on
        xlabel("time (s)")
        ylabel(sprintf("%s (%s)", opts.SignalName, u))
        title(sprintf("%s — %g s segments", opts.SignalName, opts.SegmentSec))
        hold off
    end
end

% =========================================================================

function s = computeStats(x, pcts)
%COMPUTESTATS  RMS, mean, std, max and percentiles of a vector.
%   Percentiles come from sorting, so no Statistics Toolbox is required.
    x      = x(:);
    s.rms  = sqrt(mean(x.^2));
    s.mean = mean(x);
    s.std  = std(x);
    s.max  = max(x);

    xs     = sort(x);
    s.pct  = zeros(1, numel(pcts));
    for j = 1:numel(pcts)
        s.pct(j) = xs(min(numel(xs), max(1, ceil(pcts(j)/100 * numel(xs)))));
    end
end