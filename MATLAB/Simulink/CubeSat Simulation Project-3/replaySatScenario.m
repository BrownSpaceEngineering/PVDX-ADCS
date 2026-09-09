function [scenario, viewer] = replaySatScenario(source, opts)
%REPLAYSATSCENARIO  Open the Satellite Scenario Viewer on a past run.
%
%   The Satellite Scenario Playback block builds its scenario inside its
%   own StopFcn, and only when the model is started interactively -- sim()
%   never leaves the four exSatScen* variables where the callback can find
%   them, and a run reloaded from the Simulation Data Inspector has no
%   callback to fire at all. This rebuilds the scenario from the logged
%   trajectory instead, so any saved run can be played back.
%
%   replaySatScenario("../results/PD_with_w_true.mldatx")
%   replaySatScenario(simOut)
%   replaySatScenario(simOut, SampleTime=10, Play=false)
%
%   Call with no signal available to check first:
%
%       loadRunData("../results/PD_with_w_true.mldatx")
%
%   The run must contain the four trajectory signals the playback block is
%   fed: utc_JD, X_ecef, V_ecef and q_ecef2b, tapped off the vehicle state
%   bus. Note that runCubeSatMonteCarlo comments the playback block out and
%   logs only ACSOut and omega_body, so its runs cannot be replayed -- log
%   the bus signals, or mark those four lines for logging, before the run
%   you want to visualise.
%
%   INPUT
%     source      Simulink.SimulationOutput, Simulink.SimulationData.Dataset,
%                 or the path to a saved SDI file (.mldatx)
%
%   NAME-VALUE
%     RunIndex        Which run to read from a multi-run .mldatx file.
%                     -1 selects the last run in the file.           (-1)
%     SampleTime      Scenario step, s. The logged trajectory is resampled
%                     onto this uniform grid.                        (60)
%     StartTime       UTC epoch to use when the run has no utc_JD
%                     signal.                                       (NaT)
%     TStart, TEnd    Trim the replay to this slice of the run.  (0, inf)
%     PositionUnits   "m" or "km" for the logged position.          ("m")
%     Name            Satellite label in the viewer.           ("CubeSat")
%     Play            Start playback once the viewer opens.        (true)
%     UseBlock        Build through the playback block, so its mask
%                     settings (FOV, 3-D model, triad) are honoured. Falls
%                     back to a plain scenario if the block is missing or
%                     errors.                                      (true)
%     Block           Path to the playback block.
%                                 ("asbCubeSat/Satellite Scenario Playback")
%     TimeSignal, PositionSignal, VelocitySignal, AttitudeSignal
%                     Override signal-name lookup when the run logged them
%                     under different names.
%
%   OUTPUTS
%     scenario   satelliteScenario object, also left in the base workspace
%                as exSatScenario to match the block's own behaviour
%     viewer     satelliteScenarioViewer handle
%
%   See also loadRunData, attitudeStats, runCubeSatMonteCarlo.

    arguments
        source
        opts.RunIndex        (1,1) double = -1
        opts.SampleTime      (1,1) double {mustBePositive} = 60
        opts.StartTime       (1,1) datetime = NaT
        opts.TStart          (1,1) double = 0
        opts.TEnd            (1,1) double = inf
        opts.PositionUnits   {mustBeMember(opts.PositionUnits, ["m" "km"])} = "m"
        opts.Name            {mustBeTextScalar} = "CubeSat"
        opts.Play            (1,1) logical = true
        opts.UseBlock        (1,1) logical = true
        opts.Block           {mustBeTextScalar} = "asbCubeSat/Satellite Scenario Playback"
        opts.TimeSignal      {mustBeTextScalar} = ""
        opts.PositionSignal  {mustBeTextScalar} = ""
        opts.VelocitySignal  {mustBeTextScalar} = ""
        opts.AttitudeSignal  {mustBeTextScalar} = ""
    end

    % ---- project, so this runs from a cold MATLAB ------------------------
    try
        currentProject;
    catch
        here = fileparts(mfilename("fullpath"));
        openProject(fullfile(here, "CubeSatSimulationProject.prj"));
    end

    % ---- resolve the source once ----------------------------------------
    % One pass only: re-opening a .mldatx per signal would push a duplicate
    % run into SDI on every call.
    [names, fetch] = openSource(source, opts.RunIndex);

    iTime = pickSignal(names, opts.TimeSignal, ...
                       ["utc_JD" "utcJD" "exSatScenTimeData" "JD"], "time", false);
    iPos  = pickSignal(names, opts.PositionSignal, ...
                       ["X_ecef" "exSatScenPositionData" "r_ecef" "X_eci" "X_icrf"], "position", true);
    iVel  = pickSignal(names, opts.VelocitySignal, ...
                       ["V_ecef" "exSatScenVelocityData" "v_ecef" "V_eci" "V_icrf"], "velocity", false);
    iAtt  = pickSignal(names, opts.AttitudeSignal, ...
                       ["q_ecef2b" "exSatScenAttitudeData" "q_ecef2body" "q_i2b"], "attitude", false);

    [t, pos] = getSignal(fetch, iPos, 3);
    if opts.PositionUnits == "km"
        pos = pos * 1000;
    end

    vel = [];
    if ~isnan(iVel)
        [~, vel] = getSignal(fetch, iVel, 3);
        if opts.PositionUnits == "km"
            vel = vel * 1000;
        end
    end

    att = [];
    if ~isnan(iAtt)
        [~, att] = getSignal(fetch, iAtt, 4);
    end

    jd = [];
    if ~isnan(iTime)
        [~, jd] = getSignal(fetch, iTime, 1);
    elseif isnat(opts.StartTime)
        error("replaySatScenario:NoEpoch", ...
              "The run has no utc_JD signal, so the UTC epoch is unknown. " + ...
              "Pass StartTime=datetime(...) or log utc_JD.");
    end

    % ---- frames, taken from whichever signal name matched ----------------
    posFrame = frameOf(names{iPos}, "ecef");
    attFrame = "ecef";
    if ~isnan(iAtt)
        attFrame = frameOf(names{iAtt}, "ecef");
    end

    % ---- uniform resample; satelliteScenario needs an even grid ----------
    keep = t >= opts.TStart & t <= opts.TEnd;
    if nnz(keep) < 2
        error("replaySatScenario:EmptySlice", ...
              "TStart/TEnd leave fewer than two samples of a run spanning %g to %g s.", ...
              t(1), t(end));
    end
    t = t(keep);
    pos = pos(keep,:);
    if ~isempty(vel), vel = vel(keep,:); end
    if ~isempty(att), att = att(keep,:); end
    if ~isempty(jd),  jd  = jd(keep,:);  end

    % A variable-step run repeats a time at every solver reset, which
    % interp1 rejects outright. Keep the first sample of each instant.
    [t, first] = unique(t);
    pos = pos(first,:);
    if ~isempty(vel), vel = vel(first,:); end
    if ~isempty(att), att = att(first,:); end
    if ~isempty(jd),  jd  = jd(first,:);  end

    tGrid = (t(1) : opts.SampleTime : t(end)).';
    if numel(tGrid) < 2
        error("replaySatScenario:SampleTimeTooLarge", ...
              "SampleTime=%g s exceeds the %g s span of the run.", ...
              opts.SampleTime, t(end)-t(1));
    end
    posGrid = interp1(t, pos, tGrid, "linear");
    velGrid = [];
    if ~isempty(vel), velGrid = interp1(t, vel, tGrid, "linear"); end
    attGrid = [];
    if ~isempty(att)
        attGrid = interp1(t, att, tGrid, "linear");
        % Componentwise interpolation of a quaternion leaves it slightly off
        % the unit sphere; pointAt wants it normalised.
        attGrid = attGrid ./ vecnorm(attGrid, 2, 2);
    end

    if isempty(jd)
        epoch = opts.StartTime;
        if ~isempty(epoch.TimeZone)
            epoch.TimeZone = "UTC";   % express the instant in UTC,
            epoch.TimeZone = "";      % then drop the zone, as the scenario wants
        end
        times = epoch + seconds(tGrid - tGrid(1));
        jdGrid = juliandate(times);
    else
        jdGrid = interp1(t, jd, tGrid, "linear");
        times  = datetime(jdGrid, ConvertFrom="juliandate");
    end

    % ---- preferred path: let the block build it, mask settings and all ---
    scenario = [];
    if opts.UseBlock
        scenario = buildViaBlock(opts.Block, tGrid, jdGrid, posGrid, velGrid, attGrid);
    end

    % ---- fallback: a plain scenario from the same numbers ----------------
    if isempty(scenario)
        scenario = satelliteScenario(times(1), times(end), opts.SampleTime);
        if isempty(velGrid)
            tt = timetable(times(:), posGrid);
        else
            tt = timetable(times(:), posGrid, velGrid);
        end
        sat = satellite(scenario, tt, ...
                        CoordinateFrame=posFrame, Name=opts.Name);
        if ~isempty(attGrid)
            pointAt(sat, timetable(times(:), attGrid), ...
                    CoordinateFrame=attFrame, Format="quaternion", ...
                    ExtrapolationMethod="nadir");
        end
    end

    % The block leaves the scenario in the base workspace under this name;
    % keep that contract so the mask's mission-analysis advice still holds.
    assignin("base", "exSatScenario", scenario);

    viewer = satelliteScenarioViewer(scenario, ShowDetails=true);
    if opts.Play
        play(scenario);
    end

    fprintf("Replaying %d samples, %s to %s UTC, %g s step.\n", ...
            numel(tGrid), string(times(1)), string(times(end)), opts.SampleTime);
end

% -------------------------------------------------------------------------

function scenario = buildViaBlock(blk, tGrid, jdGrid, pos, vel, att)
%BUILDVIABLOCK  Call the playback block's own scenario builder.
%   Returns [] if the block is unavailable or the call fails, so the caller
%   can fall back. The builder reads mask settings off the block, so the
%   model has to be loaded and the four inputs have to look like the To
%   Workspace timeseries it normally receives.
    scenario = [];
    if isempty(vel) || isempty(att)
        return   % builder takes all four; the fallback handles partial runs
    end
    try
        mdl = extractBefore(string(blk), "/");
        if ~bdIsLoaded(mdl)
            load_system(mdl);
        end
        if getSimulinkBlockHandle(char(blk)) <= 0
            return
        end
        scenario = Aero.spacecraft.analysis.createSatelliteScenarioFromSim( ...
            char(blk), ...
            timeseries(jdGrid, tGrid), ...
            timeseries(pos,    tGrid), ...
            timeseries(vel,    tGrid), ...
            timeseries(att,    tGrid));
    catch err
        warning("replaySatScenario:BlockBuildFailed", ...
                "Could not build through %s (%s). Falling back to a plain " + ...
                "satelliteScenario.", blk, err.message);
        scenario = [];
    end
end

function idx = pickSignal(names, override, candidates, role, required)
%PICKSIGNAL  Locate a signal by name, then by case-insensitive candidate.
    idx = NaN;
    if strlength(override) > 0
        idx = find(strcmp(names, override), 1);
        if isempty(idx)
            error("replaySatScenario:SignalNotFound", ...
                  "No signal named ""%s"". Available:\n  %s", ...
                  override, strjoin(names, sprintf("\n  ")));
        end
        return
    end
    for c = candidates
        hit = find(strcmpi(names, c), 1);
        if ~isempty(hit)
            idx = hit;
            return
        end
    end
    if required
        error("replaySatScenario:SignalNotFound", ...
              "No %s signal found (looked for %s). The run holds:\n  %s\n" + ...
              "Log the trajectory bus, or pass PositionSignal= explicitly.", ...
              role, strjoin(candidates, ", "), strjoin(names, sprintf("\n  ")));
    end
    warning("replaySatScenario:MissingSignal", ...
            "No %s signal in this run; continuing without it.", role);
end

function [t, y] = getSignal(fetch, idx, width)
    ts = fetch(idx);
    t  = ts.Time(:);
    y  = squeeze(ts.Data);
    if size(y,1) ~= numel(t) && size(y,2) == numel(t)
        y = y.';     % a matrix signal squeezes to width-by-N; put time down
    end
    if size(y,2) ~= width
        error("replaySatScenario:BadWidth", ...
              "Expected a %d-wide signal, got %d columns.", width, size(y,2));
    end
end

function frame = frameOf(name, default)
    n = lower(string(name));
    if contains(n, "ecef") || contains(n, "fixed")
        frame = "ecef";
    elseif contains(n, "eci") || contains(n, "icrf") || contains(n, "inertial")
        frame = "inertial";
    else
        frame = default;
    end
end

function [names, fetch] = openSource(source, runIndex)
    if isa(source, "Simulink.SimulationOutput")
        [names, fetch] = fromDataset(source.logsout);
    elseif isa(source, "Simulink.SimulationData.Dataset")
        [names, fetch] = fromDataset(source);
    elseif ischar(source) || isstring(source)
        [names, fetch] = fromSDIFile(string(source), runIndex);
    else
        error("replaySatScenario:BadSource", ...
              "source must be a SimulationOutput, a Dataset, or a path to " + ...
              "a .mldatx file. Got a %s.", class(source));
    end
end

function [names, fetch] = fromDataset(ds)
    n     = ds.numElements;
    names = cell(n,1);
    for k = 1:n
        names{k} = ds.getElement(k).Name;
    end
    fetch = @(k) ds.getElement(k).Values;
end

function [names, fetch] = fromSDIFile(file, runIndex)
    if ~isfile(file)
        error("replaySatScenario:FileNotFound", "No such file: %s", file);
    end
    Simulink.sdi.load(file);

    ids = Simulink.sdi.getAllRunIDs;
    if isempty(ids)
        error("replaySatScenario:NoRuns", "%s contains no runs.", file);
    end
    if runIndex < 0
        id = ids(end);
    elseif runIndex > numel(ids)
        error("replaySatScenario:BadRunIndex", ...
              "RunIndex %d requested but the file holds %d runs.", ...
              runIndex, numel(ids));
    else
        id = ids(runIndex);
    end

    runObj = Simulink.sdi.getRun(id);
    n      = runObj.SignalCount;
    names  = cell(n,1);
    for k = 1:n
        names{k} = runObj.getSignalByIndex(k).Name;
    end
    fetch = @(k) runObj.getSignalByIndex(k).Values;
end
