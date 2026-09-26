function [t, y, names] = loadRunData(source, signalName, opts)
%LOADRUNDATA  Pull one logged signal out of a simulation run as (t, y).
%
%   Accepts a live simulation output, a logsout Dataset, or a saved
%   Simulation Data Inspector file, and returns a clean time vector and
%   data matrix with the singleton dimensions removed, the orientation
%   fixed, and the startup transient trimmed off.
%
%   [t, y] = loadRunData(simOut, "att_err")
%   [t, y] = loadRunData("PD_with_w_true.mldatx", "att_err")
%   [t, y] = loadRunData(simOut, "att_err", TStart=300)
%
%   Call with no signal name to list what the source actually contains:
%
%       loadRunData(simOut)
%
%   INPUTS
%     source      Simulink.SimulationOutput, Simulink.SimulationData.Dataset,
%                 or the path to a saved SDI file (.mldatx)
%     signalName  Name of the logged signal. Omit to list available names.
%
%   NAME-VALUE
%     TStart      Discard samples before this time. Use it to skip the
%                 convergence transient before computing statistics.  (0)
%     TEnd        Discard samples after this time.                  (inf)
%     RunIndex    Which run to read from a multi-run .mldatx file.
%                 -1 selects the last run in the file.               (-1)
%
%   OUTPUTS
%     t       N-by-1 time vector
%     y       N-by-M data, time running down the rows
%     names   All signal names available in the source
%
%   See also runCubeSatMonteCarlo, Error.

    arguments
        source
        signalName    {mustBeTextScalar} = ""
        opts.TStart   (1,1) double = 0
        opts.TEnd     (1,1) double = inf
        opts.RunIndex (1,1) double = -1
    end

    % Resolve whatever we were handed into a name list and a fetch handle.
    if isa(source, "Simulink.SimulationOutput")
        [names, fetch] = fromDataset(source.logsout);
    elseif isa(source, "Simulink.SimulationData.Dataset")
        [names, fetch] = fromDataset(source);
    elseif ischar(source) || isstring(source)
        [names, fetch] = fromSDIFile(string(source), opts.RunIndex);
    else
        error("loadRunData:BadSource", ...
              "source must be a SimulationOutput, a Dataset, or a path " + ...
              "to a .mldatx file. Got a %s.", class(source));
    end

    % No signal requested: report what is on offer and stop.
    if strlength(signalName) == 0
        fprintf("%d signals available:\n", numel(names));
        fprintf("  %s\n", names{:});
        t = []; y = [];
        return
    end

    idx = find(strcmp(names, signalName), 1);
    if isempty(idx)
        error("loadRunData:SignalNotFound", ...
              "No signal named ""%s"". Available:\n  %s", ...
              signalName, strjoin(names, sprintf("\n  ")));
    end

    ts = fetch(idx);
    t  = ts.Time(:);
    y  = squeeze(ts.Data);

    % A quaternion logged as a matrix signal comes back 4-by-1-by-N and
    % squeezes to 4-by-N, which silently breaks any statistic taken down
    % the columns. Put time back on dimension 1.
    if ~ismatrix(y)
        error("loadRunData:Unsupported", ...
              "Signal ""%s"" is %s after squeeze; this helper handles " + ...
              "scalar and vector signals only.", ...
              signalName, join(string(size(y)), "-by-"));
    end
    if size(y,1) ~= numel(t) && size(y,2) == numel(t)
        y = y.';
    end

    keep = t >= opts.TStart & t <= opts.TEnd;
    t    = t(keep);
    y    = y(keep, :);
end

% -------------------------------------------------------------------------

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
        error("loadRunData:FileNotFound", "No such file: %s", file);
    end
    Simulink.sdi.load(file);

    ids = Simulink.sdi.getAllRunIDs;
    if isempty(ids)
        error("loadRunData:NoRuns", "%s contains no runs.", file);
    end
    if runIndex < 0
        id = ids(end);
    elseif runIndex > numel(ids)
        error("loadRunData:BadRunIndex", ...
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
