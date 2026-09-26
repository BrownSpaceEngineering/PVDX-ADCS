function kept = pruneSDI(inFile, outFile, opts)
%PRUNESDI  Keep a few runs from a bloated SDI file and drop the rest.
%
%   An .mldatx is a snapshot of a whole SDI session, so a file grows every
%   time you save after more simulations.  Monte Carlo sweeps are the usual
%   culprit: every trial lands in SDI as its own run.  This loads the file
%   into a cleared session, deletes what you do not want, and writes a new
%   file.  The original is never modified.
%
%   pruneSDI(in, out)                          dry run, keeps nothing, lists all
%   pruneSDI(in, out, Name="285")              keep runs whose name contains "285"
%   pruneSDI(in, out, MinDuration=10000)       keep runs at least this long
%   pruneSDI(in, out, Name="285", DryRun=false)   actually write the file
%
%   NAME-VALUE
%     Name         Keep runs whose name contains this text.          ("")
%     MinDuration  Keep runs whose last sample is at least this, s.   (0)
%     Index        Keep these run indices outright.                   ([])
%     DryRun       Report only, write nothing.                     (true)
%
%   A run is kept if it matches ANY of the criteria given.
%
%   OUTPUT
%     kept   table of the runs that survived: index, name, duration
%
%   See also readmeStats, loadRunData.

    arguments
        inFile              {mustBeTextScalar}
        outFile             {mustBeTextScalar}
        opts.Name           {mustBeTextScalar} = ""
        opts.MinDuration    (1,1) double = 0
        opts.Index          (1,:) double = []
        opts.DryRun         (1,1) logical = true
    end

    if ~isfile(inFile)
        error("pruneSDI:NotFound", "No such file: %s", inFile);
    end
    d = dir(inFile);
    fprintf("\nLoading %s (%.2f GB) ...\n", inFile, d.bytes/2^30);

    Simulink.sdi.clear;                    % start from an empty session
    Simulink.sdi.load(inFile);
    ids = Simulink.sdi.getAllRunIDs;
    n   = numel(ids);
    fprintf("%d runs loaded.\n\n", n);

    idx  = (1:n)';
    name = strings(n,1);
    dur  = nan(n,1);
    keep = false(n,1);

    for k = 1:n
        r       = Simulink.sdi.getRun(ids(k));
        name(k) = string(r.Name);
        if r.SignalCount > 0
            tv = r.getSignalByIndex(1).Values.Time;
            if ~isempty(tv), dur(k) = tv(end); end
        end

        keep(k) = (strlength(opts.Name) > 0 && contains(name(k), opts.Name)) ...
               || (opts.MinDuration > 0 && dur(k) >= opts.MinDuration) ...
               || ismember(k, opts.Index);
    end

    kept = table(idx(keep), name(keep), dur(keep), ...
                 VariableNames=["index" "name" "duration_s"]);

    fprintf("  keeping %d of %d runs, dropping %d\n\n", ...
            sum(keep), n, sum(~keep));
    if isempty(kept)
        warning("pruneSDI:NothingKept", ...
                "No run matched. Nothing would be written.");
    else
        disp(kept);
    end

    if opts.DryRun
        fprintf("\nDry run - nothing written. Re-run with DryRun=false to save.\n\n");
        return
    end

    % ---- delete the rest and write the new file --------------------------
    drop = find(~keep);
    for k = drop(:)'
        Simulink.sdi.deleteRun(ids(k));
    end

    outDir = fileparts(outFile);
    if strlength(outDir) > 0 && ~isfolder(outDir)
        mkdir(outDir);
    end
    Simulink.sdi.save(outFile);

    dOut = dir(outFile);
    fprintf("\nWrote %s (%.2f GB, was %.2f GB)\n\n", ...
            outFile, dOut.bytes/2^30, d.bytes/2^30);
end
