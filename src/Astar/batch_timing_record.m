function batch_timing_record(name,startNanoseconds)
% Separate monotonic clock: do not reset the source's unnamed tic/toc timer.
global BATCH_TIMINGS
elapsed=(double(java.lang.System.nanoTime)-startNanoseconds)/1e9;
BATCH_TIMINGS{end+1}=struct('name',name,'seconds',elapsed,...
    'configuration',getenv('MULTIDIAM_CONFIG'));
end
