function dds_step_stats(name,value)
% Accumulate all calls in the current run. Runner resets global to struct().
global DDS_STEP_STATS
% Common unit counters: keep exact totals while bypassing diagnostic dispatch.
% All max_* names retain the original maximum semantics, including arity one.
if nargin==1 && ischar(name) && ~isempty(name) && name(1)~='m' && isfield(DDS_STEP_STATS,name)
    DDS_STEP_STATS.(name)=DDS_STEP_STATS.(name)+1;
    return;
end
if ~isstruct(DDS_STEP_STATS), DDS_STEP_STATS = struct(); end
if nargin < 2, value = 1; end
if strcmp(name,'max_proposed_displacement') || strcmp(name,'max_accepted_displacement')
    if isscalar(value) && isreal(value) && isfinite(value) && value>=0
        if ~isfield(DDS_STEP_STATS,name), DDS_STEP_STATS.(name)=0; end
        DDS_STEP_STATS.(name)=max(DDS_STEP_STATS.(name),value);
    end
    return;
end
if strcmp(name,'reject')
    if ~isfield(DDS_STEP_STATS,'guard_rejections'), DDS_STEP_STATS.guard_rejections = struct(); end
    key = matlab.lang.makeValidName(value);
    if ~isfield(DDS_STEP_STATS.guard_rejections,key), DDS_STEP_STATS.guard_rejections.(key) = 0; end
    DDS_STEP_STATS.guard_rejections.(key) = DDS_STEP_STATS.guard_rejections.(key) + 1;
    name = 'guard_failures'; value = 1;
elseif strcmp(name,'exception')
    DDS_STEP_STATS.last_exception_identifier = value.identifier;
    DDS_STEP_STATS.last_exception_message = value.message;
    name = 'caught_exceptions'; value = 1;
end
if ~isfield(DDS_STEP_STATS,name), DDS_STEP_STATS.(name) = 0; end
DDS_STEP_STATS.(name) = DDS_STEP_STATS.(name) + value;
end
