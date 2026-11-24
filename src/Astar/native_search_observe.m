function [fit,code,result]=native_search_observe(varargin)
global MULTIRES_CALLS
alg=getenv('MULTIRES_ALG'); data=varargin{2}; isSearch=nargin==2;
diameters=case1_pipe_diameters();
[inputRequest,expectedDiameter]=case1_request_identity(data.S0,data.E0,diameters);
policy=case1_constraint_policy();
assert(data.BD==expectedDiameter && data.BD1==expectedDiameter && data.dsafe==policy.energy_dsafe && data.dmax==15,'Controls:Physical','Physical diameter, energy dsafe or field range changed.');
assert(data.energy_dsafe==policy.energy_dsafe && data.physical_collision_margin==policy.physical_collision_margin && strcmp(data.constraint_policy.id,policy.id),...
    'Controls:ConstraintPolicy','Energy dsafe and added physical collision margin must remain separate.');
assert(data.unit_BD==expectedDiameter,'Controls:Physical','Selected request unit_BD does not match its physical diameter.');
assert(varargin{1}.numAgent==6 && varargin{1}.maxIteration==0,'Controls:Budget','Candidate budget or optimization iterations changed.');
assert(isequal(data.dir,[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]),'Controls:Directions','Direction set changed.');
if strcmp(alg,'GFDHS')
    assert(data.planning_step_max==str2double(getenv('GFDHS_STEP_MAX')) && ismember(data.planning_step_max,[1 5 10 15 20 25]) && data.energy_range_device==15,'Controls:GFDHS','GFDHS sweep parameters changed.');
    assert(data.hard_clearance==policy.physical_collision_margin && data.step_tuning.guard_enabled && data.step_tuning.max_iterations==3000 && data.step_tuning.max_backtracks==128 && data.step_tuning.max_state_visits==8,'Controls:Guard','Physical collision margin, fixed guard or work budgets changed.');
else
    assert(data.unit==str2double(getenv('MULTIRES_H')),'Controls:Resolution','Grid size changed unexpectedly.');
    assert(data.interface_clearance==policy.interface_clearance,'Controls:ConstraintPolicy','Grid interface adds an unexpected physical gap.');
end
if isSearch && ~strcmp(alg,'GFDHS')
    expected=[data.interface_info.grid_start;data.interface_info.grid_goal]/data.unit;
    actual=[data.S;data.E];
    assert(isequal(actual,expected)||isequal(actual,flipud(expected)),'Interface:Overwritten','Search endpoints differ from repaired endpoints.');
    for endpoint=actual'
        idx=round(endpoint');assert(data.map_er(idx(1),idx(2),idx(3))~=88,'Interface:Blocked','Actual search endpoint is blocked.');
    end
end
if strcmp(alg,'GFDHS') && isSearch
    ports=[data.S0;data.E0];stubLength=2*expectedDiameter;
    [~,si]=min(sum(abs(ports-data.S),2));[~,ei]=min(sum(abs(ports-data.E),2));
    assert(si~=ei && abs(sum(abs(data.S-ports(si,:)))-stubLength)<1e-6 && abs(sum(abs(data.E-ports(ei,:)))-stubLength)<1e-6,'Repair:Ports','Cannot identify endpoint stubs');
    assert(norm(data.dir(data.dirSE(1),:)-(data.S-ports(si,:))/stubLength,Inf)<1e-6 && norm(data.dir(data.dirSE(2),:)-(ports(ei,:)-data.E)/stubLength,Inf)<1e-6,'Repair:Headings','Search endpoint headings must match oriented true-port stubs');
end
batch_search_event('start',alg,data,isSearch,numel(MULTIRES_CALLS)+1);
clock=tic;failure='';fit=Inf;code=[];result=struct();
try
    if strcmp(alg,'GFDHS')
        [fit,code,result]=aimFcn_1(varargin{:});
    elseif strcmp(alg,'GBH')
        [fit,code,result]=aimFcncell(varargin{:});
    else
        [fit,code,result]=aimFcn_astar(varargin{:});
    end
catch ex
    failure=getReport(ex,'extended','hyperlinks','off');
end
c=struct('is_search',isSearch,'start',data.S,'goal',data.E,...
    'diameter',data.BD,'diameter_BD1',data.BD1,'input_request_id',inputRequest,...
    'stub_length',2*expectedDiameter,'stub_length_rule','2*physical_diameter','dsafe',data.dsafe,...
    'energy_dsafe',data.dsafe,'physical_collision_margin',data.physical_collision_margin,'constraint_policy_id',policy.id,...
    'field_range',data.dmax,'grid_unit',data.unit,...
    'directions',data.dir,'code',code,'fitness',fit,'seconds',toc(clock),'failure',failure);
if strcmp(alg,'GFDHS')
    c.hard_clearance=data.hard_clearance;
    c.planning_step_max=data.planning_step_max;c.energy_range_device=data.energy_range_device;c.dmax_legacy=data.dmax;c.step_tuning=data.step_tuning;c.endpoint_headings=data.dirSE;c.endpoint_headings_verified=isSearch;
    global DDS_STEP_STATS
    c.dds_step_stats_cumulative=DDS_STEP_STATS;
    for stepMetric={'max_proposed_displacement','max_accepted_displacement'}
        if isfield(DDS_STEP_STATS,stepMetric{1})
            assert(DDS_STEP_STATS.(stepMetric{1})<=data.planning_step_max+1e-6,'StepSweep:Displacement','Recorded planning displacement exceeds configured maximum.');
        end
    end
end
if isfield(data,'S0'),c.physical_start=data.S0;end
if isfield(data,'E0'),c.physical_goal=data.E0;end
c.configuration=getenv('MULTIDIAM_CONFIG');c.pipe_diameters_input_order=diameters;
c.call_origin='new_parameter_setting';c.retained_prefix_count=0;
if ~strcmp(alg,'GFDHS'),c.interface_clearance=data.interface_clearance;end
if isfield(result,'failure_reason'),c.failure_reason=result.failure_reason;end
if isfield(result,'failure_point'),c.failure_point=result.failure_point;end
if isfield(result,'distance_to_goal'),c.distance_to_goal=result.distance_to_goal;end
if isfield(result,'searchedVoxels'),c.search_points_unique=size(unique(result.searchedVoxels,'rows'),1);else,c.search_points_unique=NaN;end
global CASE1_LAST_SEARCH_COUNTS NATIVE_ASTAR_COUNTS
if strcmp(alg,'Astar')
    c.native_astar_cumulative_counts=NATIVE_ASTAR_COUNTS;
    c.work_counters=[];c.work_counter_definition='Astar cumulative constructor/evaluation/trace rows; per-call unique trace recorded separately';
else
    c.work_counters=CASE1_LAST_SEARCH_COUNTS;
    c.work_counter_definition='start plus emitted neighbors; main-loop visits; raw proposals. Duplicates retained; stages differ by algorithm';
end
MULTIRES_CALLS{end+1}=c;
batch_search_event('end',alg,data,isSearch,numel(MULTIRES_CALLS));
if strcmp(alg,'GFDHS') && isfield(result,'failure_reason') && strcmp(result.failure_reason,'runtime_exception')
    error('StepSweep:NativeRuntimeException','Native objective swallowed a runtime exception; recorded diagnostics: %s',jsonencode(c.dds_step_stats_cumulative));
end
% Keep the failing native call visible; never replace seeds or silently retry.
if ~isempty(failure)
    if isSearch && strcmp(ex.identifier,'GBH:Exhausted')
        fit=Inf;code=[];result=struct('failure_reason','search_exhausted');
    else
        rethrow(ex);
    end
end
end
