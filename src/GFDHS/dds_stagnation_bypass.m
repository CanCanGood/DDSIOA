function [s,active,proposal,target]=dds_stagnation_bypass(s,p,goal,data,cfg)
% Opt-in main-constructor recovery. Constant-size state; no RNG or point log.
% A complete two-leg bypass must pass physical checks before activation.
active=false;proposal=[];target=goal;tol=cfg.state_tolerance;
if isempty(s)
    s=struct('count',0,'first',sum(abs(goal-p)),'best',sum(abs(goal-p)), ...
        'lo',p,'hi',p,'waypoints',zeros(0,3), ...
        'global_best',sum(abs(goal-p)),'no_progress',0,'last_recovery',0,'visits',0);
end
s.visits=s.visits+1;distance=sum(abs(goal-p));
if distance<s.global_best-cfg.state_tolerance
    s.global_best=distance;s.no_progress=0;
else
    s.no_progress=s.no_progress+1;
end
while ~isempty(s.waypoints) && all(abs(p-s.waypoints(1,:))<=tol)
    s.waypoints(1,:)=[];
    if isempty(s.waypoints)
        dds_step_stats('bypass_completed');s=reset_window(s,p,goal);
    end
end
if isempty(s.waypoints)
    s.count=s.count+1;s.best=min(s.best,sum(abs(goal-p)));
    s.lo=min(s.lo,p);s.hi=max(s.hi,p);
    if s.count>=96
        localStall=(s.first-s.best<cfg.base_planning_step) && max(s.hi-s.lo)<=180;
        globalStall=s.no_progress>=96 && s.visits-s.last_recovery>=96;
        stalled=localStall || globalStall;
        if globalStall,dds_step_stats('global_stagnation_triggers');end
        s=reset_window(s,p,goal);
        if stalled
            dds_step_stats('stagnation_triggers');s.last_recovery=s.visits;
            s.waypoints=make_bypass(p,goal,data);
            if isempty(s.waypoints),dds_step_stats('bypass_unavailable');
            else
                dds_step_stats('bypass_started');
                trace_event(p,goal,s.waypoints);
            end
        end
    end
end
if ~isempty(s.waypoints)
    target=s.waypoints(1,:);delta=target-p;
    if nnz(abs(delta)>tol)~=1
        s.waypoints=zeros(0,3);dds_step_stats('bypass_aborts');return;
    end
    proposal=p+sign(delta).*min(abs(delta),cfg.base_planning_step);
    active=true;dds_step_stats('bypass_steps');
end
end

function s=reset_window(s,p,goal)
s.count=0;s.first=sum(abs(goal-p));s.best=s.first;s.lo=p;s.hi=p;
end

function best=make_bypass(p,goal,data)
best=zeros(0,3);bestCost=Inf;
margin=max(data.BD/2+data.hard_clearance,1.2*data.BD1/2)+1e-3;
lo=data.dds_device_lo-margin;hi=data.dds_device_hi+margin;
[~,axes]=sort(abs(goal-p),'descend');
for axis=axes
    direction=sign(goal(axis)-p(axis));if direction==0,continue;end
    others=setdiff(1:3,axis);
    overlaps=all(p(others)>=lo(:,others) & p(others)<=hi(:,others),2);
    if direction>0
        distance=lo(:,axis)-p(axis);far=hi(:,axis);
    else
        distance=p(axis)-hi(:,axis);far=lo(:,axis);
    end
    ids=find(overlaps & distance>=0 & distance<abs(goal(axis)-p(axis)));
    [~,order]=sort(distance(ids));ids=ids(order);
    if isempty(ids),continue;end
    % Only the nearest blocking device on this goal-directed ray.
    id=ids(1);
    for sideAxis=others
        for side=[lo(id,sideAxis),hi(id,sideAxis)]
            a=p;a(sideAxis)=side;b=a;b(axis)=far(id);
            if ~clear_leg(p,a,data) || ~clear_leg(a,b,data),continue;end
            cost=sum(abs(a-p))+sum(abs(b-a))+sum(abs(goal-b));
            if cost<bestCost,best=[a;b];bestCost=cost;end
        end
    end
end
end

function ok=clear_leg(a,b,data)
ok=dds_clearance_segment(a,b,data);
if ~ok,return;end
u=potential_energy(a,b,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1);
ok=isscalar(u)&&isfinite(u);
if ~ok,return;end
u=checkpipe(a,b,data,data.dmax,data.BD1);
ok=isscalar(u)&&isfinite(u)&&u<=8;
end

function trace_event(p,goal,waypoints)
% Optional diagnostic, at most 16 events across all constructor retries.
global DDS_BYPASS_TRACE
if isstruct(DDS_BYPASS_TRACE) && isfield(DDS_BYPASS_TRACE,'events') && numel(DDS_BYPASS_TRACE.events)<16
    DDS_BYPASS_TRACE.events{end+1}=struct('point',p,'distance',sum(abs(goal-p)),'waypoints',waypoints);
end
end
