function [nextP,g] = dds_step_candidates(baseP,nowP,endP,nearGoal,jumpFlag,flagrec,data,cfg,g)
% Check complete native proposals, with smaller local retries on failure.
% data.dmax and the energy sampling interval are never changed here.
dds_step_stats('candidate_batches');
dds_step_stats('raw_candidates',size(baseP,1));
if ~cfg.guard_enabled
    nextP = baseP;
    return;
end
clearOrigin=false;
if cfg.planning_max_enabled && ~nearGoal && ~isempty(baseP)
    % Expansion is allowed only when both original point checks are exactly
    % zero. Segment checking below is still mandatory even for free endpoints.
    physicalEnergy=potential_energy(nowP,nowP,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1);
    tabuEnergy=potential_energy(nowP,nowP,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1,flagrec);
    dds_step_stats('planning_origin_checks',2);
    clearOrigin=isscalar(physicalEnergy) && isreal(physicalEnergy) && physicalEnergy==0 && ...
        isscalar(tabuEnergy) && isreal(tabuEnergy) && tabuEnergy==0;
    if clearOrigin
        dds_step_stats('clear_space_batches');
    else
        dds_step_stats('nonzero_origin_batches');
    end
end
while true
    nextP = zeros(0,3);
    for k = 1:size(baseP,1)
        p = baseP(k,:);
        delta = p-nowP;
        if any(~isfinite(delta)) || nnz(abs(delta)>cfg.state_tolerance) ~= 1
            dds_step_stats('invalid_candidates'); continue;
        end
        if cfg.planning_max_enabled
            if ~isreal(delta) || nnz(delta~=0)~=1
                dds_step_stats('invalid_candidates'); continue;
            end
            nativeLength=sum(abs(delta));
            motionAxis=find(delta~=0,1);
            if ~nearGoal && clearOrigin
                proposedLength=cfg.planning_step_max;
                if proposedLength>nativeLength, dds_step_stats('clear_space_expansions'); end
            else
                proposedLength=min(nativeLength,cfg.planning_step_max);
                if nearGoal, dds_step_stats('near_goal_native_candidates'); end
            end
            if proposedLength<nativeLength, dds_step_stats('planning_max_clamps'); end
            % The direction is exactly axis-aligned; assigning its magnitude
            % avoids rescaling huge finite native proposals by tiny ratios.
            delta(motionAxis)=sign(delta(motionAxis))*proposedLength;
            remaining=endP(motionAxis)-nowP(motionAxis);
            if remaining*delta(motionAxis)>0 && abs(delta(motionAxis))>abs(remaining)
                delta(motionAxis)=remaining;
                dds_step_stats('goal_axis_clamps');
            end
        else
            % Preserve v3's native candidate/counter behavior when no field is set.
            dds_step_stats('max_proposed_displacement',sum(abs(delta)));
        end
        if isfinite(g.local_step_cap)
            delta = delta*min(1,g.local_step_cap/sum(abs(delta)));
        end
        if cfg.planning_max_enabled
            dds_step_stats('max_proposed_displacement',sum(abs(delta)));
        end
        p = nowP+delta;
        axis = find(abs(delta)>cfg.state_tolerance,1);
        margin = max(1e-8,cfg.state_tolerance);
        clipped = min(max(p(axis),data.mapsize2(axis)+margin),data.mapsize1(axis)-margin);
        if clipped ~= p(axis), dds_step_stats('boundary_clips'); end
        p(axis) = clipped;
        if contains_array(endP,p,nowP), p = endP; end
        if sum(abs(p-nowP)) <= cfg.state_tolerance
            dds_step_stats('zero_moves_rejected'); continue;
        end
        if cfg.collision_backtracking_enabled
            % Resolve this direction now, even when another direction is clear.
            % Near-goal projections use the same collision-driven halving rule.
            [p,accepted,~]=dds_step_collision_retry(nowP,p,flagrec,data,cfg);
            if accepted, nextP(end+1,:)=p; end %#ok<AGROW>
            continue;
        end
        dds_step_stats('segment_checks');
        energy = potential_energy(nowP,p,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1);
        if ~isscalar(energy) || ~isfinite(energy)
            dds_step_stats('segment_rejections'); continue;
        end
        energy = potential_energy(p,p,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1,flagrec);
        if ~isscalar(energy) || ~isfinite(energy)
            dds_step_stats('tabu_rejections'); continue;
        end
        nextP(end+1,:) = p; %#ok<AGROW>
    end
    if ~isempty(nextP), return; end
    dds_step_stats('empty_candidates');
    % Every individual ray already exhausted its own bounded smaller trials.
    % Re-running all rays under a smaller shared cap would duplicate that work.
    if cfg.collision_backtracking_enabled, return; end
    if nearGoal, return; end
    [g,reason] = dds_step_guard('reduce',g,cfg);
    if ~isempty(reason), return; end
end
end
