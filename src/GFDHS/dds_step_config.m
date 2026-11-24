function cfg = dds_step_config(data)
% Configure planning steps and search limits.
cfg = struct('guard_enabled',false, ...
    'max_iterations',3000,'max_backtracks',128,'max_state_visits',8, ...
    'min_planning_step',1,'recovery_hold',24,'state_tolerance',1e-6, ...
    'collision_backtracking_enabled',true,'collision_max_halvings',16,'tabu_escape_enabled',true);
if isfield(data,'step_tuning') && isstruct(data.step_tuning)
    names = fieldnames(cfg);
    for k = 1:numel(names)
        name = names{k};
        if isfield(data.step_tuning,name)
            value = data.step_tuning.(name);
            if (isnumeric(value) || islogical(value)) && isscalar(value) && isfinite(value)
                cfg.(name) = double(value);
            end
        end
    end
end
cfg.guard_enabled = logical(cfg.guard_enabled);
cfg.collision_backtracking_enabled = cfg.guard_enabled && logical(cfg.collision_backtracking_enabled);
cfg.tabu_escape_enabled = cfg.guard_enabled && logical(cfg.tabu_escape_enabled);
cfg.collision_max_halvings = max(0,floor(cfg.collision_max_halvings));
% data.dmax remains constant throughout an objective call, including recovery.
cfg.base_planning_step = max(1,data.dmax);
% Optional motion-only maximum: the potential field and goal-lock radius
% continue to use data.dmax. Invalid/absent values retain the v3 behavior.
cfg.planning_max_enabled = false;
cfg.planning_step_max = NaN;
if cfg.guard_enabled && isfield(data,'planning_step_max')
    candidateMax=data.planning_step_max;
    if isnumeric(candidateMax) && isscalar(candidateMax) && isreal(candidateMax) && isfinite(candidateMax) && candidateMax>0
        cfg.planning_max_enabled=true;
        cfg.planning_step_max=double(candidateMax);
        cfg.base_planning_step=cfg.planning_step_max;
    end
end
cfg.min_planning_step = min(cfg.base_planning_step,max(1e-6,cfg.min_planning_step));
cfg.max_iterations = max(1,floor(cfg.max_iterations));
cfg.max_backtracks = max(1,floor(cfg.max_backtracks));
cfg.max_state_visits = max(1,floor(cfg.max_state_visits));
cfg.recovery_hold = max(1,floor(cfg.recovery_hold));
cfg.state_tolerance = max(1e-10,cfg.state_tolerance);
end
