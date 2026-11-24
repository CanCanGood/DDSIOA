function [g,reason] = dds_step_guard(action,g,cfg,point,direction)
% Bounded state revisits, construction work and backtracking; no extra RNG.
reason = '';
switch action
    case 'init'
        g = struct('iterations',0,'backtracks',0, ...
            'local_step_cap',Inf,'recovery_moves',0, ...
            'visits',[]);
        if cfg.guard_enabled
            % Numeric canonical keys preserve the legacy quantized state.
            % Preserve signed-zero spelling in the original sprintf keys.
            % Column 7 holds the count; grow only when needed.
            g.visits=zeros(32,7);g.visit_count=0;
            g.key_cap_raw=NaN;g.key_cap=NaN;
        end
    case 'visit'
        g.iterations = g.iterations + 1;
        dds_step_stats('iterations');
        if ~cfg.guard_enabled, return; end
        if g.iterations > cfg.max_iterations
            reason = 'iteration_budget'; return;
        end
        if numel(point) ~= 3 || any(~isfinite(point)) || ~isscalar(direction) || ~isfinite(direction)
            reason = 'invalid_state'; return;
        end
        quant = round(point/cfg.state_tolerance);
        if ~isequal(g.key_cap_raw,g.local_step_cap)
            g.key_cap_raw=g.local_step_cap;
            g.key_cap=str2double(sprintf('%.12g',g.local_step_cap));
        end
        if direction~=fix(direction)
            direction=str2double(sprintf('%d',direction));
        end
        zero_sign=sum((quant(:)'==0 & 1./quant(:)'<0).*[1 2 4]);
        key=[quant(:)',direction,g.key_cap,zero_sign];
        slot=find(all(g.visits(1:g.visit_count,1:6)==key,2),1);
        if isempty(slot)
            g.visit_count=g.visit_count+1;slot=g.visit_count;
            if slot>size(g.visits,1)
                g.visits(2*size(g.visits,1),7)=0;
            end
            g.visits(slot,:)=[key,1];count=1;
        else
            count=g.visits(slot,7)+1;g.visits(slot,7)=count;
        end
        if count > 1, dds_step_stats('repeated_states'); end
        if count > cfg.max_state_visits
            [g,reason] = dds_step_guard('reduce',g,cfg);
            if ~isempty(reason), reason = 'repeated_state_budget'; end
        end
    case 'reduce'
        if isinf(g.local_step_cap)
            g.local_step_cap = max(cfg.min_planning_step,cfg.base_planning_step/2);
        elseif g.local_step_cap > cfg.min_planning_step
            g.local_step_cap = max(cfg.min_planning_step,g.local_step_cap/2);
        else
            reason = 'minimum_step'; return;
        end
        g.recovery_moves = 0;
        dds_step_stats('step_reductions');
    case 'backtrack'
        g.backtracks = g.backtracks + 1;
        dds_step_stats('backtracks');
        if cfg.guard_enabled && g.backtracks > cfg.max_backtracks
            reason = 'backtrack_budget';
        end
    case 'move'
        dds_step_stats('accepted_moves');
        if ~cfg.guard_enabled, return; end
        g.recovery_moves = g.recovery_moves + 1;
        if g.recovery_moves >= cfg.recovery_hold
            if isfinite(g.local_step_cap)
                g.local_step_cap = 2*g.local_step_cap;
                if g.local_step_cap>=cfg.base_planning_step, g.local_step_cap=Inf; end
                dds_step_stats('recoveries');
            end
            g.recovery_moves = 0;
        end
end
end
