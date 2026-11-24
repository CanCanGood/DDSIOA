function [point,accepted,info] = dds_step_collision_retry(origin,proposal,flagrec,data,cfg)
% Test each trial from the same origin. Never commit a colliding position.
% No RNG or change to data.dmax. Equipment and domain are checked first.
point = proposal;
accepted = false;
info = struct('halvings',0,'reason','','initial_length',NaN,'accepted_length',NaN);
delta = proposal-origin;
if numel(origin)~=3 || numel(proposal)~=3 || ~isreal(delta) || ...
        any(~isfinite(delta)) || nnz(abs(delta)>cfg.state_tolerance)~=1
    info.reason='invalid_candidate';
    dds_step_stats('invalid_candidates');
    return;
end
info.initial_length=sum(abs(delta));
dds_step_stats('collision_tested_candidates');
everCollided=false;
while true
    point=origin+delta;
    if sum(abs(delta))<=cfg.state_tolerance
        info.reason='zero_step';
        dds_step_stats('zero_moves_rejected');
        return;
    end
    dds_step_stats('segment_checks');
    energy=potential_energy(origin,point,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1);
    physicalCollision=isscalar(energy) && isreal(energy) && energy==Inf;
    if physicalCollision
        info.reason='full_segment_collision';
        dds_step_stats('full_segment_collision');
        dds_step_stats('segment_rejections');
    elseif ~isscalar(energy) || ~isreal(energy) || ~isfinite(energy)
        info.reason='invalid_segment_energy';
        dds_step_stats('segment_rejections');
        dds_step_stats('invalid_candidate_energy');
        return;
    else
        energy=potential_energy(point,point,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1,flagrec);
        tabuCollision=isscalar(energy) && isreal(energy) && energy==Inf;
        if cfg.tabu_escape_enabled && ~isempty(flagrec)
            % A temporary search box may touch/enclose the current position.
            % Permit physical-clear tangent/escape progress, never entry or
            % movement deeper into that box, even between the endpoints.
            allowed=dds_tabu_escape_allowed(origin,point,flagrec,data.BD1,cfg.state_tolerance);
            if allowed && tabuCollision
                tabuCollision=false;energy=0;
                dds_step_stats('tabu_tangent_or_escape_moves');
            elseif ~allowed
                tabuCollision=true;energy=Inf;
                dds_step_stats('tabu_segment_rejections');
            end
        end
        if tabuCollision
            info.reason='tabu_endpoint_collision';
            dds_step_stats('tabu_endpoint_collision');
            dds_step_stats('tabu_rejections');
        elseif ~isscalar(energy) || ~isreal(energy) || ~isfinite(energy)
            info.reason='invalid_endpoint_energy';
            dds_step_stats('tabu_rejections');
            dds_step_stats('invalid_candidate_energy');
            return;
        else
            accepted=true;
            info.reason='accepted';
            info.accepted_length=sum(abs(delta));
            if info.halvings>0
                dds_step_stats('collision_recovered_candidates');
                dds_step_stats('collision_halving_successes');
            else
                dds_step_stats('collision_clear_candidates');
            end
            return;
        end
    end
    if ~everCollided
        dds_step_stats('collision_candidates');
        everCollided=true;
    end
    % A trial below the minimum is not created. A native sub-minimum proposal
    % may still be accepted above if clear, but may not be halved after collision.
    if info.halvings>=cfg.collision_max_halvings || sum(abs(delta))/2<cfg.min_planning_step
        dds_step_stats('collision_exhausted_candidates');
        return;
    end
    delta=delta/2;
    info.halvings=info.halvings+1;
    dds_step_stats('collision_halvings');
end
end
