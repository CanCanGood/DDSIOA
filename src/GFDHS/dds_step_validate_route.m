function [reason,dirpath] = dds_step_validate_route(dirpath,sP,endP,data,cfg)
% Verify exactly the route represented by the returned direction encoding.
reason = '';
if numel(dirpath)<5 || mod(numel(dirpath)-3,2)~=0 || any(~isfinite(dirpath))
    reason = 'invalid_route_encoding'; return;
end
if isfield(data,'hard_clearance') && isfield(data,'dds_prune_constructor') && data.dds_prune_constructor && ~dds_junction_enabled(data)
    [raw,decoded]=dds_junction_decode(sP,dirpath(4:end),data);
    if ~decoded,reason='invalid_route_encoding';return;end
    [raw,~]=dds_main_loop_prune(raw,1e-7);
    [raw,~]=dds_main_shortcut(raw,data);
    [encoding,~,encoded]=dds_junction_encode(raw,data);
    if ~encoded,reason='invalid_route_encoding';return;end
    dirpath=[sP,encoding];
end
point = sP;
for k=4:2:numel(dirpath)
    direction = dirpath(k); distance = dirpath(k+1);
    if direction~=floor(direction) || direction<1 || direction>size(data.dir,1) || distance<0
        reason = 'invalid_route_encoding'; return;
    end
    next = point + data.dir(direction,:)*distance;
    if distance>cfg.state_tolerance
        dds_step_stats('final_segment_checks');
        energy = potential_energy(point,next,data.mapsize1,data.mapsize2,data,data.dmax,data.BD1);
        if ~isscalar(energy) || ~isfinite(energy)
            reason = 'invalid_route_segment'; return;
        end
    end
    point = next;
end
if any(abs(point-endP)>cfg.state_tolerance)
    reason = 'incomplete_route';
end
end
