function [path,info]=dds_main_shortcut(path,data)
% Native-constructor simplification in direction/length space. Every new
% connector is checked continuously; the caller re-scores the final encoding.
% This is not applicable to branch graphs or internal crossover bridges.
orders=perms(1:3);tests=0;accepted=0;limit=3000;
[~,path,ok]=dds_junction_encode(path,data);assert(ok);
while tests<limit
    changed=false;n=size(path,1);
    for span=n-1:-1:2
        for first=1:n-span
            last=first+span;a=path(first,:);b=path(last,:);
            for order=1:6
                candidate=a;
                for axis=orders(order,:)
                    if candidate(end,axis)~=b(axis)
                        q=candidate(end,:);q(axis)=b(axis);candidate(end+1,:)=q; %#ok<AGROW>
                    end
                end
                if size(candidate,1)>=span+1,continue;end
                tests=tests+1;clear=true;
                for edge=1:size(candidate,1)-1
                    if ~dds_clearance_segment(candidate(edge,:),candidate(edge+1,:),data),clear=false;break;end
                    % Use the same native admissibility predicate as final scoring.
                    u=potential_energy(candidate(edge,:),candidate(edge+1,:),data.mapsize1,data.mapsize2,data,data.dmax,data.BD1);
                    if ~isscalar(u) || ~isreal(u) || ~isfinite(u),clear=false;break;end
                end
                if clear
                    proposed=[path(1:first-1,:);candidate;path(last+1:end,:)];
                    [~,proposed,encoded]=dds_junction_encode(proposed,data);
                    if encoded && size(proposed,1)<n
                        path=proposed;accepted=accepted+1;changed=true;break;
                    end
                end
                if tests>=limit,break;end
            end
            if changed||tests>=limit,break;end
        end
        if changed||tests>=limit,break;end
    end
    if ~changed,break;end
end
info=struct('connector_tests',tests,'accepted_shortcuts',accepted,'limit',limit);
end
