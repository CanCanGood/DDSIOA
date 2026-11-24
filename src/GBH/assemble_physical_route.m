function r=assemble_physical_route(code,data)
info=data.interface_info;p=info.grid_start;points=p;
for k=1:2:numel(code)
    p=p+data.dir(code(k),:)*code(k+1)*data.unit;points(end+1,:)=p;
end
assert(max(abs(points(end,:)-info.grid_goal))<1e-7,'Interface:Endpoint','Selected route endpoint mismatch.');
points=[info.start_connector;points(2:end,:);flipud(info.end_connector(1:end-1,:))];
out=points(1,:);
for k=2:size(points,1)
    q=points(k,:);if all(abs(q-out(end,:))<1e-7),continue;end
    while size(out,1)>1
        a=out(end,:)-out(end-1,:);b=q-out(end,:);
        if max(abs(a/norm(a)-b/norm(b)))<1e-7,out(end,:)=[];else,break;end
    end
    out(end+1,:)=q;
end
assert(all(sum(abs(diff(out))>1e-7,2)==1),'Interface:Axis','Non-axis-aligned geometry.');
r=struct('points',out,'length',sum(sum(abs(diff(out)),2)),'bends',size(out,1)-2,...
    'endpoints_verified',max(abs(out(1,:)-info.physical_start))<1e-7 && max(abs(out(end,:)-info.physical_goal))<1e-7);
end
