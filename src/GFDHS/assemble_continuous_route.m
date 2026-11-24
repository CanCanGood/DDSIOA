function r=assemble_continuous_route(code,data)
% GFDHS codes already contain model-unit lengths; no grid multiplier.
p=data.Start;points=p;
for k=1:2:numel(code)
    p=p+data.dir(code(k),:)*code(k+1);points(end+1,:)=p;
end
assert(max(abs(points(end,:)-data.E))<1e-7,'Interface:Endpoint','Continuous selected route endpoint mismatch.');
assert(all(sum(abs(diff(points))>1e-7,2)<=1),'Interface:Axis','Non axial route');
r=struct('points',points,'length',sum(sum(abs(diff(points)),2)),'bends',size(points,1)-2,'endpoints_verified',true);
end
