function allowed=dds_tabu_escape_allowed(origin,point,flags,diameter,tolerance)
% Search-memory boxes are distinct from equipment. The caller must first
% certify the complete segment against equipment and domain constraints.
allowed=false;
if isempty(flags),allowed=true;return;end
if size(flags,2)<6||any(~isfinite(flags(:,1:6)),'all'),return;end
if numel(origin)~=3||numel(point)~=3||any(~isfinite([origin point])),return;end
axis=find(abs(point-origin)>tolerance);
if numel(axis)~=1,return;end
tol=min(tolerance,1e-8);
lo=min(flags(:,1:3),flags(:,4:6))-0.6*diameter;
hi=max(flags(:,1:3),flags(:,4:6))+0.6*diameter;
for k=1:size(lo,1)
    at_start=signed_clearance(origin,lo(k,:),hi(k,:));
    % For an axis-aligned segment, the minimum of max(lo-p,p-hi) occurs
    % where its moving coordinate is closest to the box midpoint.
    middle=origin;
    middle(axis)=min(max((lo(k,axis)+hi(k,axis))/2,min(origin(axis),point(axis))),max(origin(axis),point(axis)));
    deepest=signed_clearance(middle,lo(k,:),hi(k,:));
    if at_start>tol
        if deepest<=tol,return;end % No entry, including outside-to-outside crossings.
    elseif deepest<at_start-tol
        return; % No deeper penetration at any position along the move.
    end
end
allowed=true;
end

function value=signed_clearance(point,lo,hi)
value=max([lo-point,point-hi]);
end
