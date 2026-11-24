function [encoding,points,ok]=dds_junction_encode(points,data)
% Canonical direction-length encoding; retain the native initial direction.
tol=1e-7;
if isfield(data,'junction_tolerance'), tol=data.junction_tolerance; end
encoding=[]; ok=false;
if isempty(points) || size(points,2)~=3 || any(~isfinite(points(:))), return; end
clean=points(1,:);
for i=2:size(points,1)
    delta=points(i,:)-clean(end,:); axis=find(abs(delta)>tol);
    if isempty(axis), continue; end
    if numel(axis)~=1, return; end
    next=clean(end,:); next(axis)=points(i,axis);
    unit=zeros(1,3); unit(axis)=sign(delta(axis));
    direction=find(all(data.dir==unit,2),1);
    if isempty(direction), return; end
    distance=abs(delta(axis));
    if ~isempty(encoding) && encoding(end-1)==direction
        encoding(end)=encoding(end)+distance;
        clean(end,:)=next;
    else
        encoding=[encoding,direction,distance]; %#ok<AGROW>
        clean(end+1,:)=next; %#ok<AGROW>
    end
end
if isempty(encoding)
    encoding=[data.dirSE(1),0];
elseif floor((encoding(1)+1)/2)~=floor((data.dirSE(1)+1)/2)
    encoding=[data.dirSE(1),0,encoding];
end
points=clean; ok=true;
end
