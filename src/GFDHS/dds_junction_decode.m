function [points,ok]=dds_junction_decode(start,encoding,data)
points=reshape(start,1,3); ok=false;
if ~isnumeric(encoding) || ~isreal(encoding) || mod(numel(encoding),2) || any(~isfinite(encoding)), return; end
for k=1:2:numel(encoding)
    direction=encoding(k); distance=encoding(k+1);
    if direction~=floor(direction) || direction<1 || direction>size(data.dir,1) || distance<0, return; end
    next=points(end,:)+data.dir(direction,:)*distance;
    if distance>0, points(end+1,:)=next; end %#ok<AGROW>
end
ok=true;
end
