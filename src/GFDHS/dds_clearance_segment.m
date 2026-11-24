function ok=dds_clearance_segment(a,b,data)
% Exact Euclidean distances for axis-aligned centerline segments and AABBs.
% Parent contact is handled by the dedicated true-junction validator.
ok=true;
if ~isfield(data,'hard_clearance'),return;end
c=data.hard_clearance;r=data.BD/2;t=1e-7;
lo=min(a,b);hi=max(a,b);
if any(lo<r+c-t)||any(hi>data.mapsize-r-c+t),ok=false;return;end
if isfield(data,'dds_device_lo')
    bl=data.dds_device_lo;bh=data.dds_device_hi;
else
    bl=[data.mapx1(:),data.mapy1(:),data.mapz1(:)];
    bh=[data.mapx2(:),data.mapy2(:),data.mapz2(:)];
end
gap=max(0,max(bl-hi,lo-bh));
if any(sum(gap.^2,2)<(r+c-t)^2),ok=false;return;end
if isfield(data,'hard_external') && ~isempty(data.hard_external)
    if isfield(data,'dds_external_lo')
        gap=max(0,max(data.dds_external_lo-hi,lo-data.dds_external_hi));
        threshold2=data.dds_external_threshold2;
    else
        e=data.hard_external;
        gap=max(0,max(e(:,1:3)-hi,lo-e(:,4:6)));
        threshold2=(r+e(:,7)/2+c-t).^2;
    end
    if any(sum(gap.^2,2)<threshold2),ok=false;end
end
end
