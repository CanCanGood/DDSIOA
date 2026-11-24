function data=dds_prepare_geometry(data)
% Objective-local immutable geometry: rebuild at each call, never global.
% Prepare before the potential function appends temporary tabu obstacles.
if ~isfield(data,'hard_clearance'),return;end
data.dds_device_lo=[data.mapx1(:),data.mapy1(:),data.mapz1(:)];
data.dds_device_hi=[data.mapx2(:),data.mapy2(:),data.mapz2(:)];
% Preserve the exact multiply/divide order and decimal convention of legacy
% energy evaluation. These are separate from the actual-clearance geometry.
data.dds_energy_diameter=data.BD1;
data.dds_energy_boxes=round([data.dds_device_lo-1.2*data.BD1/2, ...
    data.dds_device_hi+1.2*data.BD1/2],10);
data.dds_pipe_boxes=round([data.mapx11(:)-1.1*data.BD1/2, ...
    data.mapy11(:)-1.1*data.BD1/2,data.mapz11(:)-1.1*data.BD1/2, ...
    data.mapx21(:)+1.1*data.BD1/2,data.mapy21(:)+1.1*data.BD1/2, ...
    data.mapz21(:)+1.1*data.BD1/2],10);
if isfield(data,'hard_external') && ~isempty(data.hard_external)
    data.dds_external_lo=data.hard_external(:,1:3);
    data.dds_external_hi=data.hard_external(:,4:6);
    data.dds_external_threshold2=(data.BD/2+data.hard_external(:,7)/2+data.hard_clearance-1e-7).^2;
end
end
