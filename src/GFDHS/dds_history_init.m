function history=dds_history_init(point,limit,rolling)
% Membership-only storage; order is deliberately not exposed to callers.
history=struct('points',zeros(min(32,limit),3),'count',1, ...
    'cursor',1,'limit',limit,'rolling',rolling);
history.points(1,:)=point;
end
