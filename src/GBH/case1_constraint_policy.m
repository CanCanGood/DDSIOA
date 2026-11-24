function policy=case1_constraint_policy()
% Define energy and physical-clearance settings separately.
% Native diameter envelopes and finite-radius collision checks stay active.
policy=struct('id','case1_no_extra_clearance_v1',...
    'energy_dsafe',3,'physical_collision_margin',0,'interface_clearance',0,...
    'contact_allowed',true,'geometry_tolerance',1e-7);
end
