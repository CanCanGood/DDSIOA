function [goalDelta,branchDistance]=dds_candidate_distances(origin,candidates,goal,targets)
% Batch row-wise Manhattan distances without changing terminal summation order.
current=h(origin,goal);
if size(goal,1)==1
    goalDelta=sum(abs(candidates-goal),2)-current;
else
    goalDelta=zeros(size(candidates,1),1);
    for k=1:size(candidates,1),goalDelta(k)=h(candidates(k,:),goal)-current;end
end
branchDistance=zeros(size(candidates,1),1);
for k=1:size(targets,1)
    branchDistance=branchDistance+sum(abs(candidates-targets(k,:)),2);
end
end
