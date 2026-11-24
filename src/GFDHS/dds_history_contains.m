function found=dds_history_contains(history,point)
% Exact point membership, including legacy NaN-not-equal behavior.
found=any(all(history.points(1:history.count,:)==point,2));
end
