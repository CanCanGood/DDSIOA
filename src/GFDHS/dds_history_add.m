function history=dds_history_add(history,point)
% Preserve original capacity semantics, including main-history saturation.
if history.count<history.limit
    history.count=history.count+1;
    if history.count>size(history.points,1)
        history.points(min(history.limit,2*size(history.points,1)),3)=0;
    end
    history.points(history.count,:)=point;
elseif history.rolling
    history.points(history.cursor,:)=point;
    history.cursor=mod(history.cursor,history.limit)+1;
end
end
