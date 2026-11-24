function r=full_port_route(legacy,trueStart,trueGoal,diameter)
stubLength=2*diameter;
assert(abs(sum(abs(legacy.points(1,:)-trueStart))-stubLength)<1e-7 && abs(sum(abs(legacy.points(end,:)-trueGoal))-stubLength)<1e-7,'Interface:StubLength','Native stub must equal 2*physical_diameter; never extend or move true ports.');
points=[trueStart;legacy.points;trueGoal];out=points(1,:);
for k=2:size(points,1)
    q=points(k,:);if all(abs(q-out(end,:))<1e-7),continue;end
    while size(out,1)>1
        a=out(end,:)-out(end-1,:);b=q-out(end,:);
        if max(abs(a/norm(a)-b/norm(b)))<1e-7,out(end,:)=[];else,break;end
    end
    out(end+1,:)=q;
end
assert(all(sum(abs(diff(out))>1e-7,2)==1),'Interface:Axis','Non-axis-aligned full route.');
r=struct('points',out,'length',sum(sum(abs(diff(out)),2)),'bends',size(out,1)-2,...
    'diameter',diameter,'stub_length',stubLength,'stub_length_rule','2*physical_diameter',...
    'true_start',trueStart,'true_goal',trueGoal,'offset_start',legacy.points(1,:),'offset_goal',legacy.points(end,:),...
    'endpoints_verified',isequal(out(1,:),trueStart)&&isequal(out(end,:),trueGoal));
end
