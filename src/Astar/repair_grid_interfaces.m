function [data,info]=repair_grid_interfaces(data,pipes,diameter,physicalCollisionMargin)
% Keep the native obstacle grid intact. Attach exact physical endpoints to
% free lattice nodes beyond their OWN straight interface caps.
policy=case1_constraint_policy();
assert(physicalCollisionMargin==policy.interface_clearance && data.physical_collision_margin==physicalCollisionMargin && data.dsafe==policy.energy_dsafe,...
    'Interface:ConstraintPolicy','Interface gap is geometric; energy dsafe remains separate.');
assert(size(data.Start,1)==1,'Interface:MultipleStarts','Case1 expects one start per request.');
lo=double([pipes.mapx11(:),pipes.mapy11(:),pipes.mapz11(:)]);
hi=double([pipes.mapx21(:),pipes.mapy21(:),pipes.mapz21(:)]);
originalStart=double(data.Start);originalEnd=double(data.E);
[gs,sc,si]=attach(originalStart,double(data.S0),data,lo,hi,diameter,physicalCollisionMargin);
[ge,ec,ei]=attach(originalEnd,double(data.E0),data,lo,hi,diameter,physicalCollisionMargin);
info=struct('physical_start',originalStart,'physical_goal',originalEnd,...
    'grid_start',gs,'grid_goal',ge,'start_connector',sc,'end_connector',ec,...
    'start_own_stub',si,'end_own_stub',ei,'grid_unchanged',true,...
    'energy_dsafe',data.dsafe,'physical_collision_margin',physicalCollisionMargin,...
    'interface_clearance',physicalCollisionMargin,'constraint_policy_id',policy.id,...
    'connector_length',sum(sum(abs(diff(sc)),2))+sum(sum(abs(diff(ec)),2)));
data.Start=gs;data.E=ge;data.interface_info=info;
end
function [best,route,own]=attach(p,port,data,lo,hi,d,c)
t=1e-7;h=data.unit;delta=p-port;axis=find(abs(delta)>t);
assert(numel(axis)==1,'Interface:Direction','Interface must have one known outgoing axis.');
sgn=sign(delta(axis));
own=find(all(p>=lo-t & p<=hi+t,2)&all(port>=lo-t & port<=hi+t,2));
assert(numel(own)==1,'Interface:Ownership','Expected one unambiguous own interface stub.');
others=setdiff(1:size(lo,1),own);best=[];route=[];bestScore=[Inf Inf];
axesOther=setdiff(1:3,axis);
% Fixed geometric neighborhood, not repeated stochastic search. Prefer fewer
% connector bends, then smaller connector length, with stable tie breaking.
for layer=0:4
    qbase=round(p/h)*h;
    if sgn>0,qbase(axis)=ceil((hi(own,axis)+d/2+c)/h)*h+layer*h;
    else,qbase(axis)=floor((lo(own,axis)-d/2-c)/h)*h-layer*h;end
    for u=-1:1
        for v=-1:1
            q=qbase;q(axesOther)=q(axesOther)+[u v]*h;
            idx=round(q/h);
            if any(idx<1)||any(idx>size(data.map_er)),continue;end
            if data.map_er(idx(1),idx(2),idx(3))==88,continue;end
            for order=[axesOther;fliplr(axesOther)]'
                pts=p;now=p;now(axis)=q(axis);pts(end+1,:)=now;
                for a=order',now(a)=q(a);if any(abs(now-pts(end,:))>t),pts(end+1,:)=now;end;end
                ok=true;
                for k=2:size(pts,1)
                    a=pts(k-1,:);b=pts(k,:);lower=min(a,b);upper=max(a,b);
                    if any(lower<d/2+c-t)||any(upper>data.mapsize-d/2-c+t),ok=false;break;end
                    dl=[data.mapx1(:),data.mapy1(:),data.mapz1(:)];dh=[data.mapx2(:),data.mapy2(:),data.mapz2(:)];
                    gap=max(0,max(dl-upper,lower-dh));
                    if any(sum(gap.^2,2)<(d/2+c-t)^2),ok=false;break;end
                    gap=max(0,max(lo(others,:)-upper,lower-hi(others,:)));
                    if any(sum(gap.^2,2)<(d/2+c-t)^2),ok=false;break;end
                    % Only the axial first segment may continue its own stub.
                    if k>2
                        gap=max(0,max(lo(own,:)-upper,lower-hi(own,:)));
                        if sum(gap.^2)<(d/2+c-t)^2,ok=false;break;end
                    end
                end
                score=[size(pts,1)-2,sum(sum(abs(diff(pts)),2))];
                if ok && (score(1)<bestScore(1)||(score(1)==bestScore(1)&&score(2)<bestScore(2)))
                    best=q;route=pts;bestScore=score;
                end
            end
        end
    end
end
if isempty(best)
    record_grid_interface_failure(p,port,data,lo,hi,d,c,own,axis,sgn,axesOther);
end
assert(~isempty(best),'Interface:NoConnector','No valid lattice connector in fixed interface neighborhood.');
end
