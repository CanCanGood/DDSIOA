function [pruned,report]=dds_main_loop_prune(path,tolerance)
% Erase closed walks from one physical orthogonal main-route centerline.
% Never use this function on a true branch graph or on a direction encoding.
% Intersections are exact: tolerance NEVER joins merely nearby centerlines.
if nargin<2,tolerance=1e-7;end
assert(isnumeric(path)&&isreal(path)&&ismatrix(path)&&size(path,2)==3&& ...
    size(path,1)>=1&&all(isfinite(path(:))), ...
    'DDSLoopPrune:InvalidPath','Expected a nonempty finite real Nx3 centerline.');
assert(isnumeric(tolerance)&&isreal(tolerance)&&isscalar(tolerance)&& ...
    isfinite(tolerance)&&tolerance>0&&tolerance<=1e-7, ...
    'DDSLoopPrune:InvalidTolerance','Require 0 < tolerance <= 1e-7.');
path=double(path);
assert(all(sum(diff(path,1,1)~=0,2)<=1),'DDSLoopPrune:NonOrthogonal', ...
    'Segments must be exactly axis aligned; small off-axis components are not silently projected.');
original_length=route_length(path);pruned=simplify_trace(path);
initial_vertex_count=size(pruned,1);witnesses=cell(1,max(0,initial_vertex_count-1));
loop_count=0;intersection_tests=0;
while size(pruned,1)>1
    found=false;n=size(pruned,1);
    % Earliest closing segment first; then earliest retained segment.
    for later=2:n-1
        for earlier=1:later-1
            intersection_tests=intersection_tests+1;
            [hit,point,kind]=intersection(pruned(earlier,:),pruned(earlier+1,:), ...
                pruned(later,:),pruned(later+1,:));
            if ~hit,continue;end
            closed_length=sum(abs(point-pruned(earlier+1,:)))+ ...
                route_length(pruned(earlier+1:later,:))+sum(abs(pruned(later,:)-point));
            if closed_length<=0,continue;end % Ordinary adjacent elbow only.
            before=pruned;before_length=route_length(before);
            proposed=simplify_trace([before(1:earlier,:);point;before(later+1:end,:)]);
            after_length=route_length(proposed);
            assert(size(proposed,1)<n,'DDSLoopPrune:NoProgress', ...
                'Every erased positive closed walk must remove at least one vertex.');
            assert(isequal(proposed(1,:),path(1,:))&&isequal(proposed(end,:),path(end,:)), ...
                'DDSLoopPrune:ChangedEndpoint','Loop erasure changed a physical endpoint.');
            slack=tolerance+32*eps(max(1,before_length))*n;
            assert(after_length<=before_length+slack,'DDSLoopPrune:IncreasedLength', ...
                'Loop erasure must not increase centerline length.');
            loop_count=loop_count+1;
            witnesses{loop_count}=struct('iteration',loop_count, ...
                'earlier_segment',earlier,'later_segment',later, ...
                'segment_index_scope','current simplified path before this erasure', ...
                'earlier_endpoints',before(earlier:earlier+1,:), ...
                'later_endpoints',before(later:later+1,:),'intersection',point, ...
                'intersection_kind',kind,'closed_walk_length',closed_length, ...
                'length_before',before_length,'length_after',after_length, ...
                'vertices_before',n,'vertices_after',size(proposed,1));
            pruned=proposed;found=true;break;
        end
        if found,break;end
    end
    if ~found,break;end
    assert(loop_count<initial_vertex_count,'DDSLoopPrune:IterationBound', ...
        'Strict vertex decrease did not bound the loop-erasure iterations.');
end
if loop_count>0,witnesses=[witnesses{1:loop_count}];else,witnesses=struct([]);end
final_length=route_length(pruned);
report=struct('method','deterministic_exact_orthogonal_closed_walk_erasure', ...
    'scope','One independent main-route physical polyline; not a true branch graph.', ...
    'tolerance',tolerance,'intersection_uses_tolerance',false, ...
    'input_vertices',size(path,1),'initial_simplified_vertices',initial_vertex_count, ...
    'output_vertices',size(pruned,1),'input_length',original_length, ...
    'output_length',final_length,'removed_length',original_length-final_length, ...
    'loop_count',loop_count,'intersection_tests',intersection_tests, ...
    'maximum_erasure_iterations',max(0,initial_vertex_count-1), ...
    'geometry_changed',loop_count>0,'representation_changed',~isequal(path,pruned), ...
    'exact_endpoints_preserved',isequal(path(1,:),pruned(1,:))&&isequal(path(end,:),pruned(end,:)), ...
    'retained_geometry','Every output segment follows the original centerline trace; no offset, snap, or new shortcut.', ...
    'limitations','Does not certify equipment clearance, finite-radius self clearance, port tangency, other pipes, or global optimality.', ...
    'witnesses',witnesses);
end

function value=route_length(path)
value=sum(sum(abs(diff(path,1,1))));
end

function path=simplify_trace(points)
% Remove exactly duplicate vertices and merge only forward collinear edges.
path=zeros(size(points));count=0;
for k=1:size(points,1)
    point=points(k,:);
    if count>0&&isequal(path(count,:),point),continue;end
    count=count+1;path(count,:)=point;
    while count>=3
        first=path(count-1,:)-path(count-2,:);second=path(count,:)-path(count-1,:);
        if isequal(find(first~=0),find(second~=0))&&dot(first,second)>0
            path(count-1,:)=path(count,:);count=count-1;
        else,break;end
    end
end
path=path(1:count,:);
end

function [hit,point,kind]=intersection(a,b,c,d)
% Both input segments are positive and exactly axis aligned.
axis1=find(a~=b);axis2=find(c~=d);point=zeros(1,3);kind='';hit=false;
assert(isscalar(axis1)&&isscalar(axis2),'DDSLoopPrune:InternalSegment', ...
    'Unexpected degenerate segment after simplification.');
if axis1==axis2
    fixed=(1:3)~=axis1; % Exact complementary axes; no general set operation.
    if ~isequal(a(fixed),c(fixed)),return;end
    lower=max(min(a(axis1),b(axis1)),min(c(axis1),d(axis1)));
    upper=min(max(a(axis1),b(axis1)),max(c(axis1),d(axis1)));
    if lower>upper,return;end
    point=a;
    % The first point of overlap along the earlier segment erases backtracks.
    if b(axis1)>a(axis1),point(axis1)=lower;else,point(axis1)=upper;end
    hit=true;
    if lower<upper,kind='collinear_overlap';else,kind='collinear_point';end
else
    fixed=6-axis1-axis2; % Distinct axes in {1,2,3}: the remaining axis.
    if ~isequal(a(fixed),c(fixed)),return;end
    point=a;point(axis1)=c(axis1);
    if point(axis1)<min(a(axis1),b(axis1))||point(axis1)>max(a(axis1),b(axis1))|| ...
            point(axis2)<min(c(axis2),d(axis2))||point(axis2)>max(c(axis2),d(axis2)),return;end
    hit=true;kind='perpendicular_point';
end
end
