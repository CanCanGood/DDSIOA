function [inputRequest,diameter]=case1_request_identity(start,goal,diameters)
% Exact immutable original endpoint pairs, independent of source sorting.
starts=[740 50 529;1055 582 529;415 136 45];
goals=[249 600 138;434 1336 638;249 631 23];
match=all(starts==start,2)&all(goals==goal,2);
reverse=all(starts==goal,2)&all(goals==start,2);
inputRequest=find(match|reverse);
assert(numel(inputRequest)==1,'Multidiameter:RequestIdentity','Physical endpoints do not identify exactly one original Case1 request.');
diameter=diameters(inputRequest);
end
