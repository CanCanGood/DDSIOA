function fig=plot_case1(result_file)
%PLOT_CASE1 Display the obstacle boxes, physical ports and routed polylines.
root=fileparts(mfilename('fullpath'));
saved=load(result_file,'result');r=saved.result;
scene=jsondecode(fileread(fullfile(root,'data','case1_scene.json')));
fig=figure('Color','w');ax=axes(fig);hold(ax,'on');
faces=[1 2 3 4;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8];
for k=1:size(scene.obstacles,1)
    b=scene.obstacles(k,1:6)*scene.coordinate_unit_m;
    v=[b(1) b(2) b(3);b(4) b(2) b(3);b(4) b(5) b(3);b(1) b(5) b(3);...
       b(1) b(2) b(6);b(4) b(2) b(6);b(4) b(5) b(6);b(1) b(5) b(6)];
    patch(ax,'Vertices',v,'Faces',faces,'FaceColor',[.65 .7 .75],'FaceAlpha',.18,'EdgeColor',[.6 .6 .6],'HandleVisibility','off');
end
colors=lines(3);
for k=1:numel(r.routes)
    route=r.routes{k};p=route.points_m;id=route.input_request_id;
    plot3(ax,p(:,1),p(:,2),p(:,3),'-','LineWidth',2,'Color',colors(id,:),'DisplayName',sprintf('Request %d',id));
    scatter3(ax,p([1 end],1),p([1 end],2),p([1 end],3),35,colors(id,:),'filled','HandleVisibility','off');
end
domain=scene.workspace_native*scene.coordinate_unit_m;
xlim(ax,[0 domain(1)]);ylim(ax,[0 domain(2)]);zlim(ax,[0 domain(3)]);
axis(ax,'equal');grid(ax,'on');view(ax,3);xlabel(ax,'x (m)');ylabel(ax,'y (m)');zlabel(ax,'z (m)');
legend(ax,'Location','bestoutside');title(ax,sprintf('Case 1: %s, spatial setting %.2f m',r.config.method,r.config.resolution_m));
end
