function mininode= closest_point_on_path(data_S, dirpath,data)
% 计算路径上距离给定点最近的点
% 输入：
% 路径起点
node=dirpath{1}(1:3);
dis=h(data_S,node);
mindis=dis;
mininode=node;
for m=1:length(dirpath)
if length(dirpath{m})<2
    node=dirpath{m}(1);
end
node=dirpath{m}(1:3);

xdir =dirpath{m}(4:end);    % 方向向量
detaP1=data_S-node;
detaP=[];
for count=1:3
vi=[0,0,0];
vi(count)=sign(detaP1(count));
v=find(ismember(data.dir, vi, 'rows'));
if isempty(v)
    v=0;
end
detaP=[detaP,v,abs(detaP1(count))];
end
for i=1:2:length(xdir)
  if xdir(i)==detaP(fix((xdir(i)+1)/2)*2-1)
     if xdir(i+1)<=detaP(fix((xdir(i)+1)/2)*2)
         dis=dis-xdir(i+1);
         node=node+data.dir( xdir(i),:)*xdir(i+1);
         if dis<mindis
             mindis=dis;
             mininode=node;
         end
         detaP(fix((xdir(i)+1)/2)*2)=detaP(fix((xdir(i)+1)/2)*2)-xdir(i+1);
         dis=h(data_S,node);
     else
         dis=dis-detaP(fix((xdir(i)+1)/2)*2);
         if dis<mindis
             mindis=dis;
             mininode=node+data.dir(xdir(i),:)*detaP(fix((xdir(i)+1)/2)*2);
         end
         node=node+data.dir( xdir(i),:)*xdir(i+1);
         detaP(fix((xdir(i)+1)/2)*2)=abs(detaP(fix((xdir(i)+1)/2)*2)-xdir(i+1));
         detaP(fix((xdir(i)+1)/2)*2-1)=xdir(i)+mod(xdir(i),2)*2-1;
         dis=h(data_S,node);
     end



  else
       node=node+data.dir( xdir(i),:)*xdir(i+1);
       detaP(fix((xdir(i)+1)/2)*2)=detaP(fix((xdir(i)+1)/2)*2)+xdir(i+1);
       detaP(fix((xdir(i)+1)/2)*2-1)=xdir(i)+mod(xdir(i),2)*2-1;
       dis=h(data_S,node);
  end



end
end