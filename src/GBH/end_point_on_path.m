function [node,xdir1]= end_point_on_path(datas,dirpath,data)
% 计算路径上距离给定点最近的点
% 输入：
% 路径起点
xdir =dirpath;    % 方向向量
tempstep=10;
xdir1=xdir;
curve=datas; 
for i=1:2:length(xdir)
    %% 计算路径中每两个连续点的欧氏距离，并将它加到总路径长度之上。
    %%在这里，norm(path(i,:)-path(i+1,:))表示的是点path(i,:)和点path(i+1,:)的欧氏距离。
    temp=xdir(i:i+1);
    % if i~=1&&temp(2)<2*data.BD1&&i~=length(dirpath)-1
    %     curve_len=curve_len+10000000;
    % end
    tempP1=curve;
    tempP2=tempP1+data.dir(temp(1),:)*tempstep;
    ij=tempstep;
    %% 累计节点能量值固定步长
    while ij<temp(2)
        if temp(2)<tempstep
            break;
        end
    
    [u1,~]=potential_energy(tempP1,tempP2, data.mapsize1,data.mapsize2, data, data.dmax,data.BD);
    [u2,~]= checkpipe(tempP1,tempP2, data, data.dmax,data.BD1);
if u1~=inf&&u2<9
  node=tempP1;
  return;
end
    tempP1=tempP2;
    tempP2=tempP1+data.dir(temp(1),:)*tempstep;
    ij=ij+tempstep;
    xdir1(2)=xdir1(2)- ij;
    end
    %% 判断是否发生弯曲并累计弯曲次数
    curve=curve+data.dir(temp(1),:)*temp(2);
    [u1,~]=potential_energy(tempP1,curve, data.mapsize1,data.mapsize2, data, data.dmax,data.BD);
    [u2,~]= checkpipe(tempP1,curve(end,:), data, data.dmax,data.BD1);
if u1~=inf&&u2<9
  node=tempP1;
  return;
end
        xdir1(1:2)=[];
end
 tempP1=curve;
 node=tempP1;
 xdir1=[dirpath(1),0];