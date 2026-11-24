function [U,U2]= potential_energy(point1,point2, map_bounds1,map_bounds2, data, dmax,BD,path)
% 输入参数：
%   point: 点的坐标 [x, y, z]
%   map_bounds: 地图边界 [a, b, c]（本函数未使用该参数）
%   obstacles: N×6矩阵，每行表示一个障碍物的[xmin, ymin, zmin, xmax, ymax, zmax]
%   dmax: 障碍物影响的最大距离
% 输出：
%   U: 障碍物对点的总势能（点不在障碍物内且距离超过dmax时为0）
if isfield(data,'energy_range_device')
    dmax=data.energy_range_device; % Explicit sensitivity trial; goal range unchanged.
end
if isfield(data,'hard_clearance') && ~dds_clearance_segment(point1,point2,data)
    U=Inf;U2=Inf;return;
end
decimal_places = 10; % 保留4位小数
if nargin==8&&size(path,2)>=6
data.mapx1=[data.mapx1;min(path(:,1),path(:,4))];
data.mapx2=[data.mapx2;max(path(:,1),path(:,4))];
data.mapy1=[data.mapy1;min(path(:,2),path(:,5))];
data.mapy2=[data.mapy2;max(path(:,2),path(:,5))];
data.mapz1=[data.mapz1;min(path(:,3),path(:,6))];
data.mapz2=[data.mapz2;max(path(:,3),path(:,6))];
data.obstacles=[data.obstacles;path(:,end-4:end)];
end

x = point2(1);
y = point2(2);
z = point2(3);
index1=find(point1~=point2);
if point1(index1)>point2(index1)
    point2=point1+point2;
    point1=point2-point1;
    point2=point2-point1;
end

x1 = point1(1)-dmax-1.2*BD/2;
y1 = point1(2)-dmax-1.2*BD/2;
z1 = point1(3)-dmax-1.2*BD/2;
x2 = point2(1)+dmax+1.2*BD/2;
y2 = point2(2)+dmax+1.2*BD/2;
z2 = point2(3)+dmax+1.2*BD/2;

% Same ascending obstacle IDs as the six-set intersection, without sorting.
index=find(data.mapx1(:)<=x2 & data.mapx2(:)>=x1 & ...
    data.mapy1(:)<=y2 & data.mapy2(:)>=y1 & ...
    data.mapz1(:)<=z2 & data.mapz2(:)>=z1);
U = 0;
U2=-100;


 if point2(1) >= map_bounds1(1) || point1(1) <= map_bounds2(1)|| ...
       point2(2)  >= map_bounds1(2)||point1(2) <= map_bounds2(2) ||...
       point2(3) >= map_bounds1(3) || point1(3) <= map_bounds2(3)
            U = Inf;
        return;
 end


if isempty(index),U2=0;return;end
point1=round(point1,decimal_places);point2=round(point2,decimal_places);
prepared=isfield(data,'dds_energy_boxes') && data.dds_energy_diameter==BD;
for i = 1:length(index)
    if prepared && index(i)<=size(data.dds_energy_boxes,1)
        box=data.dds_energy_boxes(index(i),:);
        xmin=box(1);ymin=box(2);zmin=box(3);xmax=box(4);ymax=box(5);zmax=box(6);
    else
    xmin = data.mapx1(index(i))-1.2*BD/2; ymin = data.mapy1(index(i))-1.2*BD/2; zmin =data.mapz1(index(i))-1.2*BD/2;
    xmax = data.mapx2(index(i))+1.2*BD/2; ymax = data.mapy2(index(i))+1.2*BD/2; zmax =data.mapz2(index(i))+1.2*BD/2;

xmin = round(xmin, decimal_places);
ymin = round(ymin, decimal_places);
zmin = round(zmin, decimal_places);
xmax = round(xmax, decimal_places);
ymax = round(ymax, decimal_places);
zmax = round(zmax, decimal_places);
    end
% 检查直线是否穿过障碍物
    if ~(point2(1) < xmin || point1(1)> xmax  ||point2(2) < ymin ||point1(2) > ymax || point2(3) < zmin  || point1(3) > zmax)
        U = Inf;
        return;
    end
    
    % 计算各轴上的距离分量
    dx = min(abs([xmin - x, x - xmax]));
    dy = min(abs([ymin - y, y - ymax]));
    dz = min(abs([zmin - z, z - zmax]));
    d1 = [dx,dy,dz];
    positive_values = d1(d1 > 0);
    d=min(positive_values);
    pipestyle=data.pipestyle;
    obstaclestlye=data.obstacles(index(i),4);
    obstaclestlye1=judgeobstaclestlye(pipestyle,obstaclestlye);
    emin=data.obstacles(index(i),1);
    emax=data.obstacles(index(i),2);
    emaxmin=data.obstacles(index(i),3);
    dbest=data.obstacles(index(i),5);
    dsafe=data.dsafe;
        if dmax==dsafe
        dmax=1.5*dsafe;
    end
    % 计算势能贡献
    if d <= dmax
        if obstaclestlye1==1%吸引
        U = U + -emin/(dmax-dsafe)*abs(d-dsafe)+emin; % 线性势能函数
        U2=max(U2,-emin/(dmax-dsafe)*(d-dsafe)+emin);
        else
         if obstaclestlye1==2%强排斥
           U = U + log(d/dmax)/log(0.1)+emax;
           U2=max(U2,log(d/dmax)/log(0.1)+emax);
         else
         if d>dbest
             U = U + emaxmin;
             U2=max(U2,emaxmin);
         else
             U = U + log(d/dbest)/log(0.1)+emaxmin;
             U2=max(U2,log(d/dbest)/log(0.1)+emaxmin);
         end
         end
        end
    end
end

if U2==-100
    U2=0;
end
end
