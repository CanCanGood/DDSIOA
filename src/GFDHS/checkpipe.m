function [U,U2,n]=checkpipe(point1,point2, data, dmax,BD)
% 输入参数：
%   point: 点的坐标 [x, y, z]
%   map_bounds: 地图边界 [a, b, c]（本函数未使用该参数）
%   obstacles: N×6矩阵，每行表示一个障碍物的[xmin, ymin, zmin, xmax, ymax, zmax]
%   dmax: 障碍物影响的最大距离
% 输出：
%   U: 障碍物对点的总势能（点不在障碍物内且距离超过dmax时为0）
if isfield(data,'energy_range_pipe')
    dmax=data.energy_range_pipe;
elseif ~isfield(data,'unified_dmax') || ~data.unified_dmax
    dmax=25; % Default interaction range.
end
n=0;
decimal_places = 10; % 保留4位小数
x = point2(1);
y = point2(2);
z = point2(3);
index1=find(point1~=point2);
if point1(index1)>point2(index1)
    point2=point1+point2;
    point1=point2-point1;
    point2=point2-point1;
end
x1 = point1(1)-dmax;
y1 = point1(2)-dmax; 
z1 = point1(3)-dmax;
x2 = point2(1)+dmax;
y2 = point2(2)+dmax;
z2 = point2(3)+dmax;

% Preserve strict inequalities and ascending IDs of the original intersection.
index=find(data.mapx11(:)<x2 & data.mapx21(:)>x1 & ...
    data.mapy11(:)<y2 & data.mapy21(:)>y1 & ...
    data.mapz11(:)<z2 & data.mapz21(:)>z1);
U = 0;
U2=0;
if isempty(index),return;end
point1=round(point1,decimal_places);point2=round(point2,decimal_places);
prepared=isfield(data,'dds_pipe_boxes') && data.dds_energy_diameter==BD;
for i = 1:length(index)
    if prepared && index(i)<=size(data.dds_pipe_boxes,1)
        box=data.dds_pipe_boxes(index(i),:);
        xmin=box(1);ymin=box(2);zmin=box(3);xmax=box(4);ymax=box(5);zmax=box(6);
    else
    xmin = data.mapx11(index(i))-1.1*BD/2; ymin = data.mapy11(index(i))-1.1*BD/2; zmin =data.mapz11(index(i))-1.1*BD/2;
    xmax = data.mapx21(index(i))+1.1*BD/2; ymax = data.mapy21(index(i))+1.1*BD/2; zmax =data.mapz21(index(i))+1.1*BD/2;

xmin = round(xmin, decimal_places);
ymin = round(ymin, decimal_places);
zmin= round(zmin, decimal_places);
xmax = round(xmax, decimal_places);
ymax = round(ymax, decimal_places);
zmax = round(zmax, decimal_places);
    end
% 检查直线是否穿过障碍物
    if ~(point2(1) <= xmin || point1(1)>= xmax  ||point2(2) <= ymin ||point1(2) >= ymax || point2(3) <= zmin  || point1(3) >= zmax)
        U =max(U,0) +10;
        continue;
    end
    
    % 计算各轴上的距离分量
    dx = min(abs([xmin - x, x - xmax]));
    dy = min(abs([ymin - y, y - ymax]));
    dz = min(abs([zmin - z, z - zmax]));
    d1 = [dx,dy,dz];
    positive_values = d1(d1 > 0);
    d=min(positive_values);
    pipestyle=data.pipestyle;
    obstaclestlye=data.pipestyle1(index(i),4);
    obstaclestlye1=judgeobstaclestlye(pipestyle,obstaclestlye);
    emin=data.pipestyle1(index(i),1);
    emax=data.pipestyle1(index(i),2);
    emaxmin=data.pipestyle1(index(i),3);
    dbest=data.pipestyle1(index(i),5);
    dsafe=data.BD1;
    if dmax==dsafe
        dmax=1.5*dsafe;
    end
    % 计算势能贡献
    if index(i)<=2*data.numpipe
        U =max(U,0);
    else
    if U>=10
    if d <= dmax
        if obstaclestlye1==1%吸引
        U2=max(U2,-emin/(dmax-dsafe)*(d-dsafe)+emin);
        n=n+1;
        else
         if obstaclestlye1==2%强排斥
           U =max( U,U + log(d/dmax)+emax);
           U2=max(U2,log(d/dmax)+emax);
         else
         if d>dbest
             U =max( U, U + emaxmin);
             U2=max(U2,emaxmin);
         else
             U =max( U, U + log(d/dbest)+emaxmin);
             U2=max(U2,log(d/dbest)+emaxmin);
             if log(d/dbest)+emaxmin<0
                n=n+1; 
             end
         end

         end
        end
    end
    else
        if d <= dmax
        if obstaclestlye1==1%吸引
        U = U + emin/(dmax-dsafe)*abs(d-dsafe); % 线性势能函数
        U2=max(U2,-emin/(dmax-dsafe)*(d-dsafe));
         n=n+1;
        else
         if obstaclestlye1==2%强排斥
           U = U + log(d/dmax)+emax;
           U2=max(U2,log(d/dmax)+emax);
         else
         if d>dbest
             U = U + emaxmin;
             U2=max(U2,emaxmin);
         else
             U = U + log(d/dbest)+emaxmin;
             U2=max(U2,log(d/dbest)+emaxmin);
                 if log(d/dbest)+emaxmin<0
                n=n+1; 
             end
         end
         end
        end
    end    
    end
    end
end
end


