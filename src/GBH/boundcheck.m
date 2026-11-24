function [U]= boundcheck(point1,point2, map_bounds1,map_bounds2)
% 输入参数：
%   point: 点的坐标 [x, y, z]
%   map_bounds: 地图边界 [a, b, c]（本函数未使用该参数）
%   obstacles: N×6矩阵，每行表示一个障碍物的[xmin, ymin, zmin, xmax, ymax, zmax]
%   dmax: 障碍物影响的最大距离
% 输出：
%   U: 障碍物对点的总势能（点不在障碍物内且距离超过dmax时为0）
index1=find(point1~=point2);
if point1(index1)>point2(index1)
    point2=point1+point2;
    point1=point2-point1;
    point2=point2-point1;
end


U = 0;

 if point2(1) >= map_bounds1(1) || point1(1) <= map_bounds2(1)|| ...
       point2(2)  >= map_bounds1(2)||point1(2) <= map_bounds2(2) ||...
       point2(3) >= map_bounds1(3) || point1(3) <= map_bounds2(3)
            U = Inf;
        return;
 end