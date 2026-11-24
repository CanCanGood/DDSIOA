function aabb = computePipeAABB(curve1, curve2, BD)
    % 输入参数说明：
    % curve1: 管道中心线的第一个拐点坐标 [x1, y1, z1]
    % curve2: 管道中心线的第二个拐点坐标 [x2, y2, z2]
    % BD: 管道的直径
    % 输出：AABB包围盒 [xmin, ymin, zmin, xmax, ymax, zmax]
    count = find(abs(curve1-curve2)==0);
    index1=find(curve1~=curve2);
if curve1(index1)>curve2(index1)
    curve2=curve1+curve2;
    curve1=curve2-curve1;
    curve2=curve2-curve1;
end
curve1(count)=curve1(count)-BD/2;
curve2(count)=curve2(count)+BD/2;   
aabb=[curve1,curve2];
end




