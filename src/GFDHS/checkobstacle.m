function ob= checkobstacle(point1,point2,data,dmax,BD)
index1=find(point1~=point2);
point=point1;
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

indexx1=find(data.mapx1<x2);
indexx2=find(data.mapx2>x1);
indexy1=find(data.mapy1<y2);
indexy2=find(data.mapy2>y1);
indexz1=find(data.mapz1<z2);
indexz2=find(data.mapz2>z1);
index=intersect(intersect(intersect(intersect(intersect(indexx1, indexx2), indexy1),indexy2),indexz1),indexz2);
dmin=1000;
ob=[];
for i = 1:length(index)
    xmin = data.mapx1(index(i))-1.2*BD/2; ymin = data.mapy1(index(i))-1.2*BD/2; zmin =data.mapz1(index(i))-1.2*BD/2;
    xmax = data.mapx2(index(i))+1.2*BD/2; ymax = data.mapy2(index(i))+1.2*BD/2; zmax =data.mapz2(index(i))+1.2*BD/2;
    d1=min(abs(point-[xmin,ymin,zmin]));
    d2=min(abs(point-[xmax,ymax,zmax]));

if dmin>min(d1,d2)
    
    dmin=min(d1,d2);
    xvaluemin=[data.mapx1(index(i)),data.mapx2(index(i)),point1(1)+1.2*BD/2];
    xmin=min(xvaluemin);
    yvaluemin=[data.mapy1(index(i)),data.mapy2(index(i)),point1(2)+1.2*BD/2];
    ymin=min(yvaluemin);
    zvaluemin=[data.mapz1(index(i)),data.mapz2(index(i)),point1(3)+1.2*BD/2];
    zmin=min(zvaluemin);
    xvaluemax=[data.mapx1(index(i)),data.mapx2(index(i)),point2(1)-1.2*BD/2];
    xmax=max(xvaluemax);
    yvaluemax=[data.mapy1(index(i)),data.mapy2(index(i)),point2(2)-1.2*BD/2];
    ymax=max(yvaluemax);
    zvaluemax=[data.mapz1(index(i)),data.mapz2(index(i)),point2(3)-1.2*BD/2];
    zmax=max(zvaluemax);
    ob=[xmin,ymin,zmin,xmax,ymax,zmax];
end
end
if isempty(ob)
    ob=[dmax,dmax,dmax,dmax,dmax,dmax];
end


