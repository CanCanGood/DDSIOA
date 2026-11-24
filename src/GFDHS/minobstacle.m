function ob= minobstacle(point1,point2,data,dmax,BD,path)
if nargin==6&&size(path,2)>=6
data.mapx1=[data.mapx1;min(path(:,1),path(:,4))];
data.mapx2=[data.mapx2;max(path(:,1),path(:,4))];
data.mapy1=[data.mapy1;min(path(:,2),path(:,5))];
data.mapy2=[data.mapy2;max(path(:,2),path(:,5))];
data.mapz1=[data.mapz1;min(path(:,3),path(:,6))];
data.mapz2=[data.mapz2;max(path(:,3),path(:,6))];
data.obstacles=[data.obstacles;path(:,end-4:end)];
end
point=point1;
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


    ob=[xmax-point(1);point(1)-xmin;ymax-point(2);point(2)-ymin;zmax-point(3);point(3)-zmin];
    index1=find(ob<=0);
    for j=1:length(index1)
    ob(index1(j)+mod((index1(j)),2)*2-1)=abs(ob(index1(j)));
    ob(index1(j))=BD;
    end
tolerance = 0.0001;
ob(ob <= tolerance & ob >= 0) = 1;
end
end
if isempty(ob)
    ob=[dmax;dmax;dmax;dmax;dmax;dmax];
end
