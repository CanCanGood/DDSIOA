function [flag]=chongdie(zb1,zb2,zb3,zb4,dz)
%包含关系
flag=false;
if (min(zb1,zb2)<=min(zb3,zb4)&&max(zb1,zb2)>=max(zb3,zb4))||(min(zb1,zb2)>=min(zb3,zb4)&&max(zb1,zb2)<=max(zb3,zb4))
    flag=true;
    return;
end
%交叉关系
if (min(zb1,zb2)<max(zb3,zb4)&&max(zb1,zb2)>max(zb3,zb4))||(max(zb1,zb2)>min(zb3,zb4)&&max(zb1,zb2)<max(zb3,zb4))
    flag=true;
    return;
end
%分离关系
if min(abs(min(zb1,zb2)-max(zb3,zb4)),abs(min(zb3,zb4)-max(zb1,zb2)))<dz
    flag=true;
    return;
end
end