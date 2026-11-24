function [flag]=chonghe(xb1,xb2,xb3,xb4,yb1,yb2,yb3,yb4,dz,dz2)
flag=false;
%完全不重合
if((max(xb1,xb2)<min(xb3,xb4)||max(xb3,xb4)<min(xb1,xb2))&&(max(yb1,yb2)<min(yb3,yb4)||max(yb3,yb4)<min(yb1,yb2)))
    return;
end
%重合
d=0;
x_overlap = max(0, min(max(xb1,xb2), max(xb3,xb4)) - max(min(xb1,xb2),  min(xb3,xb4)));

% 计算两条线段在y轴上的投影的重合部分长度
y_overlap = max(0, min(max(yb1,yb2), max(yb3,yb4)) - max(min(yb1,yb2),  min(yb3,yb4)));

% 返回两条线段在x轴和y轴上的投影的重合部分长度相加
overlap_length = x_overlap + y_overlap;
length=max((abs(xb1-xb2)+abs(yb1-yb2)),(abs(xb3-xb4)+abs(yb3-yb4)));
if x_overlap==0
    d=-min(max(xb1,xb2), max(xb3,xb4))+max(min(xb1,xb2),min(xb3,xb4));
end
if y_overlap==0
    d=-min(max(yb1,yb2), max(yb3,yb4))+max(min(yb1,yb2),min(yb3,yb4));
end
loss=overlap_length/length-d/(dz*length);
if loss>dz2
    flag=true;
end

end