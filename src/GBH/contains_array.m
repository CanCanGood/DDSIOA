function is_on_orthogonal_segment = contains_array(point, endpoint1, endpoint2)
% 判断一个点是否在由两个端点定义的三维正交线段上
%
% 输入：
%   point:      要判断的点的三维坐标 [x, y, z]
%   endpoint1:  线段的第一个端点的三维坐标 [x, y, z]
%   endpoint2:  线段的第二个端点的三维坐标 [x, y, z]
%
% 输出：
%   is_on_orthogonal_segment:  如果点在线段上，则为 true；否则为 false

% 1. 检查正交性并确定平行轴
if (endpoint1(1) ~= endpoint2(1)) && (endpoint1(2) == endpoint2(2)) && (endpoint1(3) == endpoint2(3))
    parallel_axis = 'x';
elseif (endpoint1(1) == endpoint2(1)) && (endpoint1(2) ~= endpoint2(2)) && (endpoint1(3) == endpoint2(3))
    parallel_axis = 'y';
elseif (endpoint1(1) == endpoint2(1)) && (endpoint1(2) == endpoint2(2)) && (endpoint1(3) ~= endpoint2(3))
    parallel_axis = 'z';
elseif (endpoint1(1) == endpoint2(1)) && (endpoint1(2) == endpoint2(2)) && (endpoint1(3) == endpoint2(3))
    parallel_axis = 'all';
   
else
    % 线段不是正交的
    is_on_orthogonal_segment = false;
    return;
end

% 2. 检查共线性和范围
switch parallel_axis
    case 'x'
        if (point(2) == endpoint1(2)) && (point(3) == endpoint1(3))
            min_x = min(endpoint1(1), endpoint2(1));
            max_x = max(endpoint1(1), endpoint2(1));
            is_on_orthogonal_segment = (point(1) >= min_x) && (point(1) <= max_x);
        else
            is_on_orthogonal_segment = false;
        end
    case 'y'
        if (point(1) == endpoint1(1)) && (point(3) == endpoint1(3))
            min_y = min(endpoint1(2), endpoint2(2));
            max_y = max(endpoint1(2), endpoint2(2));
            is_on_orthogonal_segment = (point(2) >= min_y) && (point(2) <= max_y);
        else
            is_on_orthogonal_segment = false;
        end
    case 'z'
        if (point(1) == endpoint1(1)) && (point(2) == endpoint1(2))
            min_z = min(endpoint1(3), endpoint2(3));
            max_z = max(endpoint1(3), endpoint2(3));
            is_on_orthogonal_segment = (point(3) >= min_z) && (point(3) <= max_z);
        else
            is_on_orthogonal_segment = false;
        end
    case 'all'
        if (point(2) == endpoint1(2)) && (point(3) == endpoint1(3)) && (point(1) == endpoint1(1))
           
            is_on_orthogonal_segment =1;
        else
            is_on_orthogonal_segment = false;
        end
end
end