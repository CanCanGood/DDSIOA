function [markedGroups, filteredGroups] = filterAndMarkGroups(data)
    % 输入: data - n*3 的二维数组，三行为一组
    % 输出:
    %   markedGroups - 带有标记的原始组
    %   filteredGroups - 筛选后的直径最大或曼哈顿距离最长的组
    
    % 初始化
    numGroups = size(data, 1) / 3; % 计算组数
    groups = reshape(data', [3, 3, numGroups]); % 将数据重塑为 3x3x组数
    groups = permute(groups, [2, 1, 3]); 
    % 存储已处理的组索引和最终结果
    markedGroups = zeros(size(data,1),size(data,2)+1); % 用于标记原始组
    filteredGroups = [];
    usedIndices = []; % 记录已处理的组索引
    
    % 标记号
    currentMark = 1;
    
    % 遍历所有组，寻找起点和终点有重复的组
    for i = 1:numGroups
        if ismember(i, usedIndices)
            continue; % 如果该组已经处理，跳过
        end
        
        % 当前组的起点和终点
        startPoint_i = groups(1, :, i);
        endPoint_i = groups(2, :, i);
        
        % 存储与当前组重复的所有组
        candidateGroups = groups(:, :, i);
        candidateIndices = i;
        
        for j = i+1:numGroups
            if ismember(j, usedIndices)
                continue;
            end
            
            startPoint_j = groups(1, :, j);
            endPoint_j = groups(2, :, j);
            
            % 判断起点或终点是否有重复
            if isequal(startPoint_i, startPoint_j) || isequal(startPoint_i, endPoint_j) || ...
               isequal(endPoint_i, startPoint_j) || isequal(endPoint_i, endPoint_j)
                candidateGroups = cat(3, candidateGroups, groups(:, :, j));
                candidateIndices = [candidateIndices, j];
            end
        end
        
        % 标记所有找到的重复组
        for idx = candidateIndices
            markedGroups((idx-1)*3+1:idx*3, :) = [groups(:, :, idx), currentMark * ones(3, 1)];
        end
        
        % 筛选直径最大或曼哈顿距离最长的组
        maxDiameter = -inf;
        maxDistance = -inf;
        bestGroup = [];
        
        for k = 1:size(candidateGroups, 3)
            diameter = candidateGroups(3, 1, k); % 直径（第三行第一列）
            
            % 计算曼哈顿距离
            startPoint = candidateGroups(1, :, k);
            endPoint = candidateGroups(2, :, k);
            manhattanDistance = sum(abs(startPoint - endPoint));
            
            % 筛选条件：直径最大，若直径相同则曼哈顿距离最长
            if diameter > maxDiameter || (diameter == maxDiameter && manhattanDistance > maxDistance)
                maxDiameter = diameter;
                maxDistance = manhattanDistance;
                bestGroup = candidateGroups(:, :, k);
            end
        end
        
        % 存储筛选后的最佳组
        filteredGroups = cat(3, filteredGroups, bestGroup);
        usedIndices = [usedIndices, candidateIndices];
        
        % 更新标记号
        currentMark = currentMark + 1;
    end
end
