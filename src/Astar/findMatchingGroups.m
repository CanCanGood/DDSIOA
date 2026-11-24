function fourthColumnValues = findMatchingGroups(markedGroups, selected_data)
    % 输入:
    % markedGroups - n x 4 的二维数组
    % selected_data - n x 4 的二维数组，用于比较的标准
    
    % 输出:
    % fourthColumnValues - 匹配组的第四列值列表

    % 确定组数
    numGroups = size(markedGroups, 1) / 3;
    
    % 初始化输出数组
    fourthColumnValues = [];
    compareGroup = selected_data;   
    % 遍历每个组
    for i = 1:numGroups
        % 提取当前组的数据
        groupStartIdx = (i-1) * 3 + 1;
        groupEndIdx = i * 3;
        currentGroup = markedGroups(groupStartIdx:groupEndIdx, :);
        
        % 比较前三列是否相等
        if isequal(currentGroup(:, 1:3), compareGroup(:, 1:3))
            % 如果相等，收集该组的第四列数据
            fourthColumnValues = currentGroup(1, 4);
            return 
        end
    end
end
