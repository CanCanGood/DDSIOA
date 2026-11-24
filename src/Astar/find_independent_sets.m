function independent_sets = find_independent_sets(fenzu)
    % 获取每个 fenzutemp 的元素数量（三行为一组）
    numElements = cellfun(@(x) size(x, 1) / 3, {fenzu.fenzutemp});

    % 记录原始顺序索引
    originalIdx = 1:length(fenzu);

    % 对 fenzu 按元素数量升序排序（仅用于比较）
    [~, sortIdx] = sort(numElements, 'ascend');
    fenzu_sorted = fenzu(sortIdx);

    % 遍历字段，删除元素少的字段中的重复元素
    for i = 1:length(fenzu_sorted)-1
        for j = i+1:length(fenzu_sorted)
            % 获取两个字段的数据
            smallGroup = fenzu_sorted(i).fenzutemp;
            largeGroup = fenzu_sorted(j).fenzutemp;

            % 将两个字段的数据按三行一组进行重塑
            smallGroup_reshaped = reshape(smallGroup', 9,  [])';
            largeGroup_reshaped = reshape(largeGroup', 9,  [])';

            % 找出重复的元素（按行比较）
            [~, ia, ~] = intersect(smallGroup_reshaped, largeGroup_reshaped, 'rows');

            % 删除 smallGroup 中的重复元素
            if ~isempty(ia)
                smallGroup_reshaped(ia, :) = [];
            end

            % 更新 smallGroup 回到 fenzu_sorted 中
            fenzu_sorted(i).fenzutemp = reshape(smallGroup_reshaped', 3, [])';
        end
    end

    % 删除空字段
    fenzu_sorted = fenzu_sorted(~cellfun(@isempty, {fenzu_sorted.fenzutemp}));

    % 恢复原始顺序（根据 sortIdx 映射）
    [~, revertIdx] = sort(sortIdx);  
    temp = cell(1, length(fenzu));  
    for k = 1:length(fenzu_sorted)
        temp{revertIdx(k)} = fenzu_sorted(k);
    end
    % 删除空单元（即原来空的）
    independent_sets = [temp{~cellfun(@isempty, temp)}];

    return;
end
