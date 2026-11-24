function result = judgeobstaclestlye(pipestyle, obstaclestlye)
    % 计算两个参数的绝对差值
    diff = abs(pipestyle - obstaclestlye);
    
    % 根据差值判断输出结果
    if diff < 5
        result = 1;
    elseif diff > 5
        result = 2;
    else
        result = 3;
    end
end