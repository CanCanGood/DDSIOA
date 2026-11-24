%%是优化问题的目标函数，aimFcn_1，即option.fobj，主要用于实现路径查找或路径优化的功能，它具体实现了一个搜索算法，这个算法会根据当前位置和可能的下一步位置来确定其路径
function [fit,x,result] = aimFcn_astar(option,data,x3,x1,x2)
% aimFcn_astar
% 3D A* for voxel map with:
% - 6-neighborhood
% - Manhattan heuristic
% - bend-aware heuristic / turn penalty
% - hard obstacle from data.map_er == 88
% - soft energy from data.map_er ~= 88
%
% Inputs:
%   option, data, x3, x1, x2 : external inputs
%
% Required fields in data:
%   data.S       : start voxel coordinate [x,y,z]
%   data.E       : end voxel coordinate [x,y,z]
%   data.map_er  : 3D potential map, 88=hard obstacle, others=soft energy
%   data.dir     : 6x3 direction table
%
% Outputs:
%   fit
%   x
%   result

    %#ok<*NASGU>
    %#ok<*INUSD>

    % -----------------------------
    % 参数读取
    % -----------------------------
    map_er = data.map_er;
    sP     = round(data.S(:))';
    endP   = round(data.E(:))';

    if ~isfield(data,'dir') || isempty(data.dir)
        data.dir = [ 1, 0, 0; ...
                    -1, 0, 0; ...
                     0, 1, 0; ...
                     0,-1, 0; ...
                     0, 0, 1; ...
                     0, 0,-1 ];
    end
    dirs = data.dir;

    % -----------------------------
    % option参数，可按需从外部传入
    % -----------------------------
    w_step = getOpt(option,'w_step',1.0);          % 单步长度代价
    w_turn = getOpt(option,'w_turn',0.8);          % 转弯惩罚
    w_soft = getOpt(option,'w_soft',1.0);          % 势能代价权重
    w_fit_len   = getOpt(option,'w_fit_len',1.0);  % fit中长度权重
    w_fit_curve = getOpt(option,'w_fit_curve',1.0);% fit中弯头权重
    w_fit_energy= getOpt(option,'w_fit_energy',1.0);% fit中势能权重

    % 是否把起点终点强制裁剪到地图内
    clipToMap = getOpt(option,'clipToMap',true);

    mapSize = size(map_er);
    if numel(mapSize) < 3
        error('data.map_er 必须是三维矩阵');
    end

    if clipToMap
        sP   = clampPoint(sP,mapSize);
        endP = clampPoint(endP,mapSize);
    end

    % -----------------------------
    % 初始化输出
    % -----------------------------
    fit = inf;
    x = [];
    result = struct();

    % -----------------------------
    % 基本检查
    % -----------------------------
    if ~isInside(sP,mapSize) || ~isInside(endP,mapSize)
        result = fillFailResult(fit, sP);
        return;
    end

    if map_er(sP(1),sP(2),sP(3)) == 88
        result = fillFailResult(fit, sP);
        return;
    end

    if map_er(endP(1),endP(2),endP(3)) == 88
        result = fillFailResult(fit, sP);
        return;
    end

    if isequal(sP,endP)
        path1 = sP;
        [curve, xdir, xdir1, curve_count] = encodePath(path1, dirs);
        dis_len = 0;
        curve_len = size(curve,1);
        penetration_len = 0;
        energy_total = 0;
        palh1 = 0;

        fit = w_fit_len*dis_len + w_fit_curve*curve_count + w_fit_energy*energy_total;
        x = curve;

        result.fit = fit;
        result.length = dis_len;
        result.path = curve;
        result.path1 = path1;
        result.path2 = xdir;
        result.endnode = sP;
        result.dirpath = [sP, xdir];
        result.xdir = xdir;
        result.xdir1 = xdir1;
        result.dis_len = dis_len;
        result.curve_len = curve_len;
        result.curve = curve;
        result.penetration_len = penetration_len;
        result.energy_total = energy_total;
        result.curve_count = curve_count;
        result.palh1 = palh1;
        result.searchedVoxels =searchedVoxels;
    if exist('case1_iterations','var'), result.search_iterations=case1_iterations; end
        return;
    end

    % -----------------------------
    % A* 状态定义
    % 状态 = (x,y,z,lastDir)
    % lastDir:
    %   1~6 代表上一步方向
    %   7   代表起点无方向
    % -----------------------------
    START_DIR = 7;
    numDirStates = 7;

    gScore = inf([mapSize, numDirStates]);
    fScore = inf([mapSize, numDirStates]);
    closed = false([mapSize, numDirStates]);

    parentX   = zeros([mapSize, numDirStates],'int32');
    parentY   = zeros([mapSize, numDirStates],'int32');
    parentZ   = zeros([mapSize, numDirStates],'int32');
    parentDir = zeros([mapSize, numDirStates],'int32');

    sx = sP(1); sy = sP(2); sz = sP(3);
    gx = endP(1); gy = endP(2); gz = endP(3);

    gScore(sx,sy,sz,START_DIR) = 0;
    fScore(sx,sy,sz,START_DIR) = heuristicLB(sP,endP,START_DIR,w_turn);

    % openList: [x y z dir]
    openList = int32([sx,sy,sz,START_DIR]);
    inOpen = false([mapSize, numDirStates]);
    inOpen(sx,sy,sz,START_DIR) = true;

    found = false;
    finalState = [];
     searchedVoxels=[sx,sy,sz];
    % -----------------------------
    % A* 主循环
    % -----------------------------
    case1_timer=tic; case1_iterations=0;
    while ~isempty(openList)
        case1_iterations=case1_iterations+1;
        % No elapsed-time cutoff: explicitly requested by the user.
        % 取 f 最小节点
        bestIdx = 1;
        bestF = inf;
        bestH = inf;
        bestG = -inf;

        for i = 1:size(openList,1)
            cx = openList(i,1);
            cy = openList(i,2);
            cz = openList(i,3);
            cd = openList(i,4);

            curF = fScore(cx,cy,cz,cd);
            curG = gScore(cx,cy,cz,cd);
            curH = curF - curG;

            if (curF < bestF) || ...
               (curF == bestF && curH < bestH) || ...
               (curF == bestF && curH == bestH && curG > bestG)
                bestF = curF;
                bestH = curH;
                bestG = curG;
                bestIdx = i;
            end
        end

        current = openList(bestIdx,:);
        openList(bestIdx,:) = [];

        cx = current(1); cy = current(2); cz = current(3); cdir = current(4);
        inOpen(cx,cy,cz,cdir) = false;

        if closed(cx,cy,cz,cdir)
            continue;
        end
        closed(cx,cy,cz,cdir) = true;

        if cx == gx && cy == gy && cz == gz
            found = true;
            finalState = [cx,cy,cz,cdir];
            break;
        end

        % 扩展6邻域
        for ndir = 1:6

            nx = cx + dirs(ndir,1);
            ny = cy + dirs(ndir,2);
            nz = cz + dirs(ndir,3);
            searchedVoxels=[searchedVoxels;nx,ny,nz];
            if ~isInside([nx,ny,nz],mapSize)
                continue;
            end

            % 硬障碍
            if map_er(nx,ny,nz) == 88
                continue;
            end

            if closed(nx,ny,nz,ndir)
                continue;
            end

            % 基础步长
            stepCost = w_step;

            % 转弯代价
            if cdir == START_DIR || cdir == ndir
                turnCost = 0;
            else
                turnCost = w_turn;
            end

            % 软势能代价
            softEnergy = map_er(nx,ny,nz);
            if softEnergy == 88
                softCost = inf;
            else
                softCost = w_soft * max(0, softEnergy);
            end

            tentativeG = gScore(cx,cy,cz,cdir) + stepCost + turnCost + softCost;

            if tentativeG < gScore(nx,ny,nz,ndir)
                gScore(nx,ny,nz,ndir) = tentativeG;
                h = heuristicLB([nx,ny,nz], endP, ndir, w_turn);
                fScore(nx,ny,nz,ndir) = tentativeG + h;

                parentX(nx,ny,nz,ndir)   = cx;
                parentY(nx,ny,nz,ndir)   = cy;
                parentZ(nx,ny,nz,ndir)   = cz;
                parentDir(nx,ny,nz,ndir) = cdir;

                if ~inOpen(nx,ny,nz,ndir)
                    openList(end+1,:) = int32([nx,ny,nz,ndir]); %#ok<AGROW>
                    inOpen(nx,ny,nz,ndir) = true;
                end
            end
        end
    end

    if ~found
        result = fillFailResult(fit, sP);
        return;
    end

    % -----------------------------
    % 回溯完整路径 path1
    % -----------------------------
    [path1, dirSeq] = reconstructPath(finalState, sP, parentX, parentY, parentZ, parentDir);

    % -----------------------------
    % 计算路径指标
    % -----------------------------
    [curve, xdir, xdir1, curve_count] = encodePath(path1, dirs);

    dis_len = size(path1,1) - 1;     % 总步长（格）
    curve_len = size(curve,1);       % 拐点数量（含起终点）
    palh1 = 0;

    % penetration_len: 经过软约束网格的步数（势能>0且不是88）
    penetration_len = 0;
    energy_total = 0;

    for i = 2:size(path1,1)  % 从第二个点开始累计
        px = path1(i,1);
        py = path1(i,2);
        pz = path1(i,3);
        val = map_er(px,py,pz);

        if val ~= 88 && val > 0
            penetration_len = penetration_len + 1;
            energy_total = energy_total + val;
        end
    end

    % -----------------------------
    % 适应度
    % -----------------------------
    fit = w_fit_len   * dis_len + ...
          w_fit_curve * curve_count + ...
          w_fit_energy* energy_total;

    % 这里把 x 设为“拐点路径”
    x =xdir;
    path1 = [];
    for i = 1:size(curve,1)-1
        aabb = computePipeAABB(curve(i,:)*data.unit, curve(i+1,:)*data.unit, data.BD1);
        path1 = [path1; aabb];
    end
    % -----------------------------
    % 输出结果
    % -----------------------------
    result.fit = fit;
    result.length = dis_len;
    result.path = curve;
    result.path1 = path1;
    result.path2 = xdir;
    result.endnode = sP;          % 按你给的字段要求原样保留
    result.dirpath = [sP, xdir];
    result.xdir = xdir;
    result.xdir1 = xdir1;
    result.dis_len = dis_len;
    result.curve_len = curve_len;
    result.curve = curve;
    result.penetration_len = penetration_len;
    result.energy_total = energy_total;
    result.curve_count = curve_count;
    result.palh1 = palh1;
    result.searchedVoxels =searchedVoxels;
    if exist('case1_iterations','var'), result.search_iterations=case1_iterations; end

end

% =========================================================================
% 工具函数
% =========================================================================

function v = getOpt(option, name, defaultValue)
    v = defaultValue;
    if isstruct(option) && isfield(option, name)
        v = option.(name);
    end
end

function p = clampPoint(p, mapSize)
    p(1) = min(max(1, p(1)), mapSize(1));
    p(2) = min(max(1, p(2)), mapSize(2));
    p(3) = min(max(1, p(3)), mapSize(3));
end

function tf = isInside(p, mapSize)
    tf = p(1) >= 1 && p(1) <= mapSize(1) && ...
         p(2) >= 1 && p(2) <= mapSize(2) && ...
         p(3) >= 1 && p(3) <= mapSize(3);
end

function h = heuristicLB(pos, goal, currentDir, w_turn)
% 曼哈顿距离 + 转弯下界
    dx = goal(1) - pos(1);
    dy = goal(2) - pos(2);
    dz = goal(3) - pos(3);

    manhattan = abs(dx) + abs(dy) + abs(dz);

    needDirs = [];
    if dx > 0
        needDirs(end+1) = 1; %#ok<AGROW>
    elseif dx < 0
        needDirs(end+1) = 2; %#ok<AGROW>
    end

    if dy > 0
        needDirs(end+1) = 3; %#ok<AGROW>
    elseif dy < 0
        needDirs(end+1) = 4; %#ok<AGROW>
    end

    if dz > 0
        needDirs(end+1) = 5; %#ok<AGROW>
    elseif dz < 0
        needDirs(end+1) = 6; %#ok<AGROW>
    end

    k = numel(needDirs);

    if k == 0
        minTurnLB = 0;
    elseif currentDir == 7
        minTurnLB = max(k - 1, 0);
    else
        if any(needDirs == currentDir)
            minTurnLB = k - 1;
        else
            minTurnLB = k;
        end
    end

    h = manhattan + w_turn * minTurnLB;
end

function [path1, dirSeq] = reconstructPath(finalState, startP, parentX, parentY, parentZ, parentDir)
% 返回完整网格路径 path1, 以及每步对应的方向序列 dirSeq
    x = finalState(1);
    y = finalState(2);
    z = finalState(3);
    d = finalState(4);

    revPath = zeros(0,3);
    revDir  = zeros(0,1);

    while true
        revPath(end+1,:) = [x,y,z]; %#ok<AGROW>
        revDir(end+1,1)  = d; %#ok<AGROW>

        if x == startP(1) && y == startP(2) && z == startP(3) && d == 7
            break;
        end

        px = parentX(x,y,z,d);
        py = parentY(x,y,z,d);
        pz = parentZ(x,y,z,d);
        pd = parentDir(x,y,z,d);

        x = px; y = py; z = pz; d = pd;
    end

    path1 = flipud(revPath);
    dirSeq = flipud(revDir);

    % 去掉起点的虚拟方向7
    if ~isempty(dirSeq) && dirSeq(1) == 7
        dirSeq(1) = [];
    end
end

function [curve, xdir, xdir1, curve_count] = encodePath(path1, dirs)
% path1  : 完整路径 N x 3
% curve  : 拐点路径
% xdir   : [dirIndex, runLength]
% xdir1  : [dx,dy,dz,runLength]
% curve_count : 弯头数

    if isempty(path1)
        curve = [];
        xdir = [];
        xdir1 = [];
        curve_count = inf;
        return;
    end

    if size(path1,1) == 1
        curve = path1;
        xdir =[];
        xdir1 = zeros(0,4);
        curve_count = 0;
        return;
    end
    xdir =[];
    % 每一步的差分方向
    dP = diff(path1,1,1);   % (N-1) x 3

    dirIdSeq = zeros(size(dP,1),1);
    for i = 1:size(dP,1)
        for k = 1:size(dirs,1)
            if isequal(dP(i,:), dirs(k,:))
                dirIdSeq(i) = k;
                break;
            end
        end
    end

    % 压缩成 [方向编号, 连续长度]

    curDir = dirIdSeq(1);
    runLen = 1;
    
    for i = 2:numel(dirIdSeq)
        if dirIdSeq(i) == curDir
            runLen = runLen + 1;
        else
            if isempty(xdir)
            xdir= [curDir, runLen]; %#ok<AGROW>
            xdir1= [dirs(curDir,:), runLen]; %#ok<AGROW>
            curDir = dirIdSeq(i);
            runLen = 1;
            else
            xdir(end+1:end+2) = [curDir, runLen]; %#ok<AGROW>
            xdir1(end+1,:) = [dirs(curDir,:), runLen]; %#ok<AGROW>
            curDir = dirIdSeq(i);
            runLen = 1;
            end
        end
    end
    xdir(end+1:end+2) = [curDir, runLen]; %#ok<AGROW>
    xdir1(end+1,:) = [dirs(curDir,:), runLen];

    % 拐点提取
    curve = path1(1,:);
    cursor = path1(1,:);
    for i = 1:2:length(xdir)
        moveVec = dirs(xdir(i),:) * xdir(i+1);
        cursor = cursor + moveVec;
        curve(end+1,:) = cursor; %#ok<AGROW>
    end

    % 弯头数量 = 段数 - 1
    curve_count = max(length(xdir)/2- 1, 0);
end

function result = fillFailResult(fit, sP)
    result = struct();
    result.fit = inf;
    result.length = inf;
    result.path = [];
    result.path1 = [];
    result.path2 = [];
    result.endnode = sP;
    result.dirpath = [];
    result.xdir = [];
    result.xdir1 = [];
    result.dis_len = inf;
    result.curve_len = inf;
    result.curve = [];
    result.penetration_len = inf;
    result.energy_total = inf;
    result.curve_count = inf;
    result.palh1 =0;
end