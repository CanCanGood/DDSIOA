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

    global NATIVE_ASTAR_COUNTS
    if isempty(NATIVE_ASTAR_COUNTS),NATIVE_ASTAR_COUNTS=zeros(1,3);end
    if nargin==3
        assert(NATIVE_ASTAR_COUNTS(1)==6*(NATIVE_ASTAR_COUNTS(2)+1),'StrictSix:CandidateCount');
        NATIVE_ASTAR_COUNTS(2)=NATIVE_ASTAR_COUNTS(2)+1;
    else
        NATIVE_ASTAR_COUNTS(1)=NATIVE_ASTAR_COUNTS(1)+1;
    end
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
        return;
    end

    % -----------------------------
    % A* 状态定义
    % 状态 = (x,y,z,lastDir)
    % lastDir:
    %   1~6 代表上一步方向
    %   7   代表起点无方向
    % -----------------------------
    if nargin == 3
        % Evaluate the selected direction-length code; never launch another search.
        assert(isnumeric(x3) && isvector(x3) && mod(numel(x3),2)==0 && all(isfinite(x3)), 'StrictSix:InvalidCode');
        encoded=x3(:)';
        assert(all(encoded(1:2:end)>=1 & encoded(1:2:end)<=size(dirs,1) & encoded(1:2:end)==fix(encoded(1:2:end))), 'StrictSix:InvalidDirection');
        assert(all(encoded(2:2:end)>=0 & encoded(2:2:end)==fix(encoded(2:2:end))), 'StrictSix:InvalidLength');
        path1=zeros(sum(encoded(2:2:end))+1,3);path1(1,:)=sP;cursor=1;
        for segment=1:2:numel(encoded)
            for step_index=1:encoded(segment+1)
                cursor=cursor+1;path1(cursor,:)=path1(cursor-1,:)+dirs(encoded(segment),:);
            end
        end
        assert(isequal(path1(end,:),endP), 'StrictSix:EndpointMismatch');
        searchedVoxels=zeros(0,3); % Evaluation generates no search points.
    else
    START_DIR = 7;

    % Retain records only for discovered states.
    % A state key is the original MATLAB linear (x,y,z,lastDir) index.
    % uint64 keys avoid rounding large linear indices through double.
    strideY = uint64(mapSize(1));
    strideZ = strideY * uint64(mapSize(2));
    strideDir = strideZ * uint64(mapSize(3));
    stateLookup = containers.Map('KeyType','uint64','ValueType','uint32');

    stateCapacity = 4096;
    stateCount = 1;
    stateCoords = zeros(stateCapacity,4,'int32');
    gScore = inf(stateCapacity,1);                 % Costs remain double.
    fScore = inf(stateCapacity,1);
    closed = false(stateCapacity,1);
    parentIndex = zeros(stateCapacity,1,'uint32');
    heapPosition = zeros(stateCapacity,1,'uint32'); % Zero means not in Open.
    openOrder = zeros(stateCapacity,1,'uint64');

    % Each state has one Open entry. Updating a key retains its original
    % insertion order; this is the final tie break of the original scan.
    heap = zeros(stateCapacity,1,'uint32');
    heapCount = 1;
    nextOpenOrder = uint64(1);

    sx = sP(1); sy = sP(2); sz = sP(3);
    gx = endP(1); gy = endP(2); gz = endP(3);

    stateCoords(1,:) = int32([sx,sy,sz,START_DIR]);
    stateLookup(stateKey(sx,sy,sz,START_DIR)) = uint32(1);
    gScore(1) = 0;
    fScore(1) = heuristicLB(sP,endP,START_DIR,w_turn);
    heap(1) = uint32(1);
    heapPosition(1) = uint32(1);
    openOrder(1) = nextOpenOrder;

    found = false;
    finalIndex = uint32(0);
    % The legacy trace contains the start followed by every neighbour
    % considered, including out-of-map, blocked and already-closed points.
    traceCapacity = 65536;
    traceCount = 1;
    searchedVoxels = zeros(traceCapacity,3,'int32');
    searchedVoxels(1,:) = int32([sx,sy,sz]);
    % -----------------------------
    % A* 主循环
    % -----------------------------
    while heapCount > 0
        currentIndex = popOpen();
        current = stateCoords(currentIndex,:);

        cx = current(1); cy = current(2); cz = current(3); cdir = current(4);
        if closed(currentIndex)
            continue;
        end
        closed(currentIndex) = true;

        if cx == gx && cy == gy && cz == gz
            found = true;
            finalIndex = currentIndex;
            break;
        end

        % 扩展6邻域
        for ndir = 1:6

            nx = cx + dirs(ndir,1);
            ny = cy + dirs(ndir,2);
            nz = cz + dirs(ndir,3);
            traceCount = traceCount + 1;
            if traceCount > traceCapacity
                traceCapacity = 2 * traceCapacity;
                searchedVoxels(traceCapacity,3) = int32(0);
            end
            searchedVoxels(traceCount,:) = [nx,ny,nz];
            if ~isInside([nx,ny,nz],mapSize)
                continue;
            end

            % 硬障碍
            if map_er(nx,ny,nz) == 88
                continue;
            end

            key = stateKey(nx,ny,nz,ndir);
            if isKey(stateLookup,key)
                neighbourIndex = stateLookup(key);
            else
                neighbourIndex = uint32(0);
            end

            if neighbourIndex ~= 0 && closed(neighbourIndex)
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

            tentativeG = gScore(currentIndex) + stepCost + turnCost + softCost;

            oldG = inf;
            if neighbourIndex ~= 0
                oldG = gScore(neighbourIndex);
            end
            if tentativeG < oldG
                if neighbourIndex == 0
                    stateCount = stateCount + 1;
                    if stateCount > stateCapacity
                        previousCapacity = stateCapacity;
                        stateCapacity = 2 * stateCapacity;
                        stateCoords(stateCapacity,4) = int32(0);
                        gScore(previousCapacity+1:stateCapacity,1) = inf;
                        fScore(previousCapacity+1:stateCapacity,1) = inf;
                        closed(stateCapacity,1) = false;
                        parentIndex(stateCapacity,1) = uint32(0);
                        heapPosition(stateCapacity,1) = uint32(0);
                        openOrder(stateCapacity,1) = uint64(0);
                        heap(stateCapacity,1) = uint32(0);
                    end
                    neighbourIndex = uint32(stateCount);
                    stateLookup(key) = neighbourIndex;
                    stateCoords(neighbourIndex,:) = int32([nx,ny,nz,ndir]);
                end
                gScore(neighbourIndex) = tentativeG;
                h = heuristicLB([nx,ny,nz], endP, ndir, w_turn);
                fScore(neighbourIndex) = tentativeG + h;
                parentIndex(neighbourIndex) = currentIndex;

                if heapPosition(neighbourIndex) == 0
                    nextOpenOrder = nextOpenOrder + uint64(1);
                    openOrder(neighbourIndex) = nextOpenOrder;
                    heapCount = heapCount + 1;
                    heap(heapCount) = neighbourIndex;
                    heapPosition(neighbourIndex) = uint32(heapCount);
                    siftUp(heapCount);
                else
                    % f can remain equal after floating-point rounding;
                    % f-g can then change in either direction. Rebalance
                    % both ways rather than assuming a pure decrease-key.
                    siftUp(double(heapPosition(neighbourIndex)));
                    siftDown(double(heapPosition(neighbourIndex)));
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
    searchedVoxels = searchedVoxels(1:traceCount,:);
    [path1, dirSeq] = reconstructVisitedPath(finalIndex, sP, stateCoords, parentIndex);
    end % Search versus evaluation-only branch

    % -----------------------------
    % 计算路径指标
    % -----------------------------
    NATIVE_ASTAR_COUNTS(3)=NATIVE_ASTAR_COUNTS(3)+size(searchedVoxels,1);
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

    % Nested heap helpers modify shared arrays without passing/copying the
    % whole visited-state workspace on every operation.
    function keyValue = stateKey(px,py,pz,pd)
        keyValue = uint64(px-1) + strideY*uint64(py-1) + ...
            strideZ*uint64(pz-1) + strideDir*uint64(pd-1) + uint64(1);
    end

    function before = precedes(leftIndex,rightIndex)
        leftF = fScore(leftIndex);
        rightF = fScore(rightIndex);
        leftG = gScore(leftIndex);
        rightG = gScore(rightIndex);
        % Do not substitute the separately evaluated heuristic: the
        % original comparator calculates h as f-g, including rounding.
        leftH = leftF - leftG;
        rightH = rightF - rightG;
        before = (leftF < rightF) || ...
            (leftF == rightF && leftH < rightH) || ...
            (leftF == rightF && leftH == rightH && leftG > rightG) || ...
            (leftF == rightF && leftH == rightH && leftG == rightG && ...
             openOrder(leftIndex) < openOrder(rightIndex));
    end

    function swapHeap(leftPosition,rightPosition)
        leftState = heap(leftPosition);
        rightState = heap(rightPosition);
        heap(leftPosition) = rightState;
        heap(rightPosition) = leftState;
        heapPosition(rightState) = uint32(leftPosition);
        heapPosition(leftState) = uint32(rightPosition);
    end

    function siftUp(position)
        while position > 1
            parentPosition = floor(position/2);
            if ~precedes(heap(position),heap(parentPosition))
                break;
            end
            swapHeap(position,parentPosition);
            position = parentPosition;
        end
    end

    function siftDown(position)
        while 2*position <= heapCount
            childPosition = 2*position;
            if childPosition < heapCount && ...
                    precedes(heap(childPosition+1),heap(childPosition))
                childPosition = childPosition + 1;
            end
            if ~precedes(heap(childPosition),heap(position))
                break;
            end
            swapHeap(position,childPosition);
            position = childPosition;
        end
    end

    function poppedIndex = popOpen()
        poppedIndex = heap(1);
        heapPosition(poppedIndex) = uint32(0);
        if heapCount == 1
            heapCount = 0;
            return;
        end
        movedIndex = heap(heapCount);
        heapCount = heapCount - 1;
        heap(1) = movedIndex;
        heapPosition(movedIndex) = uint32(1);
        siftDown(1);
    end
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

function [path1, dirSeq] = reconstructVisitedPath(finalIndex, startP, stateCoords, parentIndex)
% 返回完整网格路径 path1, 以及每步对应的方向序列 dirSeq
    index = finalIndex;

    revPath = zeros(0,3);
    revDir  = zeros(0,1);

    while true
        current = stateCoords(index,:);
        x = current(1);
        y = current(2);
        z = current(3);
        d = current(4);
        revPath(end+1,:) = [x,y,z]; %#ok<AGROW>
        revDir(end+1,1)  = d; %#ok<AGROW>

        if x == startP(1) && y == startP(2) && z == startP(3) && d == 7
            break;
        end

        index = parentIndex(index);
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
