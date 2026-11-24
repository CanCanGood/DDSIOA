%%是优化问题的目标函数，aimFcn_1，即option.fobj，主要用于实现路径查找或路径优化的功能，它具体实现了一个搜索算法，这个算法会根据当前位置和可能的下一步位置来确定其路径
function [fit,x,result]=aimFcncell(option,data,x3,x1,x2)%%该函数定义有三个输入参数，xx为搜索空间，option为搜索选项，data为携带的数据信息。
global CASE1_LAST_SEARCH_COUNTS
CASE1_LAST_SEARCH_COUNTS=[double(nargin~=3),0,0]; % emitted points including start, main loop visits, raw proposals
%%
% global data %% 将数据设置为全局可用
markedGroups=data.markedGroups;
if nargin==3
    
    sP=data.S;
    searchedVoxels=[sP];             %% 当前搜索的起始位置
   endP=data.E;
     path(1,:)=data.S;
   
   for i=0:2:size(x3,2)-1
       for j=1:x3(i+2)
       path(end+1,:)= path(end,:)+data.dir(x3(i+1),:);
       end
   end    
if ~contains_arraycell(endP, path(end,:))
    path(end,:) = data.E;
end
endP=path(end,:);

else 
    
    if nargin==5
        sP=data.S;
        searchedVoxels=[sP];
        endP=data.S;
   for i=0:2:size(x1,2)-1
    sP=sP+data.dir(x1(i+1),:)*x1(i+2);
   end
   for i=0:2:size(x2,2)-1
    endP=endP+data.dir(x2(i+1),:)*x2(i+2);
   end
    else
sP=data.S;
searchedVoxels=[sP];%% 当前搜索的起始位置
endP=data.E;
    end



nowP=sP; 

% flag=data.map0; %% 创建一个地图（复制原始10层栅格地图），用于记录已经走过的位置。10层栅格地图+走过的路径点成为新障碍
flag=data.map_er; %% 创建一个地图（复制原始10层栅格地图），用于记录已经走过的位置。10层栅格地图+走过的路径点成为新障碍
flag(nowP(1),nowP(2),nowP(3))=88;%% 将起始位置标记为已经走过？？等于1是否可以？！改为其他值
path=[nowP]; %%% 初始化路径，先将起始位置加入路径
if ~isempty(data.dir1)
dir=data.dir1;
else
    dir1=randi(6);    
dir=data.dir(dir1,:);

end
searchedVoxels=[sP];
% disp("目标函数调用")
% disp(data.S);
% disp(data.E);
%% 下面的while循环是为了在搜索到终点之前让搜索过程持续进行
%%整段while的作用是为了在6个方向的节点中删除无效的，找到最好的那个
while 1
CASE1_LAST_SEARCH_COUNTS(2)=CASE1_LAST_SEARCH_COUNTS(2)+1;
    %%计算nowP所有可能的6个方向下一步的位置。data.numDir移动的方向个数6。data.dir矩阵，内容是每一个可能的移动方向对应的位移
    %%repmat(nowP,data.numDir,1)这段代码会将当前位置nowP复制为data.numDir行   的矩阵，对应了我们要向各个方向移动的基础位置
    alljp=[];
    jp_list=[];
    nextP=repmat(nowP,data.numDir,1)+data.dir;
    %%nextP是一个6*3矩阵，记录了所有可能的下一步位置
    ind=[]; %% 初始化无效步骤的索引
CASE1_LAST_SEARCH_COUNTS(3)=CASE1_LAST_SEARCH_COUNTS(3)+size(nextP,1);
    CASE1_LAST_SEARCH_COUNTS(1)=CASE1_LAST_SEARCH_COUNTS(1)+size(nextP,1);
searchedVoxels=[searchedVoxels;nextP];
    %%处理越界和已经走过的位置
    for noDir=1:size(nextP,1) %% 遍历所有可能的移动或者说搜索方向
        %%判断按此方向搜索得到的下一个点是否越界，如果越界则添加其方向索引到ind数组中
        if nextP(noDir,1)<=0 || nextP(noDir,2)<=0 || nextP(noDir,3)<=0 || nextP(noDir,1)>data.mapsize_h(1) || nextP(noDir,2)>data.mapsize_h(2) || nextP(noDir,3)>data.mapsize_h(3)
            ind=[ind;noDir];%%在ind变量（是一个列 向量）的末尾添加变量noDir的值，即把这个方向的编号noDir添加到ind列表的末尾
        else
            % 判断按此方向搜索得到的下一个点是否已被访问过（也就是在地图上是否被标记过），如果是则添加其方向索引到ind数组中
            % 其中，flag是一个三维数组，用来记录哪些点是已经走过的或者是障碍物。如果某一个点走过了或者是障碍物，就会在这个位置标记1，否则就是0！改为其他不可穿障碍134与管路径点
              if flag(nextP(noDir,1),nextP(noDir,2),nextP(noDir,3))==88
                  if contains_arraycell(endP, [nextP(noDir,1),nextP(noDir,2),nextP(noDir,3)])
                      nowP=[nextP(noDir,1),nextP(noDir,2),nextP(noDir,3)];
                      path=[path;nowP];
                  break;%%结束for循环
                  end
                 ind=[ind;noDir];
             end
            % % 当flag不为{0，1，10}中的某个值时，执行后面的操作{自由，墙，通行虚拟障碍}
            % if ~ismember(flag(nextP(noDir,1),nextP(noDir,2),nextP(noDir,3)), 0)
            %      ind=[ind;noDir];
            % end
        end
    end
    if contains_arraycell(endP, nowP)
        endP= nowP;
        break;%%结束while循环
    end       


    nextP(ind,:)=[]; %% 删除无效步骤 6*3矩阵-n*3矩阵
    %%处理路径搜索过程中遇到死路（无法继续前进）的情况：转向或后撤
    if isempty(nextP) %% 判断下一步的位置（nextP）是否为空，也就是是否有可走的路径。如果为空，说明当前位置无法再向前，即遇到了死路
%         disp(nowP);
        if size(path,1)<2, error('GBH:Exhausted','Search exhausted at the start node.'); end
        nowP=path(end-1,:); %% 从路径中取出倒数第二个位置，也就是从当前位置的前一个位置，将其设置为当前位置，用于回溯
        path(end,:)=[]; %% 从路径中删除当前的位置，也就是路径的最后一个位置
        continue; %% 跳过本次循环的后续操作，直接开始下一轮循环。也就是在回溯到上一位置后，跳过后续操作，再次尝试寻找下一步的路径。
    end
    %% 初始化 优先级得分列表 距离得分列表 方向得分列表
    pri=[];%% 初始化一个空的优先级得分列表 n*1矩阵-xx
    D=[]; %% 初始化一个空的距离得分列表 n*1矩阵
    Dr=[]; %% 初始化一个空的方向得分列表 n*3矩阵
     Dh=[];
      Dhh=0;
      energy=[];
    %% 计算优先级，距离和方向
    for i=1:length(nextP(:,1)) %% 遍历每一个可能的下一步位置
        x=nextP(i,1); %% 获取第i个下一步位置的x坐标
        y=nextP(i,2); %% 获取第i个下一步位置的y坐标
        z=nextP(i,3); %% 获取第i个下一步位置的z坐标
        
        % D(i)=max(abs(nowP(i,:)-nextP(i,:)))+max(abs(nextP(i,:)-data.E));
        % 计算第i个下一步位置的方向得分，计算当前位置到下一步位置的方向向量与已知方向向量dir的距离作为评分。 
        v=[nextP(i,:)-nowP]/sum(abs([nextP(i,:)-nowP]));
        
        D(i)=h(nextP(i,:),endP)-h(nowP,endP);
        %%计算从当前位置到候选下一步位置的方向与某个参考移动方向之间的差距
        %%使用norm函数计算向量v与参考方向dir之间的欧氏距离，然后将结果赋值给Dr(i)
        Dr(i)=norm(v-dir);
    for mm=1:3:size(markedGroups,1)
        if isequal(markedGroups(mm,4), data.fourthColumnValues)
          
           
         
            Dhh=Dhh+h(nextP(i,:),markedGroups(mm,1:3));
        
           Dhh=Dhh+h(nextP(i,:),markedGroups(mm+1,1:3));
        

        end
         Dh(i)=Dhh;

    end  
    energy(i)=data.map_er(nextP(i,1), nextP(i,2),nextP(i,3));
    end
    %% 归一化处理
    %%这段代码的作用是对三个变量（pri、D、Dr）进行归一化，然后把他们相加，并找出加和后的结果中最小值对应的索引
    temp1=myMapminmax(Dh);
    temp2=myMapminmax(D);
    temp3=myMapminmax(Dr);
    temp4=myMapminmax(energy);

    temp=0.01*temp1+0.87*temp2+0.11*temp3+0.01*temp4; %%计算以上三个归一化操作后的结果的和！？改系数？
    %%寻找temp中值等于其最小值的所有元素的索引，这个索引被存放在no变量中。
    [no]=find(temp==min(temp));%% find函数会返回所有满足条件的值的索引，返回值no是一个数组
    no=no(1); %% 只取索引数组中的第一个元素，如果有多个元素满足条件，我们只考虑第一个。

    %% 更新位置和路径
    %%这段代码的主要作用是路径规划，通过选取最优点并记录路径，实现对物体从一个位置移动到另一个位置的过程的模拟
    nextP=nextP(no,:);%% 在候选点集合nextP中选取索引为no的行作为下一步位置。
    if isempty(data.dir1)
    dir=nextP-nowP; %%计算参考方向，即当前位置到下一个位置的向量。
    end
    nowP=nextP; %%更新当前位置为下一步的位置
    %%在flag数组中，将新的当前位置（nowP）对应的元素值设置为1。数组flag用于记录已经走过的位置。
    flag(nowP(1),nowP(2),nowP(3))=88;
    path=[path;nowP];% 将新的当前位置添加到路径path的末尾。路径path用于记录移动过程中经过的所有位置。

    %% 检查是否已经找到结束位置
    %%ismember 函数用于判断 'nowP' 是否在 'data.E' 的行中存在。
    %%如果当前位置 'nowP' 已经到达结束位置（即 'nowP' 在 'data.E' 中），ismember 函数返回真 (1)。
    %%这里 "rows" 是 ismember 函数的参数，指定操作的维度为行。
    
    if contains_arraycell(endP, nowP)
        endP=nowP;
        break;%%结束while循环
    end
   
    
    %% 暂时注释掉了用于可视化搜索过程的代码
     % hold on
     % plot3(path(:,1)*data.unit_x,path(:,2)*data.unit_y,path(:,3)*data.unit_z,'-k');
     % drawnow
end
end
%% 计算路径长度
%%路径的优度是以路径的总长度（这里是所有连续两点间的欧氏距离之和）
path1=[];
fit=0;
curve_count = 0; %% 初始化弯曲次数
penetration_count = 0; %% 初始化穿墙次数
energy_total= 2000000; %%初始化能量总和
curve(1,:)=sP; 
curve_len=0;
palh1=0;
xdir1=[1,0];
xdir=[1,0];
if size(path,1)==1
    x=[];
   return;
end

x(1)=find(ismember(data.dir, path(2,:)-path(1,:), 'rows'));
if x(1)~=data.dirSE(1)
 curve_count = curve_count + 1;
end
ii=2;
for i=1:length(path(:,1))-1
    %% 计算路径中每两个连续点的欧氏距离，并将它加到总路径长度之上。
    %%在这里，norm(path(i,:)-path(i+1,:))表示的是点path(i,:)和点path(i+1,:)的欧氏距离。
    fit=fit+norm(path(i,:)-path(i+1,:));
    %% 判断是否发生弯曲并累计弯曲次数
    if i>1 && any(path(i-1,:)-path(i,:) ~= path(i,:)-path(i+1,:))
        curve(ii,:)=path(i,:);
        x((ii-1)*2)=sum(abs(curve(ii,:)-curve(ii-1,:)));
        x(ii*2-1)=find(ismember(data.dir, path(i+1,:)-path(i,:), 'rows'));
        curve_count = curve_count + 1;
        curve_distance=h(curve(ii,:),curve(ii-1,:));
        aabb = computePipeAABB(curve(ii-1,:)*data.unit,curve(ii,:)*data.unit, data.BD1);
        path1=[path1;aabb];
        ii=ii+1;
        if curve_distance<2* data.BD
        curve_len=10000000;
        end
    end
    
    %% 判断是否穿墙并累计穿墙次数
    % if data.map_h(path(i,1), path(i,2), path(i,3)) == 1
    %     penetration_count = penetration_count +1;
    % end
    %% 累计节点能量值
    energy_total = energy_total + data.map_er(path(i,1), path(i,2), path(i,3));
end
dis_len = fit * data.unit+ pi / 2 * data.BD1 * curve_count * data.unit -data.BD1 * curve_count * data.unit;
curve_len =(3^curve_count)*data.curve+curve_len;
penetration_len = penetration_count * data.unit_xyz;
if ~isequal(curve(end,:),endP)
aabb = computePipeAABB(curve(end,:)*data.unit,endP*data.unit, data.BD1);
path1=[path1;aabb];
end
%计算总能量值
% energy_total = energy_total;
%计算综合函数，注意权重之和为1，这里暂定各项权重为0.25
fit = 0.6 * dis_len + 0.4 * curve_len + 0.5 * (penetration_len) + 0.4*energy_total;
x(end+1)=sum(abs(endP-curve(end,:)));
energy_total=energy_total- 2000000;
% 下面的if语句是创建一个包含路径长度和路径信息的结果对象，如果函数有第二个输出参数的话
%%nargout 是一个内置变量，表示函数的输出参数的数量。当函数需要返回的参数数量大于 1 时，执行以下操作。
if nargout>1
    result.fit=fit;
    result.length=dis_len;
    result.path=curve;
    result.path1=path1;
    result.path2=xdir;
    result.endnode=sP;
    result.dirpath=[sP,xdir];
    result.xdir=xdir;
    result.xdir1=xdir1;
    result.dis_len = dis_len;
    result.curve_len = curve_len;
    result.curve = curve;
    result.penetration_len = penetration_len;
    result.energy_total = energy_total;
    result.curve_count=curve_count;
     result.palh1=palh1;
result.searchedVoxels =searchedVoxels;
end
end
