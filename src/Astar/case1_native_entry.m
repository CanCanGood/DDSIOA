function r=case1_native_entry()
% Generate candidate routes for the three Case 1 requests.



rng('default')%%控制随机数的生成器，这里将其设置为初始状态
global memPeak NATIVE_ASTAR_COUNTS
NATIVE_ASTAR_COUNTS=zeros(1,3);
noRng=str2double(getenv('CASE1_TRIAL_SEED'));
tt1=0;
allsearchedVoxels=zeros(0,3);
clear global MULTIRES_CALLS;
global MULTIRES_CALLS; MULTIRES_CALLS={};
global BATCH_TIMINGS; BATCH_TIMINGS={};
selectedSearchPoints={};
selectedSearchCodes={};
interfaceRecords={};physicalRoutes={};
tic;
data.curve=5000;
str_legend=[{'A*'}];
rng(noRng)%%设置随机数生成器的种子为noRng的值，也就是1
memPeak = 0;
global BATCH_RAM_SAMPLES BATCH_RAM_CLOCK BATCH_RAM_EVENTS
BATCH_RAM_SAMPLES=zeros(0,2);BATCH_RAM_EVENTS={};BATCH_RAM_CLOCK=tic;
        m = memory;
        memunit = m.MemUsedMATLAB / 1024^2;
t = timer( ...
    'ExecutionMode','fixedRate', ...
    'Period',0.2, ...
    'TimerFcn',@sampleMem);

start(t); batch_ram_event('start');
batch_read_clock=double(java.lang.System.nanoTime);
data.map1=xlsread(fullfile(fileparts(mfilename('fullpath')),'..','..','data','mapcase1.xlsx'),1);
batch_timing_record('read_excel_input',batch_read_clock);
data.unit=str2double(getenv('MULTIRES_H'));
assert(ismember(data.unit,[1 5 10 15 20 25]),'Controls:NewResolution','Unsupported Astar grid spacing.');
data.unit_xyz=data.unit;
data.mapx1=data.map1(:,1);
data.mapx2=data.map1(:,4);
data.mapy1=data.map1(:,2);
data.mapy2=data.map1(:,5);
data.mapz1=data.map1(:,3);
data.mapz2=data.map1(:,6);
data.obstacles=data.map1(:,7:11);
searchedVoxelsall=0;
data.dmax=15;
data.constraint_policy=case1_constraint_policy();
data.energy_dsafe=data.constraint_policy.energy_dsafe;
data.physical_collision_margin=data.constraint_policy.physical_collision_margin;
data.interface_clearance=data.constraint_policy.interface_clearance;
data.unit_x=1180; %x轴单位长度，其实就是栅格到像素分辨率，设置单个立方体的长度、宽度
data.unit_y=1350; %y轴单位长度
data.unit_z=720; %z轴单位长度，后面栅格化×10层，所以h高50，

data.mapsize=[data.unit_x,data.unit_y,data.unit_z];
data.mapsize1=[data.unit_x,data.unit_y,data.unit_z];
data.mapsize2=[0,0,0];
data.unit_BD=10;%%管线弯曲半径单位长度。
data.w0=0.2;
recording=struct();
    data.pipestyle1=[];
    data.mapx11=[];
    data.mapx21=[];
    data.mapy11=[];
    data.mapy21=[];
    data.mapz11=[];
    data.mapz21=[];
    mapsize_h =[];
batch_grid_alloc_clock=double(java.lang.System.nanoTime);
data.mapsize_h=floor(data.mapsize1/data.unit); %三维地图栅格化 %%%改层高数量
data.map_h=[]; % Unused in the active A* call chain; omit dense allocation.
data.map_er=[]; % First use follows the per-request double field build below.
batch_timing_record('initial_grid_allocation',batch_grid_alloc_clock);

data.numType=size(data.map1,1);
data.dir1=[];
A2=[];
toc;
tt1=tt1+toc;

tic;
disp('初始化起点、终点坐标');
data.coordinates = [
[740,50,529]; [249,600,138];[1,10,3];
[1055,582,529]; [434,1336,638];[1,10,3]; % 第一对起始坐标 热水
[415,136,45];[249,631,23];[6,10,3];
     ];
data.pipe_diameters_input_order=case1_pipe_diameters();
assert(isequal(data.pipe_diameters_input_order,[10 10 10]),'Controls:Physical','This extension preserves all Case1 physical diameters at 10.');
assert(isequal(data.coordinates(3:3:end,2)',data.pipe_diameters_input_order),'Controls:Physical','Physical input diameters changed.');
assert(all(data.coordinates(3:3:end,3)==data.energy_dsafe),'Controls:ConstraintPolicy','Native coordinate metadata must retain the energy dsafe.');
data.numpipe=size(data.coordinates,1)/3;



data.pipe = [[-2,2,0.5,1,5];[-2,2,0.5,6,5];[-2,2,0.5,11,5];];

for mm = 1:3:(size(data.coordinates,1)-3)
    type=data.coordinates(mm+2,2);
    type1=data.coordinates((mm:mm+2),:);
    for nn=mm+3:3:size(data.coordinates,1)
    type2=data.coordinates(nn+2,2);
       type1=data.coordinates((nn:nn+2),:);
    if type2>type
        data.coordinates((nn:nn+2),:)=data.coordinates((mm:mm+2),:);
        data.coordinates((mm:mm+2),:)=type1;   
        type=data.coordinates(mm+2,2);
    end
    end
end

for mm = 1:3:(size(data.coordinates,1)-3)
    type=data.coordinates(mm+2,1);
    type1=data.coordinates((mm:mm+2),:);
    for nn=mm+3:3:size(data.coordinates,1)
    type2=data.coordinates(nn+2,1);
       type1=data.coordinates((nn:nn+2),:);
    if type2>type
        data.coordinates((nn:nn+2),:)=data.coordinates((mm:mm+2),:);
        data.coordinates((mm:mm+2),:)=type1;  
        type=data.coordinates(mm+2,1);
    end
    end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
end

for mm = 1:3:(size(data.coordinates,1)-3)
    type=data.coordinates(mm+2,3);
    type1=data.coordinates((mm:mm+2),:);
    for nn=mm+3:3:size(data.coordinates,1)
    type2=data.coordinates(nn+2,3);
       type1=data.coordinates((nn:nn+2),:);
    if type2<type
        data.coordinates((nn:nn+2),:)=data.coordinates((mm:mm+2),:);
        data.coordinates((mm:mm+2),:)=type1;   
        type=data.coordinates(mm+2,3);
    end
    end
end

for mm = 1:3:(size(data.coordinates,1)-3)
    type=h(data.coordinates(mm,:),data.coordinates(mm+1,:));
    type1=data.coordinates((mm:mm+2),:);
    for nn=mm+3:3:size(data.coordinates,1)
    type2=h(data.coordinates(nn,:),data.coordinates(nn+1,:));
       type1=data.coordinates((nn:nn+2),:);
    if type2>type
        data.coordinates((nn:nn+2),:)=data.coordinates((mm:mm+2),:);
        data.coordinates((mm:mm+2),:)=type1;  
        type=h(data.coordinates(mm,:),data.coordinates(mm+1,:));
    end
    end
end

[markedGroups, filteredGroups] = filterAndMarkGroups(data.coordinates);
data.markedGroups=markedGroups;

Dz=6;
fenzuz=1;
fenzu=struct();
datapipe=struct();
datapipe1=struct();
recording1=struct();





 for mm = 1:size(filteredGroups,3)
   
      typez1=filteredGroups(1,3,mm);
      typez2=filteredGroups(2,3,mm);
    type1=filteredGroups(:,:,mm);
    fenzutemp=[];
     for mm = 1:size(filteredGroups,3)
       typez3=filteredGroups(1,3,mm);
      typez4=filteredGroups(2,3,mm);
    type2=filteredGroups(:,:,mm);
    if chongdie(typez1,typez2,typez3,typez4,Dz)
        fenzutemp=[fenzutemp;type2];
    end
     end
fenzu(fenzuz).fenzutemp=fenzutemp;
     fenzuz=fenzuz+1;
 end
fenzu= find_independent_sets(fenzu);

Dz1=1;
Dz2=0.4;
fenzuz=1;
mm1=1;

fenzutempt=struct();
 for mm1 = 1:(length(fenzu))
   for mm = 1:3:size(fenzu(mm1).fenzutemp,1) 
      typex1=fenzu(mm1).fenzutemp(mm,1);
      typex2=fenzu(mm1).fenzutemp(mm+1,1);
      typey1=fenzu(mm1).fenzutemp(mm,2);
      typey2=fenzu(mm1).fenzutemp(mm+1,2);
      type1=fenzu(mm1).fenzutemp((mm:mm+2),:);
      fenzutemp=[];
     for nn = 1:3:size(fenzu(mm1).fenzutemp,1)
      typex3=fenzu(mm1).fenzutemp(nn,1);
      typex4=fenzu(mm1).fenzutemp(nn+1,1);
      typey3=fenzu(mm1).fenzutemp(nn,2);
      typey4=fenzu(mm1).fenzutemp(nn+1,2);
      type2=fenzu(mm1).fenzutemp((nn:nn+2),:);
    if chonghe(typex1,typex2,typex3,typex4,typey1,typey2,typey3,typey4,Dz1,Dz2)
        fenzutemp=[fenzutemp;type2];
    end
     end

fenzutempt(fenzuz).fenzutemp=fenzutemp;
     fenzuz=fenzuz+1;
     end

 end
fenzutempt= find_independent_sets(fenzutempt);
for i = 1:length(fenzutempt)
    fenzu(i).fenzufinal = fenzutempt(i).fenzutemp;
end
str_legend=[{'A*'}];
aimFcn=[{@GACELL}];
bestY1=zeros(1, length(aimFcn));tt = zeros(1, length(aimFcn));
uu=1;
    path1= struct( ...
    'path', num2cell(zeros(1,length(aimFcn))));
fenzu(1).path=repmat(struct('path', []), 1, size(filteredGroups,3));
toc;
tt1=tt1+toc;
for q=1:length(fenzu)
mm=1;
mm1=1;
path=[];
group_probs =generate_sum_to_one_array(size(fenzu(q).fenzufinal,1)/3);
 selected_done=[];
while mm ~=((size(fenzu(q).fenzufinal,1))/3+1)    
    data.dir=[1,  0, 0; ...
           -1,  0, 0; ...
            0,  1, 0; ...   
            0, -1, 0; ...
            0,  0, 1; ...a
            0,  0,-1; ];
    data.numDir=length(data.dir(:,1));%%获取方向矢量的大小 6
    lb=0;%%设置下界（lower bound）为0。
    ub=1;%%设置上界（upper bound）为1。
    option.lb=lb;
    option.ub=ub;
    option.fobj=@native_search_observe;%%设置目标函数（objective function）为aimFcn_1。这表示优化问题的目标函数是aimFcn_1函数
    option.fobj0=@aimFcn_1;
    option.showIter=0;%%设置生成迭代信息的选项。showIter=0表示不显示优化过程中的每一次迭代信息？!好像没用到

    option.numAgent=6;   
    option.numAgent1=option.numAgent;  
    option.numAgent2=fix(option.numAgent*1.2);  %种群个体数 size of population
    option.numAgent3=option.numAgent;  
    option.numAgent4=option.numAgent;
    option.maxIteration=0; 
    option.maxIteration1=1; 
    option.p1_GA=0.3;  %%交叉概率-疑似
    option.p2_GA=0.6;  %%变异概率-疑似
    option.p3_GA=1; 
    
    option.fir=0.25;
    option.sec=0.5;
    option.thi=0.8;
    option.so=0.5;
    option.co=0.01;
    option.w_pso=0.5;
    option.c1_pso=0.8;
    option.c2_pso=0.8;
    option.betal=-0.3;
    option.betah=0.9;
    option.pcoss1=0.5;
    option.pcoss2=0.5;
    option.nOnlooker=0.3;
    option.a_abc=1;                   % Acceleration Coefficient Upper Bound
    option.temperature = 300;
    disp('初始化种群的位置和适应度');
    rng(noRng)%%设置随机数生成器（RNG）的种子。确保每次运行时生成的随机数都是一样的
    x=cell(1,option.numAgent);
    y=ones(option.numAgent,1);
    group_starts = 1:3:size(fenzu(q).fenzufinal,1);

  if mm==1
    datapipe1.pipestyle1=[];
    datapipe1.mapx11=[];
    datapipe1.mapx21=[];
    datapipe1.mapy11=[];
    datapipe1.mapy21=[];
    datapipe1.mapz11=[];
    datapipe1.mapz21=[];
   datapipe1.dirs=[];
   datapipe1.dire=[];
  end

start_row = group_starts(mm);
selected_data = fenzu(q).fenzufinal(start_row:start_row+2,:);
fourthColumnValues = findMatchingGroups(markedGroups, selected_data);
fenzu(q).fenzutemp=fenzu(q).fenzufinal;
     for qidian=1:length(group_starts)
        data.S= fenzu(q).fenzufinal(group_starts(qidian),:);
    
        radis=fenzu(q).fenzufinal(group_starts(qidian)+2,:);
        data.dsafe=radis(3);
        data.BD=radis(2);
        data.BD1=radis(2);
        data.pipestyle=radis(1);
   nextP=repmat(data.S,data.numDir,1)+data.dir*data.BD1*2;
    nextP1=repmat(data.S,data.numDir,1)+data.dir*data.BD1*2.5;
   for start1 = 1:size(nextP,1)
   if potential_energy(nextP(start1,:),nextP(start1,:), data.mapsize1,data.mapsize2, data, data.dmax,data.BD)~=inf
        fenzu(q).fenzutemp(group_starts(qidian),:) = nextP(start1,:);
        datapipe1.dirs=[ datapipe1.dirs,start1];
    A1=computePipeAABB(nextP1(start1,:),data.S, data.unit);  
    pipe=repmat(data.pipe(fix((radis(1)+2)/5)+1,:),size(A1,1),1);
    datapipe1.pipestyle1=[datapipe1.pipestyle1;pipe];
    datapipe1.mapx11=[datapipe1.mapx11;A1(:,1)];
    datapipe1.mapx21=[datapipe1.mapx21;A1(:,4)];
    datapipe1.mapy11=[datapipe1.mapy11;A1(:,2)];
    datapipe1.mapy21=[datapipe1.mapy21;A1(:,5)];
    datapipe1.mapz11=[datapipe1.mapz11;A1(:,3)];
    datapipe1.mapz21=[datapipe1.mapz21;A1(:,6)];

       break;
   else
       disp("error");
   end
   end
       
   data.E=fenzu(q).fenzufinal(group_starts(qidian)+1,:);
    nextP=repmat(data.E,data.numDir,1)+data.dir*data.BD1*2;
    nextP1=repmat(data.E,data.numDir,1)+data.dir*data.BD1*2.5;
    for start1 = 1:size(nextP,1)
   if potential_energy(nextP(start1,:),nextP(start1,:), data.mapsize1,data.mapsize2, data, data.dmax,data.BD)~=inf
        fenzu(q).fenzutemp(group_starts(qidian)+1,:) = nextP(start1,:);
          datapipe1.dire=[ datapipe1.dire,start1+mod(start1,2)*2-1];
    A1=computePipeAABB(nextP1(start1,:),data.E, data.unit);  
    pipe=repmat(data.pipe(fix((radis(1)+2)/5)+1,:),size(A1,1),1);
    datapipe1.pipestyle1=[datapipe1.pipestyle1;pipe];
    datapipe1.mapx11=[datapipe1.mapx11;A1(:,1)];
    datapipe1.mapx21=[datapipe1.mapx21;A1(:,4)];
    datapipe1.mapy11=[datapipe1.mapy11;A1(:,2)];
    datapipe1.mapy21=[datapipe1.mapy21;A1(:,5)];
    datapipe1.mapz11=[datapipe1.mapz11;A1(:,3)];
    datapipe1.mapz21=[datapipe1.mapz21;A1(:,6)];

       break;
   else
       disp("error");
   end
    end
     end
data.branch=fenzu(q).fenzutemp;     
selected_data =fenzu(q).fenzutemp(start_row:start_row+2,:);
selected_data0=fenzu(q).fenzufinal(start_row:start_row+2,:);
data.fourthColumnValues =fourthColumnValues; 
   data.Start =selected_data(1,:);
   data.S0=selected_data0(1,:);
   data.E =selected_data(2,:);
   data.E0=selected_data0(2,:);
   radis=selected_data(3,:);
      B=size(data.Start,1);
 


if true % One native environment construction per request at every h.
batch_grid_field_clock=double(java.lang.System.nanoTime);
tic
mapsize_h = floor(data.mapsize1 ./ data.unit);   % [Nx, Ny, Nz]
data.map_er = []; % Release the previous request field before allocation.
data.map_er = zeros(mapsize_h);           % 默认都是 0

numObs = size(data.map1,1);
dmax   = data.dmax;
unit   = data.unit;

r_idx = ceil(dmax / unit);   % 在索引空间扩多少格

for i = 1:numObs
    xmin = data.map1(i,1);
    xmax = data.map1(i,4);
    ymin = data.map1(i,2);
    ymax = data.map1(i,5);
    zmin = data.map1(i,3);
    zmax = data.map1(i,6);

    ix0 = max(1, floor(xmin/unit) - r_idx);
    ix1 = min(mapsize_h(1), ceil(xmax/unit) + r_idx);
    iy0 = max(1, floor(ymin/unit) - r_idx);
    iy1 = min(mapsize_h(2), ceil(ymax/unit) + r_idx);
    iz0 = max(1, floor(zmin/unit) - r_idx);
    iz1 = min(mapsize_h(3), ceil(zmax/unit) + r_idx);

    obstaclestlye  = data.obstacles(i,4);
    obstaclestlye1 = judgeobstaclestlye(data.pipestyle, obstaclestlye);
    emin    = data.obstacles(i,1);
    emax    = data.obstacles(i,2);
    emaxmin = data.obstacles(i,3);
    dbest   = data.obstacles(i,5);
    dsafe   = data.dsafe;

    for xg = ix0:ix1
        for yg = iy0:iy1
            for zg = iz0:iz1

                x = xg * unit;
                y = yg * unit;
                z = zg * unit;

                if x >= xmin && x <= xmax && ...
                   y >= ymin && y <= ymax && ...
                   z >= zmin && z <= zmax

                    data.map_er(xg,yg,zg) = 88;
                    continue;
                end

                dx = min(abs([xmin - x, x - xmax]));
                dy = min(abs([ymin - y, y - ymax]));
                dz = min(abs([zmin - z, z - zmax]));
                d1 = [dx, dy, dz];
                positive_values = d1(d1 > 0);
                if isempty(positive_values)
                    continue;
                end
                d = min(positive_values);

                if d > dmax
                    continue;  % 超出影响半径，势能为0，跳过
                end

                U2 = 0;
                if obstaclestlye1 == 1       % 吸引
                    U2 = emin/(dmax-dsafe) * (d - dsafe);

                elseif obstaclestlye1 == 2   % 强排斥
                    U2 = log(d/dmax)/log(0.1) + emax;

                else                         % 其它模式
                    if d > dbest
                        U2 = emaxmin;
                    else
                        U2 = log(d/dbest)/log(0.1) + emaxmin;
                    end
                end

                data.map_er(xg,yg,zg) = max(data.map_er(xg,yg,zg), U2);

            end
        end
    end
end

batch_timing_record('grid_field_build',batch_grid_field_clock);
batch_interface_clock=double(java.lang.System.nanoTime);
[data,interfaceInfo]=repair_grid_interfaces(data,datapipe1,radis(2),data.interface_clearance);
batch_timing_record('grid_interface_repair',batch_interface_clock);
interfaceRecords{end+1}=interfaceInfo;
tt1=tt1+toc;
end

    for i=1:length(aimFcn)%%遍历所有的算法。aimFcn 是一个包含了所有算法函数句柄的向量。aimFcn=[{@GA},{@PSO},{@myDBO},{@myIDBO}];
    data.curve=5000;

    if i==1||i==2
        option.numAgent=option.numAgent1;
    else
        option.numAgent=option.numAgent4;
    end
   x=cell(1,option.numAgent);
    y=ones(option.numAgent,1);
 if mm==1&&q==1
    datapipe(i).pipestyle1=[];
    datapipe(i).mapx11=[];
    datapipe(i).mapx21=[];
    datapipe(i).mapy11=[];
    datapipe(i).mapy21=[];
    datapipe(i).mapz11=[];
    datapipe(i).mapz21=[];
   datapipe(i).dirs=[];
end
   data.Start =data.interface_info.grid_start;
   data.E =data.interface_info.grid_goal;
   radis=selected_data(3,:);

   qidian=1;
   B=size(data.Start,1);
 if mm==1
    datapipe(i).pipestyle1=[datapipe(i).pipestyle1; datapipe1.pipestyle1];
    datapipe(i).mapx11=[datapipe(i).mapx11;datapipe1.mapx11];
    datapipe(i).mapx21=[datapipe(i).mapx21;datapipe1.mapx21];
    datapipe(i).mapy11=[datapipe(i).mapy11;datapipe1.mapy11];
    datapipe(i).mapy21=[datapipe(i).mapy21;datapipe1.mapy21];
    datapipe(i).mapz11=[datapipe(i).mapz11;datapipe1.mapz11];
    datapipe(i).mapz21=[datapipe(i).mapz21;datapipe1.mapz21];
 end


     for qidian=1:B+1
        if (qidian==(B+1))&& (B==size(data.Start,1))
            break;
        end
        data.S= data.Start(qidian,:);
        data.dsafe=radis(3);
        data.BD=radis(2);
        data.BD1=radis(2);
       data.dirSE=[datapipe1.dirs(mm),datapipe1.dire(mm)];
        data.pipestyle=radis(1);
  if qidian==3
        data.dirSE=data.SE3;
  end
    data.mapx11=datapipe(i).mapx11;
    data.mapx21=datapipe(i).mapx21;
    data.mapy11=datapipe(i).mapy11;
    data.mapy21=datapipe(i).mapy21;
    data.mapz11=datapipe(i).mapz11;
    data.mapz21=datapipe(i).mapz21;
    data.pipestyle1=datapipe(i).pipestyle1;
        data.rad=rand()/2;
        data.mm=mm;
        data.curve=5000;
        if mm==1
       if i==3
     for k=1:fix(0.25*option.numAgent)%%循环遍历每一个解
        data.w0=0.55;
        dir=mod(k, 6) + 1;
        data.dir1=data.dir(dir,:);  
            [y(k),x{k}]=option.fobj0(option,data);
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.25*option.numAgent)+1:fix(0.5*option.numAgent)%%循环遍历每一个解
        dir=mod(k, 6) + 1;
        data.dir1=data.dir(dir,:);  
            [y(k),x{k}]=option.fobj0(option,data);
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
       data.dir1=[];
       data.w0=0.2;
     for k=fix(0.5*option.numAgent)+1:fix(0.75*option.numAgent)%%循环遍历每一个解
        
            [y(k),x{k}]=option.fobj0(option,data);
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.75*option.numAgent)+1:option.numAgent%%循环遍历每一个解
            [y(k),x{k}]=option.fobj0(option,data);
     end
          x11=x;
          y11=y;   
       elseif i==6
     for k=1:fix(0.25*option.numAgent)%%循环遍历每一个解
        data.w0=0.55;
        dir=mod(k, 6) + 1;
        data.dir1=data.dir(dir,:);  
        if i==5||i==9||i==10
         [y(k),x{k}]=option.fobj(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit
        else
            [y(k),x{k}]=option.fobj0(option,data);
        end
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.25*option.numAgent)+1:fix(0.5*option.numAgent)%%循环遍历每一个解
        dir=mod(k, 6) + 1;
        data.dir1=data.dir(dir,:);  
        if i==5||i==9||i==10
         [y(k),x{k}]=option.fobj(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit
        else
            [y(k),x{k}]=option.fobj0(option,data);
        end
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
       data.dir1=[];
       data.w0=0.2;
     for k=fix(0.5*option.numAgent)+1:fix(0.75*option.numAgent)%%循环遍历每一个解
        
               if i==5||i==9||i==10
         [y(k),x{k}]=option.fobj(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit
        else
            [y(k),x{k}]=option.fobj0(option,data);
        end
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.75*option.numAgent)+1:option.numAgent%%循环遍历每一个解
        if i==5||i==9||i==10
         [y(k),x{k}]=option.fobj(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit
        else
            [y(k),x{k}]=option.fobj0(option,data);
        end
     end
       elseif i==2
           tic
                  data.dir1=[];
      for k=1:fix(0.25*option.numAgent)%%循环遍历每一个解
        data.w0=0.55;
            [y(k),x{k}]=option.fobj0(option,data);
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.25*option.numAgent)+1:fix(0.5*option.numAgent)%%循环遍历每一个解
            [y(k),x{k}]=option.fobj0(option,data);
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
       data.dir1=[];
       data.w0=0.2;
     for k=fix(0.5*option.numAgent)+1:fix(0.75*option.numAgent)%%循环遍历每一个解
        
            [y(k),x{k}]=option.fobj0(option,data);
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.75*option.numAgent)+1:option.numAgent%%循环遍历每一个解
            [y(k),x{k}]=option.fobj0(option,data);
     end
     x12=x;
     y12=y;
     tt(i)=tt(i)+toc;
       elseif i==1
           tic;
           ymin=Inf; result0=struct(); % Reset candidate minimum; grid-builder ymin is a coordinate.
           data.S=floor(data.S/data.unit);
           data.E=floor(data.E/data.unit);
 data.dir1=[];
      for k=1:fix(0.25*option.numAgent)%%循环遍历每一个解
        data.w0=0.55;
         
                        [y(k),x{k},result]=option.fobj(option,data);
            if y(k)<ymin
                ymin=y(k)
                 bestX{i}=x{k};
               result0=result; 
            end
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.25*option.numAgent)+1:fix(0.5*option.numAgent)%%循环遍历每一个解
           [y(k),x{k},result]=option.fobj(option,data);
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
                     
            if y(k)<ymin
                ymin=y(k)
                 bestX{i}=x{k};
               result0=result; 
            end
     end
       data.dir1=[];
       data.w0=0.2;
     for k=fix(0.5*option.numAgent)+1:fix(0.75*option.numAgent)%%循环遍历每一个解
        
            [y(k),x{k},result]=option.fobj(option,data);
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
                 if y(k)<ymin
                ymin=y(k)
                 bestX{i}=x{k};
               result0=result; 
            end
     end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
     for k=fix(0.75*option.numAgent)+1:option.numAgent%%循环遍历每一个解
            [y(k),x{k},result]=option.fobj(option,data);
                    if y(k)<ymin
                ymin=y(k)
                 bestX{i}=x{k};
               result0=result; 
            end
     end
     tt(i)=tt(i)+toc;
       else
           x=x11;
           y=y11;   
       end
        elseif  i==6
data.s1=[data.S;data.E];
for k=1:fix(0.25*option.numAgent)
      dir=mod(k, 6) + 1;
      data.w0=0.55;
        data.dir1=data.dir(dir,:);  
    data.S=data.s1(1,:);
    data.E=closest_point_on_path(data.S, path1(i).path,data);
    [y1(k),x{k},result1]=option.fobj(option,data);
    data.S=data.s1(2,:);
    data.E=closest_point_on_path(data.S, path1(i).path,data);
    [y2(k),x{k},result2]=option.fobj(option,data);
    data.S=result1.endnode;
    data.E=result2.endnode;
    [y3(k),x{k},result3]=option.fobj(option,data);
    y(k)=y1(k)+y2(k)+y3(k);

          dirpath22=result2.xdir1(end:-2:1);
          dirpath21=result2.xdir1(end-1:-2:1);
          dirpath21=dirpath21+mod(dirpath21,2)*2-1;
          dirpath2=zeros(1, length(dirpath21)*2);
          dirpath2(1:2:end) = dirpath21;
          dirpath2(2:2:end) = dirpath22;
    x{k}=[result1.xdir1,x{k},dirpath2];
 
end
for k=fix(0.25*option.numAgent)+1:fix(0.5*option.numAgent)
      dir=mod(k, 6) + 1;
        data.dir1=data.dir(dir,:);  
    data.S=data.s1(1,:);
    data.E=closest_point_on_path(data.S, path1(i).path,data);
    [y1(k),x{k},result1]=option.fobj(option,data);
    data.S=data.s1(2,:);
    data.E=closest_point_on_path(data.S, path1(i).path,data);
    [y2(k),x{k},result2]=option.fobj(option,data);
    data.S=result2.endnode;
    data.E=result1.endnode;
    [y3(k),x{k},result3]=option.fobj(option,data);
    y(k)=y1(k)+y2(k)+y3(k);

          dirpath22=result2.xdir1(end:-2:1);
          dirpath21=result2.xdir1(end-1:-2:1);
          dirpath21=dirpath21+mod(dirpath21,2)*2-1;
          dirpath2=zeros(1, length(dirpath21)*2);
          dirpath2(1:2:end) = dirpath21;
          dirpath2(2:2:end) = dirpath22;
          
         dirpath32=x{k}(end:-2:1);
         dirpath31=x{k}(end-1:-2:1);
         dirpath31=dirpath31+mod(dirpath31,2)*2-1;
         x{k}(1:2:end) = dirpath31;
         x{k}(2:2:end) = dirpath32;

    x{k}=[result1.xdir1,x{k},dirpath2];
 
end
     data.dir1=[];
     data.w0=0.2;
for k=fix(0.5*option.numAgent)+1:0.75*option.numAgent

        data.S=data.s1(1,:);
    data.E=closest_point_on_path(data.S, path1(i).path,data);
    [y1(k),x{k},result1]=option.fobj(option,data);
    data.S=data.s1(2,:);
    data.E=closest_point_on_path(data.S,path1(i).path,data);
    [y2(k),x{k},result2]=option.fobj(option,data);
    data.S=result1.endnode;
    data.E=result2.endnode;
                  [y3(k),x{k},result3]=option.fobj(option,data);
    y(k)=y1(k)+y2(k)+y3(k);

          dirpath22=result2.xdir1(end:-2:1);
          dirpath21=result2.xdir1(end-1:-2:1);
          dirpath21=dirpath21+mod(dirpath21,2)*2-1;
          dirpath2=zeros(1, length(dirpath21)*2);
          dirpath2(1:2:end) = dirpath21;
          dirpath2(2:2:end) = dirpath22;
    x{k}=[result1.xdir1,x{k},dirpath2];
end
for k=fix(0.75*option.numAgent)+1:option.numAgent
    data.S=data.s1(1,:);
    data.E=closest_point_on_path(data.S, path1(i).path,data);
    [y1(k),x{k},result1]=option.fobj(option,data);
    data.S=data.s1(2,:);
    data.E=closest_point_on_path(data.S, path1(i).path,data);
    [y2(k),x{k},result2]=option.fobj(option,data);
    data.S=result2.endnode;
    data.E=result1.endnode;
    [y3(k),x{k},result3]=option.fobj(option,data);
    y(k)=y1(k)+y2(k)+y3(k);

          dirpath22=result2.xdir1(end:-2:1);
          dirpath21=result2.xdir1(end-1:-2:1);
          dirpath21=dirpath21+mod(dirpath21,2)*2-1;
          dirpath2=zeros(1, length(dirpath21)*2);
          dirpath2(1:2:end) = dirpath21;
          dirpath2(2:2:end) = dirpath22;
          
         dirpath32=x{k}(end:-2:1);
         dirpath31=x{k}(end-1:-2:1);
         dirpath31=dirpath31+mod(dirpath31,2)*2-1;
         x{k}(1:2:end) = dirpath31;
         x{k}(2:2:end) = dirpath32;

    x{k}=[result1.xdir1,x{k},dirpath2];
 
end
    data.S=data.s1(1,:);
    data.E=data.s1(2,:);

        else
      for k=1:fix(0.25*option.numAgent)%%循环遍历每一个解
        dir=mod(k, 6) + 1;
        data.dir1=data.dir(dir,:);  
        data.w0=0.55;
         [y(k),x{k}]=option.fobj0(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit

      end
        data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
       for k=fix(0.25*option.numAgent)+1:fix(0.5*option.numAgent)%%循环遍历每一个解
        dir=mod(k, 6) + 1;
        data.dir1=data.dir(dir,:);  
         [y(k),x{k}]=option.fobj0(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
                 data.dir1=[];
                 data.w0=0.2;
     for k=fix(0.5*option.numAgent)+1:fix(0.75*option.numAgent)%%循环遍历每一个解
         [y(k),x{k}]=option.fobj0(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit
         dirpath22=x{k}(end:-2:1);
         dirpath21=x{k}(end-1:-2:1);
         dirpath21=dirpath21+mod(dirpath21,2)*2-1;
         x{k}(1:2:end) = dirpath21;
         x{k}(2:2:end) = dirpath22;
     end
             data.stemp=data.S;
        data.S=data.E;
        data.E=data.stemp;
          for k=fix(0.75*option.numAgent)+1:option.numAgent%%循环遍历每一个解
         [y(k),x{k}]=option.fobj0(option,data);%%y(i)会接收到函数的第一个返回值，也就是 fit

          end
        end
    assert(isfield(result0,'searchedVoxels'),'Multires:NoValidCandidate','All six attempts failed for this request; network is incomplete.');
    selectedPoints=unique(result0.searchedVoxels,'rows');
    selectedSearchPoints{end+1}=selectedPoints;
    selectedSearchCodes{end+1}=bestX{i};
    physicalRoutes{end+1}=assemble_physical_route(bestX{i},data);
    allsearchedVoxels=[allsearchedVoxels;selectedPoints];
    clear selectedPoints;
    disp('开始算法求解');
    bestY = zeros(length(aimFcn),length(y)); % 初始化最优目标函数值 bestY
        
        tic;
        disp(str_legend{i});
        disp(data.S);
        disp(data.E);
        rng(noRng)%%确保每次运行时所有的算法都是在相同的条件下进行，使结果具有可比性
  
     curvecount=result0.curve_count;
                    energytotal=result0.energy_total;
                    length1=result0.length;
                    palh1=result0.palh1;

        recordingtemp(i).count=curvecount;
        recordingtemp(i).energytotal=energytotal;
        recordingtemp(i).length=length1;
        recordingtemp(i).palh1=palh1;
        
    recordingtemp(i).bestFit=ymin;
    recordingtemp(i).meanFit=ymin;
    recordingtemp(i).wFit=ymin;    
         tt(i)=tt(i)+toc;%%停止计时，然后将运行时间保存到 tt 向量的第 i 个位置
    bestY1(i)=bestY1(i)+bestY(i);
    if q==1&&mm==1&&i==1
        recording=recordingtemp;
        recording(i).bestFit=recording(i).bestFit-800000;
        recording(i).meanFit=recording(i).meanFit-800000;
         recording(i).wFit=recording(i).wFit-800000;
         recording(i).x_g = {bestX{i}};
         recording(i).x_W = {bestX{i}};
        recording(i).count=recording(i).count;
        recording(i).energytotal=recording(i).energytotal;
        recording(i).length=recording(i).length+data.BD1*2;
        recording(i).palh1=recording(i).palh1;
    elseif q==1&&mm==1 
        recording(i)=recordingtemp(i);
        recording(i).bestFit=recording(i).bestFit-800000;
        recording(i).meanFit=recording(i).meanFit-800000;
         recording(i).wFit=recording(i).wFit-800000;
         recording(i).x_g = {bestX{i}};
         recording(i).x_W = {recording(i).x_W};
        recording(i).count=recording(i).count;
        recording(i).energytotal=recording(i).energytotal;
        recording(i).length=recording(i).length+data.BD1*2;
        recording(i).palh1=recording(i).palh1;
    else
    recording(i).bestFit=recordingtemp(i).bestFit+recording(i).bestFit-800000;
    recording(i).meanFit=recordingtemp(i).meanFit+recording(i).meanFit-800000;
    recording(i).wFit=recordingtemp(i).wFit+recording(i).wFit-800000;
        recording(i).count=recordingtemp(i).count+recording(i).count;
        recording(i).energytotal=recordingtemp(i).energytotal+recording(i).energytotal;
        recording(i).length=recordingtemp(i).length+recording(i).length+data.BD1*2;
    if ~isfield(recording(i), 'x_g') || isempty(recording(i).x_g)
        recording(i).x_g = {};
    end
    recording(i).x_g{end+1} = bestX{i};

    if ~isfield(recording(i), 'x_W') || isempty(recording(i).x_W)
        recording(i).x_W = {};
    end
    recording(i).x_W{end+1} = {recording(i).x_W};
    end
    if i==1
    [ybest,~,result0(i)]=option.fobj(option,data,bestX{i});
    
    else
          [ybest,~,result0(i)]=option.fobj0(option,data,bestX{i});
    end
    recording1(uu).curve=result0(i).curve;
    searchedVoxels=result0(i).searchedVoxels;
    searchedVoxelsall=searchedVoxelsall+size(searchedVoxels,1);
    uu=uu+1;
      if mm1==1
          if i==2
        fenzu(1).path(fourthColumnValues).path=result0(i).path;
        A=result0(i).path;
        aabb = computePipeAABB(fenzu(q).fenzufinal(start_row,:),data.S, data.unit);
        A1=[aabb;result0(i).path1];
        aabb = computePipeAABB(fenzu(q).fenzufinal(start_row+1,:),data.E, data.unit);
        A1=[A1;aabb];
       path1(i).path=result0(i).dirpath;  
       searchedVoxels=result0(i).searchedVoxels;
       else
               fenzu(1).path(fourthColumnValues).path=result0(i).path;
        A=result0(i).path;
        aabb = computePipeAABB(fenzu(q).fenzufinal(start_row,:),data.S*data.unit,data.unit);
        A1=[aabb;result0(i).path1];
        aabb = computePipeAABB(fenzu(q).fenzufinal(start_row+1,:),data.E*data.unit, data.unit);
        A1=[A1;aabb];
       path1(i).path=result0(i).dirpath;  
          end
      else
       
      if qidian==1
       data.SE3=[];
       A=result0(i).path(1:end-1,:);
       aabb1 = computePipeAABB(data.Start(qidian,:),data.S, data.unit);
       A1=[aabb1;result0(i).path1];
       fenzu(1).path(fourthColumnValues).path=[fenzu(1).path(fourthColumnValues).path;A];
       dirpath1=result0(i).xdir;
       data.SE3=[result0(i).dirpath(end-1)];
       data.SE4=data.dirSE(1);

                 
       
       ybest1=ybest; 
      end

      if qidian==2
          A=result0(i).path(1:end-1,:);
          aabb2 = computePipeAABB(data.Start(qidian,:),data.S, data.unit);
       A1=[aabb2;result0(i).path1];
          A=A(end:-1:1, :);  
          data.SE3=[data.SE3;result0(i).dirpath(end-1)+mod(result0(i).dirpath(end-1),2)*2-1];
          dirpath22=result0(i).xdir(end:-2:1);
          dirpath21=result0(i).xdir(end-1:-2:1);
          dirpath21=dirpath21+mod(dirpath21,2)*2-1;
          dirpath2=zeros(1, length(dirpath21)*2);
          dirpath2(1:2:end) = dirpath21;
          dirpath2(2:2:end) = dirpath22;
          if i==2||i==3||i==5
           [data.E,dirpath2]= end_point_on_path(data.E,dirpath2,data);
          end
          ybest2=ybest;
          data.E1=data.S;
          data.SE4=[data.SE4;data.dirSE(2)+mod(data.dirSE(2),2)*2-1];
      end
      if qidian==3
          fenzu(1).path(fourthColumnValues).path=[fenzu(1).path(fourthColumnValues).path;result0(i).path(1:end-1,:);A;];
          A=result0(i).path(1:end,:);
          A1=[result0(i).path1];
          ybest3=ybest+recordingtemp(i).y_p+ybest1;
          dirpath11=repmat({dirpath1},option.numAgent,1);
          dirpath21=repmat({dirpath2},option.numAgent,1); 
          dirpath3=cellfun(@(a,b,c) [a, b, c], dirpath11,recordingtemp(i).x_p', dirpath21, 'UniformOutput', false);
        xco=dirpath3;
        yco=ybest3;



      end

      end
        if mm1==1||qidian==3
            stop(t); batch_ram_event('stop');
    disp('画多线ing');

    searchedVoxels = (searchedVoxels - 1) * data.unit;
    pipe=repmat(data.pipe(fix((radis(1)+2)/5)+1,:),size(A1,1),1);
    datapipe(i).pipestyle1=[datapipe(i).pipestyle1;pipe];
    datapipe(i).mapx11=[datapipe(i).mapx11;A1(:,1)];
    datapipe(i).mapx21=[datapipe(i).mapx21;A1(:,4)];
    datapipe(i).mapy11=[datapipe(i).mapy11;A1(:,2)];
    datapipe(i).mapy21=[datapipe(i).mapy21;A1(:,5)];
    datapipe(i).mapz11=[datapipe(i).mapz11;A1(:,3)];
    datapipe(i).mapz21=[datapipe(i).mapz21;A1(:,6)];
    
    start(t); batch_ram_event('start');
        end
        A2=[A2;A1];
      end
    end
    mm=mm+1;
    mm1=1;

    end
end



    for i=1:length(aimFcn)
    end


    for i=1:length(aimFcn)
    end


    for i=1:length(aimFcn)
    end

    for i=1:length(aimFcn)
    end
    stop(t); batch_ram_event('stop'); delete(t);
peakMB = memPeak-memunit;
fprintf('Peak memory: %.2f MB\n', peakMB);


r=struct();
r.calls=MULTIRES_CALLS;r.physical_routes=physicalRoutes;r.interfaces=interfaceRecords;
r.selected_search_points=selectedSearchPoints;r.selected_search_codes=selectedSearchCodes;
r.recording=recording;r.native_search_time_s=tt;r.native_environment_time_s=tt1;
r.sampled_memory_increment_MiB=peakMB;r.timing_segments=BATCH_TIMINGS;
end
    function sampleMem(~,~)
    global memPeak
        m = memory;
        memMB = m.MemUsedMATLAB / 1024^2;
        memPeak = max(memPeak, memMB);
        global BATCH_RAM_SAMPLES BATCH_RAM_CLOCK
        BATCH_RAM_SAMPLES(end+1,:)=[toc(BATCH_RAM_CLOCK),memMB];
    end

