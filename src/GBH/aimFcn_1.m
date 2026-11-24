%%是优化问题的目标函数，aimFcn_1，即option.fobj，主要用于实现路径查找或路径优化的功能，它具体实现了一个搜索算法，这个算法会根据当前位置和可能的下一步位置来确定其路径
function [fit,xdir,result,fit1]=aimFcn_1(option,data,x3,x1,x2)%%该函数定义有三个输入参数，xx为搜索空间，option为搜索选项，data为携带的数据信息。
%%
% global data %% 将数据设置为全局可用
step=1;
dir=[];

data.mapsize1(1)=data.mapsize(1)-data.BD1;
data.mapsize1(2)=data.mapsize(2)-data.BD1;
data.mapsize1(3)=data.mapsize(3)-data.BD1;
data.mapsize2=[data.BD1,data.BD1,data.BD1];
markedGroups=data.markedGroups;
if nargin==3
    sP=data.S; 
    nowP1=sP;%% 当前搜索的起始位置
   dirpath1=[sP,x3];
   % searchedVoxels=[sP];  
   dirpath=[sP,x3(1:2)];
   nowP1=nowP1+data.dir(x3(1),:)*x3(2);
   endP=data.E;
   if length(dirpath1)>6
   for i=6:2:length(dirpath1)
     if dirpath1(i)==dirpath(end-1)
                dirpath(end)=dirpath(end)+dirpath1(i+1);
     elseif fix((dirpath1(i)+1)/2)==fix((dirpath(end-1)+1)/2)
                dirpath(end)=dirpath(end)-dirpath1(i+1);    
     elseif dirpath1(i+1)~=0
                dirpath=[dirpath,dirpath1(i),dirpath1(i+1)];
     end

     if sign(dirpath(end))<0
         dirpath(end-1)=dirpath(end-1)+mod(dirpath(end-1),2)*2-1;
         dirpath(end)=abs(dirpath(end));

     end
     nowP1=nowP1+data.dir(dirpath1(i),:)*dirpath1(i+1);
   end
   if  fix((dirpath(4)+1)/2)~=fix((data.dirSE(1)+1)/2)
         dirpath=[dirpath(1:3),data.dirSE(1),0,dirpath(4:end)];
   end
   
   if  fix((dirpath(end-1)+1)/2)~=fix((data.dirSE(2)+1)/2)
         dirpath=[dirpath,data.dirSE(2),0];
   end
   end

endnode=sP;

%% 计算路径长度
%%路径的优度是以路径的总长度（这里是所有连续两点间的欧氏距离之和）
xdir=dirpath(4:end);
xdir1=[dirpath(4),0];
fit=0;
u=0;
curve_count = 0; %% 初始化弯曲次数
penetration_count = 0; %% 初始化穿墙次数
energy_total= 0; %%初始化能量总和
curve(1,:)=sP; 
ii=2;
 path1=[];
 nowP=dirpath(1:3);
 curve_count = (length(dirpath)-5)/2;
tempstep=12;
curve_len=0;
 palh1=0;


for i=4:2:length(dirpath)-1
    %% 计算路径中每两个连续点的欧氏距离，并将它加到总路径长度之上。
    %%在这里，norm(path(i,:)-path(i+1,:))表示的是点path(i,:)和点path(i+1,:)的欧氏距离。
    temp=dirpath(i:i+1);
    fit=fit+temp(2);
    % if i~=1&&temp(2)<2*data.BD1&&i~=length(dirpath)-1
    %     curve_len=curve_len+10000000;
    % end
    tempP1=curve(end,:);
    tempP2=tempP1+data.dir(temp(1),:)*tempstep;
    ij=tempstep;
    %% 累计节点能量值固定步长
    while ij<temp(2)
        if temp(2)<tempstep
            break;
        end
    
    [u,~]=potential_energy(tempP1,tempP2, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
     energy_total = energy_total + u;
    [u,~,np]= checkpipe(tempP1,tempP2, data, data.dmax,data.BD1);
            if u>8
    u=u*1000000000;
            else
            endnode=tempP2;
            xdir1=dirpath(4:i+1);
            xdir1(end)=ij;
             
            end
         
            if u<0
palh1=palh1+tempstep;
            end
            
    energy_total = energy_total + u;
    tempP1=tempP2;
    tempP2=tempP1+data.dir(temp(1),:)*tempstep;
    ij=ij+tempstep;
    end
    %% 判断是否发生弯曲并累计弯曲次数
    curve(ii,:)=curve(end,:)+data.dir(temp(1),:)*temp(2);
        aabb = computePipeAABB(curve(ii-1,:),curve(ii,:), data.BD);
        path1=[path1;aabb];
        ii=ii+1;  
        if u<8 && temp(2)~=0
        [u,~,np]= checkpipe(tempP1,curve(end,:), data, data.dmax,data.BD1);
        if u>8
           u=u*1000000000;
        else
            endnode=curve(end,:);
            xdir1=dirpath(4:i+1);
         
        end
%                     if u<0
% palh1=palh1+tempstep;
%             end
        else
        [u,~,np]= checkpipe(tempP1,curve(end,:), data, data.dmax,data.BD1);
        if u<8
            endnode=curve(end,:);
            xdir1=dirpath(4:i+1);
         
        end  
%                     if u<0
% palh1=palh1+tempstep;
%             end
        end
        energy_total = energy_total + u;
        [u,~]=potential_energy(tempP1,curve(end,:), data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
         if u==inf
            energy_total=inf;
         end 
    energy_total = energy_total + u;

        

    %% 判断是否穿墙并累计穿墙次数

        penetration_count =0;
end


%根据要求计算各部分
dis_len = fit+ pi / 2 * data.BD * curve_count -data.BD * curve_count;
curve_len =(3^curve_count)*data.curve+curve_len;
curve_len1 =(3^curve_count)*5000+curve_len;
penetration_len = penetration_count * data.unit_xyz;
palh=0;
palh0=0;
 for mm=1:3:size(data.branch,1)
   if~(isequal(data.S, data.branch(mm,1:3))&&isequal(data.E, data.branch(mm+1,1:3)))||(isequal(data.S,data.branch(mm+1,1:3))&&isequal(data.E,data.branch(mm,1:3)))

p1=closest_point_on_path(data.branch(mm,1:3),{dirpath},data);
p2=closest_point_on_path(data.branch(mm+1,1:3),{dirpath},data);
palh=palh-h(p1,data.branch(mm,1:3))-h(data.branch(mm+1,1:3),p2);
palh0=palh-h(p1,data.branch(mm,1:3))-h(data.branch(mm+1,1:3),p2);
   end
 end 
 if size(data.branch,1)>3&&(isequal(data.S, sP)&&isequal(data.E,endP))
 palh=palh/(size(data.branch,1)-3)+palh1/2*data.plw;
 palh0=palh/(size(data.branch,1)-3)*3+palh1/2;
 end

% 计算总能量值
energy_total = energy_total;
%计算综合函数，注意权重之和为1，这里暂定各项权重为0.25
fit =  0.6 * dis_len + 0.4 * curve_len + 0.5 * (penetration_len-palh) + 0.4*energy_total;
fit1 =  0.6 * dis_len + 0.4 * curve_len1 + 0.5 * (penetration_len-palh0) + 0.4*energy_total;
if curve(end,:)~=endP
aabb = computePipeAABB(curve(end,:),endP, data.BD);
path1=[path1;aabb];
end
aabb = computePipeAABB(curve(end,:),endP-data.dir(dirpath(end-1),:)*data.BD*1.2, data.BD);
path2=[path1;aabb];
% 下面的if语句是创建一个包含路径长度和路径信息的结果对象，如果函数有第二个输出参数的话
%%nargout 是一个内置变量，表示函数的输出参数的数量。当函数需要返回的参数数量大于 1 时，执行以下操作。
if nargout>1

    result.fit=fit;
    result.length=dis_len;
    result.path=curve;
    result.path1=path1;
    result.path2=path2;
    result.endnode=endnode;
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
     % result.searchedVoxels =searchedVoxels;
end
return;

else 

    if nargin==5
        sP=data.S;
            % searchedVoxels=[sP];  
        endP=data.S;
   for i=0:2:length(x1)-1
    sP=sP+data.dir(x1(i+1),:)*x1(i+2);
   end
   dir=data.dir(x1(end-1),:);
   for i=0:2:length(x2)-1
    endP=endP+data.dir(x2(i+1),:)*x2(i+2);
   end
if potential_energy(endP,endP, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1)==inf
    endP=data.E;
end


   
    else
sP=data.S;     
    % searchedVoxels=[sP];  %% 当前搜索的起始位置
endP=data.E;
    end



nowP=sP; 
endnode=sP;
flag=0;
flagrec=[];%% 将起始位置标记为已经走过？？等于1是否可以？！改为其他值
if isequal(sP,data.S)
len=5;
dirpath=[nowP,data.dirSE(1),0];
v=data.dirSE(1);
else
  
len=0;
v=randi(6);
dirpath=[nowP,v,0];
end
path=nowP; %%% 初始化路径，先将起始位置加入路径
if ~isempty(data.dir1)
dir=data.dir1;
elseif isempty(dir)
dir1=randi(6);  
dir=data.dir(dir1,:);
end

circlenumber=1;
circledir1=[];
flagpipe=0;
pathcir=[nowP];
% disp("目标函数调用")
% disp(data.S);
% disp(data.E);
%% 下面的while循环是为了在搜索到终点之前让搜索过程持续进行
%%整段while的作用是为了在6个方向的节点中删除无效的，找到最好的那个
while 1
junpflag=0;
%%目标锁定阶段 
if all(abs(nowP-endP)<=data.dmax+data.BD1/2)

    detaP=endP-nowP;
    count = find(detaP);
    nextP=repmat(nowP,length(count),1);
    if length(count)==2&&isequal(data.E,endP)
       dir=data.dir(data.dirSE(2)+mod(data.dirSE(2),2)*2-1,:);   
    end
    for counti=1:length(count)
    nextP(counti,count(counti))=endP(count(counti));
    end
  


else


 %%循环排除阶段
if (step<1&&(potential_energy(nowP,nowP+data.dir(v,:)*data.BD1, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1)==inf)&&boundcheck(nowP,nowP+data.dir(v,:)*data.BD1, data.mapsize1,data.mapsize2)~=inf)

 step=minobstacle(nowP,nowP+data.dir(v,:)*data.BD1,data,data.dmax,data.BD1)+0.0001 ;  

 nextP=repmat(nowP,data.numDir,1)+data.dir.*step;
 ob= checkobstacle(nowP+data.dir(v,:)*data.BD1,nowP-data.dir(dirpath(end-1),:)*dirpath(end),data,data.dmax,data.BD1);
 nextP(fix((v-1)/2)*2+1:fix((v-1)/2)*2+2,:)=[];
junpflag=1;
 % flagrec=[flagrec;ob,v];
  flagrec=[ob,v];
elseif  step<1&&boundcheck(nowP,nowP+data.dir(v,:)*data.BD1,data.mapsize1,data.mapsize2)~=inf
 step=data.BD1; 

 nextP=repmat(nowP,data.numDir,1)+rand(1)*data.dir.*step;
 nextP(v,:)=nowP+data.dir(v,:)*data.BD1;

elseif length(find(abs(nowP-endP)<=data.dmax))==2
    detaP=endP-nowP;
    count = find(abs(nowP-endP)>data.dmax);
    vi=[0,0,0];
    vi(count)=sign(detaP(count));
    if potential_energy(endP,endP, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1,flagrec)==inf
    flagrec=[];
    end
    if potential_energy(nowP,nowP+vi*data.dmax, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1,flagrec)==inf&&boundcheck(nowP,nowP+vi*data.dmax,data.mapsize1,data.mapsize2)~=inf
     step=minobstacle(nowP,nowP+vi*data.dmax,data,data.dmax,data.BD1)+0.0001; 
      junpflag=1;
     nextP=repmat(nowP,data.numDir,1)+data.dir.*step;
     data.endP=endP;
     ob= checkobstacle(nowP,nowP+vi*data.dmax,data,data.dmax,data.BD1);
     nextP(count*2-1:count*2,:)=[];
     % flagrec=[flagrec;ob,v];
      flagrec=[ob,v];
    else
       if checkpipe(nowP,nowP+vi*data.dmax, data, data.dmax,data.BD1)>8
 step=2*data.BD; 
 nextP=repmat(nowP,data.numDir,1)+data.dir.*step;
       else
       step=min(max(abs(detaP)-data.BD1/2)/2,data.dmax);
       nextP=repmat(nowP,data.numDir,1)+data.dir*step;
       end
        % nextP=nowP+vi*data.dmax;
        % raddirection=rand();
        % % if raddirection<0.3
        % % nextP=repmat(nowP,data.numDir,1)+data.dir.*data.dmax;
        % % nextP(count*2,:)=nowP+vi*data.dmax;
        % % end
    end
elseif flagpipe==1&&checkpipe(nowP,nowP+data.dir(v,:)*data.BD1, data, data.dmax,data.BD1)>9
 step=2*data.dmax; 
 nextP=repmat(nowP,data.numDir,1)+rand(1)*data.dir.*step;
else


%%步长确认
  [u,u2]=potential_energy(path(end,:),nowP, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
 if u2 <0
       [~,u3]=potential_energy(nowP+data.dir(v,:)*0.01,nowP+data.dir(v,:)*0.01, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
       if u2==u3
       step=rand(1)*data.dmax;

       nextP=repmat(nowP,data.numDir,1)+data.dir*step;

       else
       step=data.dmax-data.dsafe-0.01*u2/(u3-u2); 
  if step<0
     step=abs(step);
  elseif step>10*data.dmax
      step=data.dmax;
  elseif step==0
      step=1;
 end
       nextP=repmat(nowP,data.numDir,1)+data.dir*step;

       end
 elseif u2==0
       step=data.dmax-data.dsafe;  
 
       nextP=repmat(nowP,data.numDir,1)+data.dir*step;

 else
       if u>=0
       step=rand(1)*minobstacle(nowP,nowP,data,data.dmax,data.BD1);
       step(step>=data.dmax) = data.dmax;
       nextP=repmat(nowP,data.numDir,1)+data.dir.*step;

       else
       [~,u3]=potential_energy(nowP,nowP+data.dir(v,:)*0.01, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);           
       if u2==u3
       step=rand(1)*minobstacle(nowP,nowP,data,data.dmax,data.BD1);
       nextP=repmat(nowP,data.numDir,1)+data.dir.*step;

       else
       step=data.dmax-data.dsafe-0.01*u2/(u3-u2); 
  if step<0
     step=0;
     continue;
  elseif step>10*data.dmax
      step=data.dmax;
  elseif step==0
      step=1;
 end
       nextP=repmat(nowP,data.numDir,1)+data.dir*step;
       end
 
       end
 end
end

end


 ind=[]; %% 初始化无效步骤的索引
 %%排除方向
% searchedVoxels=[searchedVoxels;nextP];
 for noDir=1:size(nextP,1) %% 遍历所有可能的移动或者说搜索方向
        %%判断按此方向搜索得到的下一个点是否越界，如果越界则添加其方向索引到ind数组中
            % 判断按此方向搜索得到的下一个点是否已被访问过（也就是在地图上是否被标记过），如果是则添加其方向索引到ind数组中
            % 其中，flag是一个三维数组，用来记录哪些点是已经走过的或者是障碍物。如果某一个点走过了或者是障碍物，就会在这个位置标记1，否则就是0！改为其他不可穿障碍134与管路径点
            [u,~]=  potential_energy(nowP,nextP(noDir,:), data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
               if u==Inf
                  if contains_array(endP, [nextP(noDir,1),nextP(noDir,2),nextP(noDir,3)],nowP)
                      nextP=[nextP(noDir,1),nextP(noDir,2),nextP(noDir,3)];
                      path(end,:)= endP;
                  break;%%结束for循环
                  end
                 ind=[ind;noDir];
               end
            [u,~]=  potential_energy(nextP(noDir,:),nextP(noDir,:), data.mapsize1,data.mapsize2, data, data.dmax,data.BD1,flagrec);
               if u==Inf
                  if contains_array(endP, [nextP(noDir,1),nextP(noDir,2),nextP(noDir,3)],nowP)
                      path(end,:)= endP;
                  break;%%结束for循环
                  end
               
                 ind=[ind;noDir];
                  
               end        
 end
   

   
  
  nextP(ind,:)=[]; %% 删除无效步骤 6*3矩阵-n*3矩阵
    %%处理路径搜索过程中遇到死路（无法继续前进）的情况：转向或后撤
    if isempty(nextP) %% 判断下一步的位置（nextP）是否为空，也就是是否有可走的路径。如果为空，说明当前位置无法再向前，即遇到了死路
%         disp(nowP);
    if contains_array(endP, nowP,path(end,:))
       if(length(dirpath)==5)
           step=0;
       end

    v=find(ismember(data.dir, fix((endP-nowP)/step), 'rows'));
    step=h(endP,nowP);
    if v==dirpath(end-1)
    dirpath(end)=dirpath(end)+step;
    elseif fix((v+1)/2)==fix((dirpath(end-1)+1)/2)
    dirpath(end)=dirpath(end)-step;    
    elseif v~=0
    dirpath=[dirpath,v,step];
    end
        path(end,:)= endP;
        break;%%结束while循环
    end  
    if size(path,1)==1
        nowP=sP;
        dirpath=[nowP,data.dirSE(1),0];
        flagrec=[];
        step=10;
        continue;
    else
        nowP=path(end-1,:); %% 从路径中取出倒数第二个位置，也就是从当前位置的前一个位置，将其设置为当前位置，用于回溯
    end  
        
        if dirpath(end)==0&&length(dirpath(end))~=5
           dirpath(end-1:end)=[];
        end
        v=find(ismember(data.dir, fix((path(end,:)-path(end-1,:))/abs(sum(path(end-1,:)-path(end,:)))), 'rows'));
       if v==dirpath(end-1)
        dirpath(end)=dirpath(end)-abs(sum(path(end-1,:)-path(end,:)));
       elseif fix((v+1)/2)==fix((dirpath(end-1)+1)/2)
        dirpath(end)=dirpath(end)+abs(sum(path(end-1,:)-path(end,:)));
       else
           dirpath=[dirpath,v+mod(v,2)*2-1,abs(sum(path(end-1,:)-path(end,:)))];
       end
        path(end,:)=[];
        if size(path,1)>1
        step=abs(sum(path(end-1,:)-path(end,:)));
        else
            step=1;
        end
        %% 从路径中删除当前的位置，也就是路径的最后一个位置
        continue; %% 跳过本次循环的后续操作，直接开始下一轮循环。也就是在回溯到上一位置后，跳过后续操作，再次尝试寻找下一步的路径。
    end
    %% 初始化 优先级得分列表 距离得分列表 方向得分列表
    pri=[];%% 初始化一个空的优先级得分列表 n*1矩阵-xx
    D=[]; %% 初始化一个空的距离得分列表 n*1矩阵
    Dr=[]; %% 初始化一个空的方向得分列表 n*3矩阵
     Dh=[];
      Dhh=0;
      energy=[];
       energypipe=[];
        w1=0.2;
        flagpipe=0;
         w2=0.87;
         w3=0.1;
         w4=0.1;

         w0=0.2;
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
        [energy(i),~]=potential_energy(nowP,nextP(i,:), data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
       [energypipe(i),~]=checkpipe(nowP,nextP(i,:), data, data.dmax,data.BD1);
       if energypipe(i)>9||any(ismember(pathcir,nextP(i,:), 'rows'))
           flagpipe=1;
       end
       if energypipe(i)<0
           energypipe(i)=-3;
       end
      
    % for mm=1:3:size(markedGroups,1)
    %     if isequal(markedGroups(mm,4), data.fourthColumnValues)&&((~isequal(data.S0, markedGroups(mm,1:3))&&~isequal(data.E0, markedGroups(mm+1,1:3)))||(~isequal(data.S0, markedGroups(mm+1,1:3))&&~isequal(data.E0, markedGroups(mm,1:3))))
    % 
    % 
    % 
    %         Dhh=Dhh+h(nextP(i,:),markedGroups(mm,1:3));
    % 
    %        Dhh=Dhh+h(nextP(i,:),markedGroups(mm+1,1:3));
    % 
    % 
    %     end
    %      Dh(i)=Dhh;
    % 
    % end 
   Dhh=0;
        for mm=1:3:size(data.branch,1)
         if~(isequal(data.S, data.branch(mm,1:3))&&isequal(data.E, data.branch(mm+1,1:3)))||(isequal(data.S,data.branch(mm+1,1:3))&&isequal(data.E,data.branch(mm,1:3)))



            Dhh=Dhh+h(nextP(i,:),data.branch(mm,1:3));

           Dhh=Dhh+h(nextP(i,:),data.branch(mm+1,1:3));


        end
         Dh(i)=Dhh;

    end 
if any(ismember(pathcir,nextP(i,:), 'rows'))
    energypipe(i)=100;  
    flagpipe=1;
end

    end
    %% 归一化处理
    %%这段代码的作用是对三个变量（pri、D、Dr）进行归一化，然后把他们相加，并找出加和后的结果中最小值对应的索引
    temp1=myMapminmax(Dh);
    temp2=myMapminmax(D);
    temp3=myMapminmax(Dr);
    temp4=myMapminmax(energy);
    temp5=myMapminmax(energypipe);
    if flagpipe==1
        index=find(energypipe>=10);
        if ~isempty(index)
      temp5(index) =energypipe(index); 
        end
        w3=0.5;  w4=0.1;
        % w0=data.w0;
    if ~isempty(data.dir1)
    w3=0.5;w2=1;
    % w0=data.w0;
    end
    end
    if junpflag==0
    temp=w0*temp1+w2*temp2+w3*temp3+w4*temp4+w1*temp5;
    else
      temp=temp5;
    end%%计算以上三个归一化操作后的结果的和！？改系数？
    %%寻找temp中值等于其最小值的所有元素的索引，这个索引被存放在no变量中。
    [no]=find(temp==min(temp));%% find函数会返回所有满足条件的值的索引，返回值no是一个数组
    if isempty(no)
    no=[1];
    end
 %% 只取索引数组中的第一个元素，如果有多个元素满足条件，我们随机考虑。
no=no(randi(length(no)));
 nextP=nextP(no,:);%% 在候选点集合nextP中选取索引为no的行作为下一步位置。
    if isempty(data.dir1)
    dir=(nextP-nowP)/h(nextP,nowP); %%计算参考方向，即当前位置到下一个位置的向量。
    end


    step=h(nextP,nowP);
    
    v=find(ismember(data.dir, fix((nextP-nowP)/step), 'rows'));
    nowP=nextP; %%更新当前位置为下一步的位置
    %%在flag数组中，将新的当前位置（nowP）对应的元素值设置为1。数组flag用于记录已经走过的位置。
    path=[path;nowP];% 将新的当前位置添加到路径path的末尾。路径path用于记录移动过程中经过的所有位置。
   pathcir=[nowP;pathcir];
   if size(pathcir,1)>=10000
   pathcir(1,:)=[];
   end
      if size(path,1)>=100
   path(1,:)=[];
   end
    %% 检查是否已经找到结束位置
    %%ismember 函数用于判断 'nowP' 是否在 'data.E' 的行中存在。
    %%如果当前位置 'nowP' 已经到达结束位置（即 'nowP' 在 'data.E' 中），ismember 函数返回真 (1)。
    %%这里 "rows" 是 ismember 函数的参数，指定操作的维度为行。
  %%方向路径记录
if v==dirpath(end-1)
dirpath(end)=dirpath(end)+step;
elseif fix((v+1)/2)==fix((dirpath(end-1)+1)/2)
dirpath(end)=dirpath(end)-step;    
elseif v~=0
dirpath=[dirpath,v,step];
end

if dirpath(end)<0
    dirpath(end)=-dirpath(end);
    dirpath(end-1)=dirpath(end-1)+mod(dirpath(end-1),2)*2-1;
end

if len==length(dirpath)
    if isempty(circledir1)||dirpath(end)>circledir1(end)
    circledir1=dirpath(end-1:end);
  
    else
        circlenumber=circlenumber+1;
    
      if circlenumber>2
       ob= computePipeAABB(nowP,nowP-data.dir(dirpath(end-1),:)*dirpath(end),data.BD1);
       flagrec=[flagrec;ob,v];
       circlenumber=0;
       circledir1=[];
      end
    end
else
    circledir1=[];
end
len=length(dirpath);




    if contains_array(endP, nowP,path(end,:))
        path(end,:)= endP;
        break;%%结束while循环
    end  


end
end

%% 计算路径长度
%%路径的优度是以路径的总长度（这里是所有连续两点间的欧氏距离之和）
xdir=dirpath(4:end);
xdir1=[dirpath(4),0];
fit=0;
u=0;
curve_count = 0; %% 初始化弯曲次数
penetration_count = 0; %% 初始化穿墙次数
energy_total= 0; %%初始化能量总和
energy_total1= 0;
curve(1,:)=sP; 
ii=2;
 path1=[];
 nowP=dirpath(1:3);
 curve_count = (length(dirpath)-5)/2;
tempstep=12;
curve_len=0;
 palh1=0;


for i=4:2:length(dirpath)-1
    %% 计算路径中每两个连续点的欧氏距离，并将它加到总路径长度之上。
    %%在这里，norm(path(i,:)-path(i+1,:))表示的是点path(i,:)和点path(i+1,:)的欧氏距离。
    temp=dirpath(i:i+1);
    fit=fit+temp(2);
    % if i~=1&&temp(2)<2*data.BD1&&i~=length(dirpath)-1
    %     curve_len=curve_len+10000000;
    % end
    tempP1=curve(end,:);
    tempP2=tempP1+data.dir(temp(1),:)*tempstep;
    ij=tempstep;
    %% 累计节点能量值固定步长
    while ij<temp(2)
        if temp(2)<tempstep
            break;
        end
    
    [u,~]=potential_energy(tempP1,tempP2, data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
     energy_total = energy_total + u;
     energy_total1= energy_total1 + u;
    [u,~,np]= checkpipe(tempP1,tempP2, data, data.dmax,data.BD1);
            if u>8
    u=u*1000000000;
            else
            endnode=tempP2;
            xdir1=dirpath(4:i+1);
            xdir1(end)=ij;
            
            end
                        if u<0
palh1=palh1+tempstep;
            end
    energy_total = energy_total + u;
    tempP1=tempP2;
    tempP2=tempP1+data.dir(temp(1),:)*tempstep;
    ij=ij+tempstep;
    end
    %% 判断是否发生弯曲并累计弯曲次数
    curve(ii,:)=curve(end,:)+data.dir(temp(1),:)*temp(2);
        aabb = computePipeAABB(curve(ii-1,:),curve(ii,:), data.BD);
        path1=[path1;aabb];
        ii=ii+1;  
        if u<8 && temp(2)~=0
        [u,~,np]= checkpipe(tempP1,curve(end,:), data, data.dmax,data.BD1);
        if u>8
           u=u*1000000000;
        else
            endnode=curve(end,:);
            xdir1=dirpath(4:i+1);
        
        end
                    if u<0
palh1=palh1+tempstep;
            end
        else
        [u,~,np]= checkpipe(tempP1,curve(end,:), data, data.dmax,data.BD1);
        if u<8
            endnode=curve(end,:);
            xdir1=dirpath(4:i+1);
        end  
                    if u<0
palh1=palh1+tempstep;
            end
        end
        energy_total = energy_total + u;
        [u,~]=potential_energy(tempP1,curve(end,:), data.mapsize1,data.mapsize2, data, data.dmax,data.BD1);
         if u==inf
            energy_total=inf;
         end 
    energy_total = energy_total + u;
    energy_total1= energy_total1 + u;
        

    %% 判断是否穿墙并累计穿墙次数

        penetration_count =0;
end


%根据要求计算各部分
dis_len = fit+ pi / 2 * data.BD * curve_count -data.BD * curve_count;
curve_len =(3^curve_count)*data.curve+curve_len;
penetration_len = penetration_count * data.unit_xyz;
palh=0;
 for mm=1:3:size(data.branch,1)
   if~(isequal(data.S, data.branch(mm,1:3))&&isequal(data.E, data.branch(mm+1,1:3)))||(isequal(data.S,data.branch(mm+1,1:3))&&isequal(data.E,data.branch(mm,1:3)))

p1=closest_point_on_path(data.branch(mm,1:3),{dirpath},data);
p2=closest_point_on_path(data.branch(mm+1,1:3),{dirpath},data);
palh=palh-h(p1,data.branch(mm,1:3))-h(data.branch(mm+1,1:3),p2);
palh0=palh-h(p1,data.branch(mm,1:3))-h(data.branch(mm+1,1:3),p2);
   end
 end 
 if size(data.branch,1)>3&&(isequal(data.S, sP)&&isequal(data.E,endP))
 palh=palh/(size(data.branch,1)-3)+palh1/2*data.plw;
 palh0=palh/(size(data.branch,1)-3)*3+palh1/2;
 end
% 计算总能量值
energy_total = energy_total;
% 计算综合函数，注意权重之和为1，这里暂定各项权重为0.25
fit = 0.6 * dis_len + 0.4 * curve_len + 0.5 * (penetration_len-palh) + 0.4*energy_total;
if curve(end,:)~=endP
aabb = computePipeAABB(curve(end,:),endP, data.BD);
path1=[path1;aabb];
end
aabb = computePipeAABB(curve(end,:),endP-data.dir(dirpath(end-1),:)*data.BD*1.2, data.BD);
path2=[path1;aabb];
% 下面的if语句是创建一个包含路径长度和路径信息的结果对象，如果函数有第二个输出参数的话
%%nargout 是一个内置变量，表示函数的输出参数的数量。当函数需要返回的参数数量大于 1 时，执行以下操作。
if nargout>1
    result.fit=fit;
    result.length=dis_len;
    result.path=curve;
    result.path1=path1;
    result.path2=path2;
    result.endnode=endnode;
    result.dirpath=[sP,xdir];
    result.xdir=xdir;
    result.xdir1=xdir1;
    result.dis_len = dis_len;
    result.curve_len = curve_len;
    result.curve = curve;
    result.penetration_len = penetration_len;
    result.energy_total = energy_total1;
    result.curve_count=curve_count;
     result.palh1=palh1;
     % result.searchedVoxels =searchedVoxels;
end
end
