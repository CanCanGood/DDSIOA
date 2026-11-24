%%定义遗传算法函数 GA，它接收四个输入参数：x (当前解)，y (当前解的函数值)，option (算法的参数)，data (地图数 recording.length(1+iter)=length1;据)，
%%并返回三个输出参数：bestY (最优解的函数值)，bestX (最优解)，recording (记录每次迭代的最优解和平均值)
function [bestY,bestX,recording]=GACELL(x,y,option,data)

    %% 遗传算法
    %% 初始化
    data.curve=5000;
    %%初始化 recording 结构体用来记录每次迭代的最优解和平均优化值
    recording.bestFit=zeros(option.maxIteration+1,1);
    recording.meanFit=zeros(option.maxIteration+1,1);
    recording.wFit=zeros(option.maxIteration+1,1);
    %% 更新记录
    %%找出初始解中的最优解，并记录下最优解的位置和值
    [y_g,position]=min(y);
    x_g=x{position};
     data.dir1=[];
     recording.searchedVoxels=[];
    %%记录每个粒子的位置和最优解
    y_p=y;
    x_p=x;
    iter=0;
         for i=1:option.numAgent
           [y_p(i),newX{i}]=option.fobj(option,data,x_p{i});
        end
         [y_g,newX{i},result]=option.fobj(option,data,x_g);
                    curvecount=result.curve_count;
                    energytotal=result.energy_total;
                    length1=result.length;
                    palh1=result.palh1;

        recording.count(1+iter)=curvecount;
        recording.energytotal(1+iter)=energytotal;
        recording.length(1+iter)=length1;
        recording.palh1(1+iter)=palh1;
        
        recording.searchedVoxels = [recording.searchedVoxels; result.searchedVoxels];
        recording.searchedVoxels = unique(recording.searchedVoxels, 'rows');
        newX=x;
    newY=y;%%记录全局的最优解和平均优化值
    recording.bestFit=y_g;
    recording.meanFit=mean(y_p);
    bestY=y_g;
    bestX=x_g;
    recording.wFit=max(y_p);
    disp(['GA,iter:0,minFit:',num2str(y_g)]);
    %% 开始更新
    %%开始循环迭代
    for iter=1:option.maxIteration
        data.rad=1-randi(option.maxIteration-iter+1)/option.maxIteration;
        disp(['GA,iter:',num2str(iter),',minFit:',num2str(y_g)]);
        %% 选择：竞标赛法选择个体
        for i=1:option.numAgent*2
            maxContestants=ceil(randi(option.numAgent)*option.p1_GA);
            index=randperm(option.numAgent,maxContestants);
            [~,position]=min(y_p(index));
            parent(i)=index(position(1));
        end

        %% 交叉变异
        for i=1:option.numAgent
            xx=x_p{parent(i)};
            xx1=x_p{parent(i+option.numAgent)};
            rand1=rand;
            if rand1<0.4 %多点交叉

                 lengthpath=length(xx);
                 lengthpath1=length(xx1);
                ii=randi(lengthpath/2)*2;  
                iiy=randi(lengthpath1/2)*2;  

            xx(ii)=rand()*(xx(ii));
           xx1(iiy)=rand()*xx1(iiy);
            sP=data.S;
            endP=data.S;
             for ix=0:2:ii-1
            sP=sP+data.dir(xx(ix+1),:)*xx(ix+2);
             end

             for ix=0:2:iiy-1
            endP=endP+data.dir(xx1(ix+1),:)*xx1(ix+2);
             end

             if h(sP,data.E)< h(endP,data.E)
             [~,xxx]=option.fobj(option,data,xx(1:ii),xx1(1:iiy),xx(1:ii));
              xx(ii)=x_p{parent(i)}(ii)-xx(ii);
              xx=[xx1(1:iiy),xxx,xx(ii-1:end)];
              
             else

            [~,xxx]=option.fobj(option,data,xx(1:ii),xx(1:ii),xx1(1:iiy));
              xx1(iiy)=x_p{parent(i+option.numAgent)}(iiy)-xx1(iiy);
              xx=[xx(1:ii),xxx,xx1(iiy-1:end)];
             end



            [newY(i),newX{i}]=option.fobj(option,data,xx);
            else%多点变异
             yy=xx;
             ii=0;
             iiy=0;
            while 1
                 lengthpath=length(xx);
               if iiy>=lengthpath-4
                   break;
               end


                ii1=randi(lengthpath/2-2-iiy/2)*2+iiy;

            while xx(ii1) ==0
                  if ii1==length(xx)-4
                    break;
                  end

               ii1=randi(lengthpath/2-2-iiy/2)*2+iiy;
           end
           if xx(ii1) ==0

                    break; 
           end

            ii=ii1;


           iiy=ii+2+randi((lengthpath-ii)/2-1)*2;
           while yy(iiy) ==0
                iiy= ii+2+randi((lengthpath-ii)/2-1)*2;
                 if iiy==length(xx)
                    break;
                end
           end
           if yy(iiy) ==0

                    break;

           end
           xx(ii)=rand()*xx(ii);
            yy(iiy)=rand()*yy(iiy);
            [~,xxx]=option.fobj(option,data,xx(1:ii),xx(1:ii),yy(1:iiy));
            yy(iiy)=xx(iiy)-yy(iiy);

           if iiy==length(xx)
            xx=[xx(1:ii),xxx,yy(iiy-1:end)];
            yy=[xx(1:ii),xxx,yy(iiy-1:end)];
            break;
           else
                    xx=[xx(1:ii),xxx,yy(iiy-1:end)];
            yy=[xx(1:ii),xxx,yy(iiy-1:end)];
            end

            end
            [newY(i),newX{i}]=option.fobj(option,data,xx);
            end


        end
        %% 重新计算适应度值 
        for i=1:option.numAgent
                if newY(i)<y_p(i)
                y_p(i)=newY(i);
                x_p{i}=newX{i};
                if y_p(i)<y_g
                    y_g=y_p(i);
                    x_g=x_p{i};
                    [newY(i),newX{i},result]=option.fobj(option,data,x_g);
                    curvecount=result.curve_count;
                    energytotal=result.energy_total;
                    length1=result.length;
                    palh1=result.palh1;
                end
            end
        end
        %% 更新记录
        recording.bestFit(1+iter)=y_g;
         recording.bestFit1(1+iter)=y_g;
        recording.meanFit(1+iter)=mean(y_p);
        recording.wFit(1+iter)=max(y_p);
        recording.count(1+iter)=curvecount;
        recording.energytotal(1+iter)=energytotal;
        recording.length(1+iter)=length1;
        recording.palh1(1+iter)=palh1;
        recording.searchedVoxels(1+iter)=result.searchedVoxels;
    end
    
          bestY=y_g;
          bestX=x_g;
    recording.bestFit(1+iter)=y_g;
    recording.meanFit(1+iter)=mean(y_p);
    recording.x_p=x_p;
    recording.y_p=y_p;
     recording.bestFit1(1+iter)=y_g;
    recording.x_g=[data.S,x_g,data.E];
            recording.count(1+iter)=curvecount;
        recording.energytotal(1+iter)=energytotal;
        recording.length(1+iter)=length1;
        recording.palh1(1+iter)=palh1;
[recording.wFit(1+iter),position] =max(y_p);
recording.x_W=x_p{position};
end

    
 