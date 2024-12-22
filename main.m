clear all;
close all;
tic
global seru product_time product demand setup Y block r_block  parents pop_size generation p_c p_mu Duedate fitness solution_set
pop_size=100;%之后可以取300 
generation=200;%之后可以取500
seru=3;
p_c=0.75;%交叉
p_mu=0.05;%变异
product_time=[8 9.5 11.5 9 10.5 16
    19 12 15 12 13 9
    16.5 16 12.5 15 13.5 10];%./10 6个产品在3个seru的加工时间
product=size(product_time,2); %6个product
demand=[80 90 70 85 88 83];
setup=[4 9 8 6 5 6]; %每种产品的设置时间
block=[2 1.5 1 3 1.5 4 
       5 2 1 2 3 2.5 
       3.5 1 2.5 1.5 3.5 4];%6个产品在3个seru中的阻塞时间
r_block=[0.1 0.2 0.3 0.15 0.25 0.2
         0.34 0.23 0.13 0.1 0.4 0.3
         0.2 0.1 0.4 0.3 0.21 0.13];%3个seru生产6个产品发生阻塞的概率
Duedate=[1920 2360 2400 1840 2380 2400];

% feasible_matrix=zeros(seru,product); 
% feasible_matrix(product_time~=0)=1;（i，j）数值为1时，表示该seru i可以装配产品j
parents=zeros(pop_size,seru*product);%用来装载解
fitness=zeros(1,pop_size);%用来装载适应度
Y=zeros(1,generation);
 
%% 生成初始种群,随机生成订单调度初始种群
init=1;
while init<=pop_size%用来生成初始种群
    ratio=rand(seru,product);%在矩阵的前三行生成调度结果（0-1）之间的随机数3*6
    for i=1:seru
        for j=1:product
            if rand<0.05
                ratio(i,j)=0;
            end
        end
    end
    parents(init,:)=reshape(ratio',seru*product,1)';%把6*3转成 18*1再转置 然后init=10 总共10行
    init=init+1;
end
%parents
%% 迭代
for g=1:generation
    allocated_set=zeros(pop_size,seru*product);%10*18
    for che=1:pop_size
        allocated=zeros(seru,product);%3*6
        ratio_row=parents(che,:);%1*18 一行
        ratio=reshape(ratio_row,product,seru)';%转置 3*6
        u=sum(ratio);
        for v=1:product
            if  u(v)==0
                ratio(1,v)=1;
            end
        end
        u=sum(ratio);
        ratio=ratio./repmat(u,seru,1);%对分配比进行标准化 求出对应占比 sum(ratio)竖着求和 得到一行 1*6 除数是3*6 每行是3个seru对应产品的ratio的总和
        allocated_quan=ratio.*repmat(demand,seru,1);%占比*需求 
        allocated_quan=round(allocated_quan);%四舍五入 整数 3*6
        allocated_sum=sum(allocated_quan);%1*6 竖着求和 分配每个产品对应需求
        for w=1:product
            if   allocated_sum(w)==0
                allocated_quan(1,w)=demand(w);
            end
        end
        for j=1:product
            sit=find(allocated_quan(:,j)~=0,1);%find(x,n)返回与X中的非零元素对应的前n个位置 化成01矩阵 找到第一个有分配projuct j的seru
            allocated_quan(sit,j)=demand(j)-allocated_sum(j)+allocated_quan(sit,j);%修正解 确保分配给seru的子批量等于该产品的需求
        end
      
        allocated_set(che,:)=reshape(allocated_quan',seru*product,1)';%10*18 把每一次得到的allocated quan分配比的值都储存在一行中
        allocated(allocated_quan~=0)=1;
        setup_time=repmat(setup,seru,1);%3*6的矩阵
        setup_time=setup_time.*allocated;
        for i=1:seru
            for j=1:product
                if setup_time(i,j)~=0
                    setup_time(i,j)=0;
                    break
                end
            end
        end         
        block_time=block.*r_block;
        processing_time=allocated_quan.*product_time+setup_time+block_time;%所占用的时间
        %交付期
        for t=2:product 
            processing_time(:,t)=processing_time(:,t)+processing_time(:,t-1); 
        end
         processing_time1=max(processing_time);
         makespan=max(processing_time1);%对矩阵求最大值
      
        flag=0;
        for j=1:product
            if processing_time1(j)>Duedate(j)
                flag=flag+1;
            end
        end
        if flag>0
            F=max(makespan)+flag*1000;%惩罚函数
        else
            F=max(makespan);
        end
        
         fitness(che)=F; 
    end
    [opt_fit,index]=min(fitness);%某次迭代过程中的最优是最大的适应度
    opt_solution=allocated_set(index,:);
    Y(g)=opt_fit; %把每次的最优适应度保存
 
    %% 选择：二元锦标赛
    select=zeros(pop_size,seru*product);
    for s=1:pop_size
        k=fix(rand*pop_size)+1;%随机找到一个具体分配结果
        kk=fix(rand*pop_size)+1;
        if fitness(k)<=fitness(kk)
            select(s,:)=parents(k,:);%如果方案K的适应度更高,就选择K
        else
            select(s,:)=parents(kk,:);
        end
    end
    parents=select;
    %% 交叉：随机找到两个父代，交换中间某一段的数值 代数交叉
    crossn=fix(pop_size/2);
    for c=1:crossn %交叉
        r_c=rand;
        if r_c<p_c
            k=fix(rand*pop_size)+1; %fix函数向0方向取整。
            kk=fix(rand*pop_size)+1;
            Lambda=rand;
            child1=Lambda*parents(k,:)+(1-Lambda)*parents(kk,:);
            child2=(1-Lambda)*parents(k,:)+Lambda*parents(kk,:);
            parents(k,:)=child1;
            parents(kk,:)=child2;
        end
    end
    %% 变异：随机找到一个父代，将某两个位置交换位置 非均匀变异
    b=2;
    for m=1:pop_size
        r_mu=rand;
        if r_mu<p_mu
            k=fix(rand*seru*product)+1;
            if rand<0.5
                parents(m,k)=parents(m,k)+(1-parents(m,k))*(1-rand^(1-g/generation)^b);
            else
                parents(m,k)=parents(m,k)-parents(m,k)*(1-rand^(1-g/generation)^b);
            end
        end
        Pa=rand;
        for i=1:seru
            for j=1:product
                if Pa<0.05
                    parents(m,product*(i-1)+j)=0;
                end
            end
        end
    end
end    
[makespan,index]=min(fitness);
allocation=allocated_set(index,:);


X=1:generation;
plot(X,Y);
hold on
xlabel('Iteration Times');
ylabel('minimum makespan');
final_makespan=makespan;
final_scheduling_result=allocation;
disp('最小的makespan为')
disp(final_makespan)
disp('此时订单调度情况为')
disp(final_scheduling_result)
toc
disp(['运行时间: ',num2str(toc)]);