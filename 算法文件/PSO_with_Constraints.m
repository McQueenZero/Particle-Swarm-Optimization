%% 清空环境
clc;clear;close all;
%% 作目标函数图像
 figure('Name','Schaffer F6 函数绘图');
 DrawOptimOBJ();
%% 目标函数与约束条件
%    objval = OptimOBJ(x) 
%    nvars = 2  （变量个数）
%    约束条件为：
%    |x(1)| <= 100;            
%    |x(2)| <= 100;                 
c1 = [1,0];
c2 = [0,1];
fun = @OptimOBJ;
cons1 = @(X)(abs(c1*X)<=100);
cons2 = @(X)(abs(c2*X)<=100);
%% 设置种群参数
dim = 2;                            % 空间维数=nvars
sizepop = 100;                      % 初始种群个数
ger = 200;                          % 最大迭代次数    
xlimit_max = 200*ones(dim,1);       % 设置位置参数限制
xlimit_min = -200*ones(dim,1);
vlimit_max = 2*ones(dim,1);         % 设置速度限制
vlimit_min = -2*ones(dim,1);
c_1 = 0.5;                          % 原始惯性权重
c_2 = 2;                            % 自我学习因子
c_3 = 2;                            % 群体学习因子 
PV = 10^10;                         % 原始惩罚值（Punishment Value）
%% 种群初始化
%  首先随机生成初始种群位置
%  然后随机生成初始种群速度
%  然后初始化个体历史最佳位置，以及个体历史最佳适应度
%  然后初始化群体历史最佳位置，以及群体历史最佳适应度
for i=1:dim
    for j=1:sizepop
        pop_x(i,j) = xlimit_min(i)+(xlimit_max(i) - xlimit_min(i))*rand;  % 初始种群的位置
        pop_v(i,j) = vlimit_min(i)+(vlimit_max(i) - vlimit_min(i))*rand;  % 初始种群的速度
    end
end                 
pbest = pop_x;                                % 每个个体的历史最佳位置
for j=1:sizepop
    if cons1(pop_x(:,j))
        if cons2(pop_x(:,j))
            fitness_pbest(j) = fun(pop_x(:,j));     % 每个个体的历史最佳适应度
        else
            fitness_pbest(j) = PV/10^8;
        end
    else
        fitness_pbest(j) = PV/10^8;
    end
end
% 初始化种群时实际惩罚值较小，达到扩大搜索空间的目的
gbest = pop_x(:,1);                           % 种群的历史最佳位置
fitness_gbest = fitness_pbest(1);             % 种群的历史最佳适应度
for j=1:sizepop
    if fitness_pbest(j) < fitness_gbest       % 如果求最小值，则为<; 如果求最大值，则为>;
        gbest = pop_x(:,j);                   % 记录种群的历史最佳位置
        fitness_gbest=fitness_pbest(j);
    end
end
 
 
%% 粒子群迭代更新
%    更新速度并对速度进行边界处理    
%    更新位置并对位置进行边界处理
%    进行约束条件判断并计算新种群各个个体的适应度
%    新适应度与个体历史最佳适应度做比较
%    个体历史最佳适应度与种群历史最佳适应度做比较
%    再次循环或结束
 
iter = 1;                        %迭代次数
record = zeros(ger, 1);          % 记录器
while iter <= ger
    for j=1:sizepop
        %    更新速度并对速度进行边界处理
        %    原始惯性权重因子叠加线性函数，使实际惯性权重随迭代次数增加而减小
        pop_v(:,j)= c_1*(-1/ger*iter+1) * pop_v(:,j) + c_2*rand*(pbest(:,j)-pop_x(:,j))+c_3*rand*(gbest-pop_x(:,j));% 速度更新
        for i=1:dim
            if  pop_v(i,j) > vlimit_max(i)
                pop_v(i,j) = vlimit_max(i);
            end
            if  pop_v(i,j) < vlimit_min(i)
                pop_v(i,j) = vlimit_min(i);
            end
        end
        
        %    更新位置并对位置进行边界处理
        pop_x(:,j) = pop_x(:,j) + pop_v(:,j);% 位置更新
        for i=1:dim           
            if  pop_x(i,j) > xlimit_max(i)
                pop_x(i,j) = xlimit_max(i);
            end
            if  pop_x(i,j) < xlimit_min(i)
                pop_x(i,j) = xlimit_min(i);
            end
        end        
  
        %    进行约束条件判断并计算新种群各个个体的适应度
        if cons1(pop_x(:,j))
            if cons2(pop_x(:,j))
                fitness_pop(j) = fun(pop_x(:,j));   % 当前个体的适应度
            else
                fitness_pop(j) = PV/(1+1.05^(-iter+ger/3));
            end
        else
            fitness_pop(j) = PV/(1+1.05^(-iter+ger/3));
        end
        % 进行适应值比例变换：
        % 迭代时原始惩罚值叠加类Sigmoid函数，相当于实际惩罚值随着代际增大而增大。
        % 迭代开始时，实际惩罚值小，抑制竞争，扩大搜索范围；
        % 迭代后期时，实际惩罚值大，鼓励竞争，加快收敛速度。
        
        %    新适应度与个体历史最佳适应度做比较
        if fitness_pop(j) < fitness_pbest(j)       % 如果求最小值，则为<; 如果求最大值，则为>; 
            pbest(:,j) = pop_x(:,j);               % 更新个体历史最佳位置            
            fitness_pbest(j) = fitness_pop(j);     % 更新个体历史最佳适应度
        end   
        
        %    个体历史最佳适应度与种群历史最佳适应度做比较
        if fitness_pbest(j) < fitness_gbest        % 如果求最小值，则为<; 如果求最大值，则为>; 
            gbest = pbest(:,j);                    % 更新群体历史最佳位置  
            fitness_gbest=fitness_pbest(j);        % 更新群体历史最佳适应度  
        end    
    end
    
    record(iter) = fitness_gbest;                  % 最小值记录
    
    iter = iter+1;
 
end
%% 迭代结果输出
figure('Name','粒子群优化算法迭代过程');
plot(record);title('收敛过程'); hold on;
disp(['最优值：',num2str(fitness_gbest)]);
disp('变量取值：');
fprintf('%.6f\t',gbest);

