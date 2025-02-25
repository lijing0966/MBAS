function [xbest,fbest,distHistory] = BAS_1(dis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%初始化部分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=length(dis);     % k表示距离矩阵dis的行列数,即，访问城市个数为k，天牛须维数为k，访问路径从第1个城市到第k个城市，并返回第1个城市
m=5;               % m表示天牛群体的个体数量 m=5
x_store=rand(k,m); % 天牛种群初始化，随机生成0~1之间的数
f_store=zeros(1,m);% 各天牛所对应的适应度值（使得目标函数最小的最优解）
X_store=zeros(k,m);% 各天牛所对应的最优路径（访问的城市序列）
for i=1:m
    x=x_store(:,i);    % 依次选取单个天牛，记为x
    [~,xI]=sort(x);    % 将坐标向量x中每个元素进行从小到大排列
    [~,xII]=sort(xI);  % 得到城市序列xII，用来计算当前天牛的适应度值（看看目标函数是不是最小的）
    v=zeros(k,k);      % 构造一个k*k的二进制数矩阵，取值0或1，用来表示是当前路径否存在第i个城市到第j个城市的路径
    for p=1:k-1        % 按照城市序列xII中的情况赋值
        v(xII(p),xII(p+1))=1;
    end
    v(xII(k),xII(1))=1;
    xbest=xII;           % 初始城市序列xII设为当前最优解
    fbest=f(v,dis);      % 函数f()是待优化的目标函数
    X_store(:,i)=xII;    % 初始化第i个天牛所对应的最优路径
    f_store(i)=fbest;    % 初始化第i个天牛所对应的最优解
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%迭代部分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta=0.95;
c=5;%ratio between step and d0
step=1;%initial step set as the largest input range $\delta$
n=500;%iterations  n=500
distHistory = zeros(1,n);          %历史值，每一次迭代都产生一个最优适应度值作为历史值存起来。

for r=1:n
    d0=step/c;            % xleft表示左须坐标，xright表示右须坐标，x表示质心坐标，用d0表示两须之间距离
    for i=1:m
        x=x_store(:,i);
        fbest=f_store(i);
        dir=rand(k,1);    % intial value，随机生成一个k行1列的矩阵，且其中的各元素范围在0到1之间（产生一个随机向量，表示须的方向）
        dir=dir/norm(dir);% norm(dir)表示2范数，就是矩阵的模sqrt(sum(abs(dir.^2)))，每个元素平方求和再开方。
                          % 进行归一化处理，让矩阵dir中每个元素值除以它的模（变得更小)
        xleft=x+dir*d0/2; % 左须坐标
        xright=x-dir*d0/2;% 右须坐标
        % 计算左须在旅行商问题上的适应度值                  
        [~,xI]=sort(xleft);
        [~,xII]=sort(xI);
        v=zeros(k,k);
        for q=1:k-1        
            v(xII(q),xII(q+1))=1;
        end
        v(xII(k),xII(1))=1;
        fleft=f(v,dis);   % 左须的气味强度
        % 计算右须在旅行商问题上的适应度值  
        [~,xI]=sort(xright);
        [~,xII]=sort(xI);
        v=zeros(k,k);
        for q=1:k-1        
            v(xII(q),xII(q+1))=1;
        end
        v(xII(k),xII(1))=1;
        fright=f(v,dis); % 右须的气味强度
        % 确定下一步的位置
        x=x-step*dir.*sign(fleft-fright); % fleft>fright，则sign()=1，向着右须方向行进距离step；fleft<fright，则sign()=-1，向着左须方向行进距离step
        x_store(:,i)=x;
        [~,xI]=sort(x);    
        [~,xII]=sort(xI); 
        v=zeros(k,k);     
        for q=1:k-1        
            v(xII(q),xII(q+1))=1;
        end
        v(xII(k),xII(1))=1;
        fnew=f(v,dis);          % 根据计算的x更新要优化的函数f
        % 如果本次迭代中，f值更小，则将xII作为xbest，f作为fbest
        if fnew<fbest 
            xbest=xII;
            fbest=fnew;
            X_store(:,i)=xII; 
            f_store(i)=fbest;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %混合遗传算法开始——随机设置两个断点，用于进行反转、交换，让x的排序更优
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        InsertionPoints=sort(ceil(k*rand(1,2)));
        I=InsertionPoints(1);
        J=InsertionPoints(2);
        setnum=3;
        Route=zeros(k,setnum);
        
        for p = 1:setnum
            Route(:,p)=x;
            sub=x(I:J);
            switch p
                case 2 
                    Route(I:J,p)=sub(end:-1:1);             % 反转突变（Inversion mutation）
                case 3  
                    Route([I J],p)=sub([end 1]);            % 交换突变（Exchange mutation）            
                otherwise
            end
        end

        for p=1:setnum
            [~,xI]=sort(Route(:,p));
            [~,xII]=sort(xI);
            v=zeros(k,k);
            for q=1:k-1        
                v(xII(q),xII(q+1))=1;
            end
            v(xII(k),xII(1))=1;
            fnew = f(v,dis);
            if fnew < fbest
                x_store(:,i)=Route(:,p);
                xbest=xII;
                fbest=fnew;
                X_store(:,i)=xII; 
                f_store(i)=fbest;                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %混合遗传算法结束
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %从m个天牛中选择最优值
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:m
        if f_store(i)<fbest
            xbest=X_store(:,i);
            fbest=f_store(i);
        end
    end
    %%%%%%%%%%%
    distHistory(r) = fbest;
    step=step*eta;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%迭代部分结束
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 算法结束
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


