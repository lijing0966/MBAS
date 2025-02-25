function [xbest,fbest,distHistory] = BAS_1(dis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ʼ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=length(dis);     % k��ʾ�������dis��������,�������ʳ��и���Ϊk����ţ��ά��Ϊk������·���ӵ�1�����е���k�����У������ص�1������
m=5;               % m��ʾ��ţȺ��ĸ������� m=5
x_store=rand(k,m); % ��ţ��Ⱥ��ʼ�����������0~1֮�����
f_store=zeros(1,m);% ����ţ����Ӧ����Ӧ��ֵ��ʹ��Ŀ�꺯����С�����Ž⣩
X_store=zeros(k,m);% ����ţ����Ӧ������·�������ʵĳ������У�
for i=1:m
    x=x_store(:,i);    % ����ѡȡ������ţ����Ϊx
    [~,xI]=sort(x);    % ����������x��ÿ��Ԫ�ؽ��д�С��������
    [~,xII]=sort(xI);  % �õ���������xII���������㵱ǰ��ţ����Ӧ��ֵ������Ŀ�꺯���ǲ�����С�ģ�
    v=zeros(k,k);      % ����һ��k*k�Ķ�����������ȡֵ0��1��������ʾ�ǵ�ǰ·������ڵ�i�����е���j�����е�·��
    for p=1:k-1        % ���ճ�������xII�е������ֵ
        v(xII(p),xII(p+1))=1;
    end
    v(xII(k),xII(1))=1;
    xbest=xII;           % ��ʼ��������xII��Ϊ��ǰ���Ž�
    fbest=f(v,dis);      % ����f()�Ǵ��Ż���Ŀ�꺯��
    X_store(:,i)=xII;    % ��ʼ����i����ţ����Ӧ������·��
    f_store(i)=fbest;    % ��ʼ����i����ţ����Ӧ�����Ž�
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta=0.95;
c=5;%ratio between step and d0
step=1;%initial step set as the largest input range $\delta$
n=500;%iterations  n=500
distHistory = zeros(1,n);          %��ʷֵ��ÿһ�ε���������һ��������Ӧ��ֵ��Ϊ��ʷֵ��������

for r=1:n
    d0=step/c;            % xleft��ʾ�������꣬xright��ʾ�������꣬x��ʾ�������꣬��d0��ʾ����֮�����
    for i=1:m
        x=x_store(:,i);
        fbest=f_store(i);
        dir=rand(k,1);    % intial value���������һ��k��1�еľ��������еĸ�Ԫ�ط�Χ��0��1֮�䣨����һ�������������ʾ��ķ���
        dir=dir/norm(dir);% norm(dir)��ʾ2���������Ǿ����ģsqrt(sum(abs(dir.^2)))��ÿ��Ԫ��ƽ������ٿ�����
                          % ���й�һ�������þ���dir��ÿ��Ԫ��ֵ��������ģ����ø�С)
        xleft=x+dir*d0/2; % ��������
        xright=x-dir*d0/2;% ��������
        % ���������������������ϵ���Ӧ��ֵ                  
        [~,xI]=sort(xleft);
        [~,xII]=sort(xI);
        v=zeros(k,k);
        for q=1:k-1        
            v(xII(q),xII(q+1))=1;
        end
        v(xII(k),xII(1))=1;
        fleft=f(v,dis);   % �������ζǿ��
        % ���������������������ϵ���Ӧ��ֵ  
        [~,xI]=sort(xright);
        [~,xII]=sort(xI);
        v=zeros(k,k);
        for q=1:k-1        
            v(xII(q),xII(q+1))=1;
        end
        v(xII(k),xII(1))=1;
        fright=f(v,dis); % �������ζǿ��
        % ȷ����һ����λ��
        x=x-step*dir.*sign(fleft-fright); % fleft>fright����sign()=1���������뷽���н�����step��fleft<fright����sign()=-1���������뷽���н�����step
        x_store(:,i)=x;
        [~,xI]=sort(x);    
        [~,xII]=sort(xI); 
        v=zeros(k,k);     
        for q=1:k-1        
            v(xII(q),xII(q+1))=1;
        end
        v(xII(k),xII(1))=1;
        fnew=f(v,dis);          % ���ݼ����x����Ҫ�Ż��ĺ���f
        % ������ε����У�fֵ��С����xII��Ϊxbest��f��Ϊfbest
        if fnew<fbest 
            xbest=xII;
            fbest=fnew;
            X_store(:,i)=xII; 
            f_store(i)=fbest;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %����Ŵ��㷨��ʼ����������������ϵ㣬���ڽ��з�ת����������x���������
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
                    Route(I:J,p)=sub(end:-1:1);             % ��תͻ�䣨Inversion mutation��
                case 3  
                    Route([I J],p)=sub([end 1]);            % ����ͻ�䣨Exchange mutation��            
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
        %����Ŵ��㷨����
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %��m����ţ��ѡ������ֵ
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
%�������ֽ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �㷨����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


