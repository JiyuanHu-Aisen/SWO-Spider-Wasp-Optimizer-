%% 蛛蜂算法（SWO）(整数取值版)
function [SW_Best,curve_draw,t]=SWOFUN1(Target_Fun,N,Nmin,CR,TR,tmax,lower_lim,upper_lim)
% 输入变量分别为：目标函数、种群规模、种群最小值、交叉概率、狩猎和交配行为之间的权衡概率、最大步长
% 最佳的 TR=0.3 CR=0.2 Nmin=20 N=100
%% 输入数据的读取
Target_Function=@(x) Target_Fun(x);
dim=size(lower_lim,2); % 返回的是元素个数，也即需要每个Wasp表示的维度
%% 1.Generation of the initial population
SWpop=repmat(lower_lim,N,1)+rand(N,dim).*(upper_lim-lower_lim);

SWpop=round(SWpop);

Best_Cost=inf;
for i=1:N
    T0=Target_Function(SWpop(i,:));
    if T0<Best_Cost
        Best_Cost=T0;
        SW_Best=SWpop(i,:);
    end
end
%% 2+3.Hunting and nesting behavior/Mating behavior
t=0;
curve_draw=zeros(1,tmax);
while(t<tmax)
    N_origin=N;
    r6=rand(1);
    k=1-t/tmax;
    SW_Best_Past=SWpop;
    %% 2.Hunting and nesting behavior
    if(r6<TR) % to fig.4
        for i=1:N           
            r3=rand(1);
            r4=rand(1);
            l=unifrnd(-1,2);
            %% Searching stage(exploration)/Following and escaping stage (exploration and exploitation)
            if(i<k*N)
                p=rand(1);
                %% 2.1Searching stage(exploration)
                if(p<k)
                    if(r3<r4)
                        rn=randn(1);
                        r1=rand(1);
                        miu1=abs(rn)*r1;
                        tworand=randperm(N,2);
                        SWpop(i,:)=SWpop(i,:)+miu1*(SWpop(tworand(1,1),:)-SWpop(tworand(1,2),:)); % eq(4)
                    else

                        B=1/(1+exp(l));
                        miu2=B*cos(2*pi*l);
                        c=randperm(N,1);
                        rr2=rand(1,dim);
                        SWpop(i,:)=SWpop(c,:)+miu2*(lower_lim+rr2.*(upper_lim-lower_lim)); % eq(6)
                    end

                %% 2.2 Following and escaping stage (exploration and exploitation)
                else
                    if(r3<r4)
                        C=(2-2*(t/tmax))*r6;
                        a=randperm(N,1);
                        rr5=rand(1,dim);
                        SWpop(i,:)=SWpop(i,:)+C*abs(2*rr5.*SWpop(a,:)-SWpop(i,:)); % eq(10)
                    else
                        vc=unifrnd(-k,k,1,dim);
                        SWpop(i,:)=SWpop(i,:).*vc; % eq(12)
                    end
                end

            %% 2.3 Nesting behavior (exploitation)
            else
                if(r3<r4)
                    SWpop(i,:)=SW_Best+cos(2*pi*l)*(SW_Best-SWpop(i,:)); % eq(16)
                else
                    rr4=rand(1,dim);
                    rr5=rand(1,dim);
                    % 生成需要的U向量
                    Upanju=(rr4>rr5);
                    U=zeros(1,dim);
                    U(Upanju)=1;
                    %***************%
                    a=randperm(N,1);
                    tworand=randperm(N,2);
                    gama=levy_flight;
                    SWpop(i,:)=SWpop(i,:)+r3*abs(gama)*(SWpop(a,:)-SWpop(i,:))+(1-r3)*U.*(SWpop(tworand(1,1),:)-SWpop(tworand(1,2),:)); % eq(17)
                end
            end

        end

    %% 3.Mating behavior
    else
        origin_sequence=linspace(1,N,N);
        for i=1:N
            l=unifrnd(-2,1);
            handle_sequence=origin_sequence(origin_sequence~=i);
            threerand=handle_sequence(randperm(N-1,3)); % 分别是a,b,c

            if Target_Function(SWpop(threerand(1),:))<Target_Function(SWpop(i,:)) % eq(23)
                vector1=SWpop(threerand(1),:)-SWpop(i,:);
            else
                vector1=-SWpop(threerand(1),:)+SWpop(i,:);
            end

            if Target_Function(SWpop(threerand(2),:))<Target_Function(SWpop(threerand(3),:)) % eq(24)
                vector2=SWpop(threerand(2),:)-SWpop(threerand(3),:);
            else
                vector2=-SWpop(threerand(2),:)+SWpop(threerand(3),:);
            end

            beta0=randn(1);
            beta1=randn(1);
            SWmale=SWpop(i,:)+exp(l)*abs(beta0)*vector1+(1-exp(l))*abs(beta1)*vector2; % eq(22)

            SWmale=round(SWmale);
            SWmale=max(SWmale,lower_lim);
            SWmale=min(SWmale,upper_lim);


            %% Crossover（交叉）
            %假定左父右母
            cutpoint=randperm(dim,1);
            if (rand(1)<CR) % eq(21)
                if(cutpoint==dim)
                    SWpop(i,:)=SWmale;
                elseif(cutpoint~=0)
                    SWpop(i,1:cutpoint)=SWmale(1,1:cutpoint);
                else
                end
            else
            end

        end
    end
%% 4.Population reduction and memory saving
%% 4.1 Memory saving
    %% 修正越限的情况以及选优
    for i=1:N
        SWpop(i,:)=max(SWpop(i,:),lower_lim);
        SWpop(i,:)=min(SWpop(i,:),upper_lim);

        SWpop=round(SWpop);

        T1=Target_Function(SWpop(i,:));
        if(T1<Best_Cost)
            Best_Cost=T1;
            SW_Best=SWpop(i,:);
        elseif(T1 > Target_Function(SW_Best_Past(i,:)) )
            SWpop(i,:)=SW_Best_Past(i,:);
        end

    end
    t=t+1;
    disp(['Iteration ' num2str(t) ': Best Cost = ' num2str(Best_Cost)]);
    curve_draw(t)=Best_Cost;

%% 4.2 Population reduction
    N=fix(Nmin+(N-Nmin)*k);

end

end

%% levy_flight生成随机步长
function gama=levy_flight
beta = 3/2;
alpha_u = ((gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2))))^(1/beta);
alpha_v = 1;
u = alpha_u*randn(1);
v = alpha_v*randn(1);
gama = u/(abs(v))^(1/beta);
end

