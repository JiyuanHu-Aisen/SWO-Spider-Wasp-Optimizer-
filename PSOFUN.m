function [ bestPosition, fitValue,BestCost ] = ...
PSOFUN( CostFun,nVar,VarMin,VarMax,MaxIt,nPop )
%% PSO Parameters 
CostFunction=@(x) CostFun(x);        % Cost Function
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient
VarSize=[1 nVar];   % Size of Decision Variables Matrix
% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    
    % Initialize Position
    particle(i).Position=round(unifrnd(VarMin,VarMax,VarSize));%�����ʼ������
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);%�ٶ�
    
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position);%Ŀ��
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;%���¸�������
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best%����ȫ������
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

%% PSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity�����ٶ�
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits�ٶ�Լ��
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position����λ��
        particle(i).Position = round(particle(i).Position + particle(i).Velocity);
        
        % Velocity Mirror Effect����Խ�޾�����
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits%����Խ�޴���
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation��Ŀ��
        particle(i).Cost = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best����ȫ������
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    w=w*wdamp;
    
end

bestPosition = GlobalBest.Position;%��¼����λ��
fitValue = GlobalBest.Cost;%��¼����Ŀ����

end

