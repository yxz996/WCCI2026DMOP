classdef RDCP9 < PROBLEM %
% <multi> <real> <large/none> <dynamic> <constrained>
% taut --- 10 --- Number of generations for static optimization
% nt   --- 10 --- Number of distinct steps

%------------------------------- Reference --------------------------------
% M. Farina, K. Deb, and P. Amato, Dynamic multiobjective optimization
% problems: Test cases, approximations, and applications, IEEE Transactions
% on Evolutionary Computation, 2004, 8(5): 425-442.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    properties
        Optimums={};   % Point sets on all Pareto fronts
        output1=NaN(1,200);
        lower_limit; % the lower limit of random sequence
        upper_limit; % the upper limit of random sequence
        random_integers; % random sequence
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.taut,obj.nt] = obj.ParameterSet(10,10);
            obj.M = 2;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
            rng(42);
            obj.lower_limit = 0;
            obj.upper_limit = floor(obj.maxFE/(obj.N/obj.taut))+ 1; % number of environments = number of changes 
            obj.random_integers = randperm(obj.upper_limit - obj.lower_limit, obj.upper_limit) + obj.lower_limit+1;
        end
        %% Evaluate solutions
        function Population = Evaluation(obj,varargin)
            PopDec     = obj.CalDec(varargin{1});
            PopObj     = obj.CalObj(PopDec);
            PopCon     = Constraint(obj,PopObj);
            % Attach the current number of function evaluations to solutions
            Population = SOLUTION(PopDec,PopObj,PopCon,zeros(size(PopDec,1),1)+obj.FE);%这里把zeros(size(PopDec,1),1)+obj.FE加到了每个解的add属性(我这里改成了加代数)
            obj.FE     = obj.FE + length(Population);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            index_t = mod(floor(obj.FE/obj.N/obj.taut), obj.upper_limit);
            Q_t = obj.random_integers(:,index_t + 1);
            t = Q_t/obj.nt;
            G = abs(sin(0.5*pi*t));
            r = 1 + floor((obj.D-1)*G);
            UnDec = PopDec;
            UnDec(:,r) = [];
            g = 1 + 10*sum((UnDec-G).^2,2);
            PopObj(:,1) =g.*PopDec(:,r);
            PopObj(:,2) =g.*(1-PopDec(:,r)) ;
        end

        %% Generate points on the Pareto front  
        function R = GetOptimum(obj,N) 
            index_t = mod(floor(0:obj.maxFE/obj.N/obj.taut), obj.upper_limit);
            Q_t = obj.random_integers(:,index_t + 1);
            t = Q_t/obj.nt;
            tt = t;
            H = sin(0.01.*pi.*tt);
            H = unique(round(H*1e6)/1e6);

            x1 = (0:1/(501-1):1)';
            X = UniformPoint(500,2);
            pf1 = X./repmat(sqrt(sum(X.^2,2)),1,2);

            for i = 1 : length(H)
                pf=[];
                t = (i-1) / obj.nt;
                G = sin(0.5*pi*t);
                pf(:,1) = x1 ;
                pf(:,2) = 1-x1;
                c1 = (pf(:,1)).^2+ (pf(:,2)).^2 - ((0.2+abs(G))).^2;
                pf(c1<0,:) = [];
                pf = [pf;(0.2+abs(G))*pf1];
                pf(NDSort(-pf,1)~=1,:) = [];

    
                obj.Optimums(i,:) = {H(i),pf};
            end
            R=cat(1,obj.Optimums{:,2});
        end

        %% Generate the image of Pareto front  
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end


        %% Display a population in the objective space 
        function DrawObj(obj,Population) 
                t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
                H      = sin(0.01.*pi.*t);
                H      = round(H*1e6)/1e6;
                change = [0;find(H(1:end-1)~=H(2:end));length(H)];
                allH   = cell2mat(obj.Optimums(:,1));
                tempStream = RandStream('mlfg6331_64','Seed',2);
                for i = 1 : length(change)-1
                    color = rand(tempStream,1,3);
                    showdata=Population(change(i)+1:change(i+1)).objs+(i-1)*0.5;
                    Draw(showdata,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                    ax=Draw([],'-','LineWidth',1,'Color',color);
                    drawPF=obj.Optimums{find(H(change(i)+1)==allH,1),2}+(i-1)*0.5;
                    Draw(drawPF,'-','LineWidth',1,'Color',[0 0 0]);
                end
        end

        %% Calculate the metric value  
        function score = CalMetric(obj,metName,Population)
                t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;
                H      = sin(0.01.*pi.*t);
                H      = round(H*1e6)/1e6;
                tt = int32(unique(Population.adds/100));
                tt=tt(end);
                change = [0;find(H(1:end-1)~=H(2:end));length(H)];
                Scores = zeros(1,length(change)-1);
                allH   = cell2mat(obj.Optimums(:,1));
                for i = 1 : length(change)-1
                    subPop    = Population(change(i)+1:change(i+1));
                    True_PF = obj.Optimums{i,2};
                    Scores(i) = feval(metName,subPop,True_PF);
                end
                score = mean(Scores);

        end

    end
end

%% Calculate constraint violations values
function PopCon = Constraint(obj,PopObj)
    index_t = mod(floor(obj.FE/obj.N/obj.taut), obj.upper_limit);
    Q_t = obj.random_integers(:,index_t + 1);
    t = Q_t/obj.nt;
    G = abs(sin(0.5*pi*t));
    c11 = PopObj(:,1).^2+ PopObj(:,2).^2 - (0.2+G).^2;
    c12 = PopObj(:,1).^(0.75+1.25*G) + PopObj(:,2).^(0.75+1.25*G) - 4.^(0.75+1.25*G);
    PopCon(:,1) = c11.*c12;
    c21 = PopObj(:,1).^2+ PopObj(:,2).^2 - 1.6.^2;
    c22 = 2.1 - (PopObj(:,1)./(1+0.15*cos(6*atan(PopObj(:,2)./PopObj(:,1)).^3).^10)).^2 - (PopObj(:,2)./(1+0.75*cos(6*atan(PopObj(:,2)./PopObj(:,1)).^3).^10)).^2;
    PopCon(:,2) = c21.*c22;
end