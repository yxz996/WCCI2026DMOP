classdef RJY6 < PROBLEM  %PS变PF变
% <multi> <real> <large/none> <dynamic>
% Benchmark dynamic MOP proposed by Farina, Deb, and Amato
% taut --- 20 --- Number of generations for static optimization
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
        taut;       % Number of generations for static optimization
        nt;         % Number of distinct steps
        Optimums;   % Point sets on all Pareto fronts
        lower_limit; % the lower limit of random sequence
        upper_limit; % the upper limit of random sequence
        random_integers; % random sequence
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.taut,obj.nt] = obj.ParameterSet(20,10);
            obj.M = 2;
            if isempty(obj.D); obj.D = 10; end
            obj.lower    = [0,-ones(1,obj.D-1)];
            obj.upper    = [1, ones(1,obj.D-1)];
            obj.encoding = ones(1,obj.D);
            rng(42);
            obj.lower_limit = 0;
            obj.upper_limit = 2*obj.nt + 1;
            obj.random_integers = randperm(obj.upper_limit - obj.lower_limit+1, obj.upper_limit) + obj.lower_limit - 1;
        end
        %% Evaluate solutions
        function Population = Evaluation(obj,varargin)
            PopDec     = obj.CalDec(varargin{1});
            PopObj     = obj.CalObj(PopDec);
            PopCon     = obj.CalCon(PopDec);%约束默认为0，除非重构函数
            % Attach the current number of function evaluations to solutions
            Population = SOLUTION(PopDec,PopObj,PopCon,zeros(size(PopDec,1),1)+obj.FE);
            obj.FE     = obj.FE + length(Population);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            %t = floor(obj.FE/obj.N/obj.taut)/obj.nt;
            index_t = mod(floor(obj.FE/obj.N/obj.taut), obj.upper_limit);
            Q_t = obj.random_integers(:,index_t + 1);
            t = Q_t/obj.nt;
            G = sin(0.5 * pi * t);
            A = 0.1;
            W = 3;
            K = 2*floor(10*abs(G));
            y = PopDec(:,2:end)- G;
            
            g = 1 + sum((4*y(:,2:end).^2 - cos(K*pi*y(:,2:end))+1),2);
            PopObj(:,1) = g.*( PopDec(:,1) + A*sin(W*pi.*PopDec(:,1)));
            PopObj(:,2) = g.*(1-PopDec(:,1) + A*sin(W*pi.*PopDec(:,1)));
        end
        %% Generate points on the Pareto front 画出truePF
        function R = GetOptimum(obj,N)
            %tt = floor((0:obj.maxFE)/obj.N/obj.taut)/obj.nt;
            index_t = mod(floor(0:obj.maxFE/obj.N/obj.taut), obj.upper_limit);
            Q_t = obj.random_integers(:,index_t + 1);
            tt = Q_t/obj.nt;
            H = sin(0.1.*pi.*tt);
            H = round(H*1e6)/1e6;
            x = linspace(0,1,N)';

            obj.Optimums = {};
            for i = 1 : length(H)
                A=0.1;
                R=[];
                R(:,1)=x + A*sin(3*pi*x);
                R(:,2)=1 - x + A*sin(3*pi*x);
                obj.Optimums(i,:) = {H(i),R};
            end
            % Combine all point sets
            R = cat(1,obj.Optimums{:,2});
        end
        %% Generate the image of Pareto front PF不变的情况下就这么做
        function R = GetPF(obj)
            R = obj.GetOptimum(1000);
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            t  = floor(0:obj.maxFE/obj.N/obj.taut)/obj.nt;
            Scores = zeros(1,length(t));
%              for i = 1:obj.N:length(Population)
%                 subPop = Population(i:i + obj.N - 1);
%                 k = ceil(i/obj.N);
%                 Scores(k) = feval(metName,subPop,obj.Optimums{k, 2});
%             end
            for i = 0:obj.N: (length(Population) - obj.N)
                if i ==0
                    subPop = Population(i+1:i + obj.N);
                    k = 1;
                    Scores(k) = feval(metName,subPop,obj.Optimums{k, 2});
                else
                    subPop = Population(i+1:i + obj.N);
                    k = floor((i + obj.N)/obj.N);
                    Scores(k) = feval(metName,subPop,obj.Optimums{k, 2});
                end
            end
            %   Score_IGD=[log10(Scores)]; %存MIGD
            %    save("C:\Users\zhang\Desktop\动态偏好\transfer-dynamics\code\PlatEMO-master\PlatEMO-master\PlatEMO\IGD\KAEPDF3"+".mat","Score_IGD");
            score = mean(Scores);
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            t      = floor(Population.adds/obj.N/obj.taut)/obj.nt;%当迭代这里t=0时，其实对于测试问题的t=1 (只不过这里的t变化几次，那里的t也变化几次)
            H      = sin(0.1.*pi.*t);
            H      = round(H*1e6)/1e6;
            change = [0;find(H(1:end-1)~=H(2:end));length(H)];%当t没变的时候，change就是[0;100],当变了的时候，length(change)就大于2了(就是为了后面画出每次变化后对应的个体)
            allH   = cell2mat(obj.Optimums(:,1));
            tempStream = RandStream('mlfg6331_64','Seed',2);
            for i = 1 : length(change)-1
                color = rand(tempStream,1,3);
                showdata=Population(change(i)+1:change(i+1)).objs+(i-1)*0.1;%因为PF不变，所以这里要每一维都要加个步长来区分
                Draw(showdata,'o','MarkerSize',5,'Marker','o','Markerfacecolor',sqrt(color),'Markeredgecolor',color,{'\it f\rm_1','\it f\rm_2',[]});
                drawPF=obj.Optimums{find(H(change(i)+1)==allH,1),2}+(i-1)*0.1;%加了一个步长
%                 Draw(Population.objs,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
                Draw(drawPF,'-','LineWidth',1,'Color',color);%PF不变，故truePF要加一定步长来区分
            end
        end
    end
end