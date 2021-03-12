function [MeanNumOfHealthy,MeanTotalInfections,MeanNumberOfSympt...
    ,MeanNumberOfAsympt,initTestingDay,MeanNumOfRecoveries,...
    MeanNumberOfCasesDetected,MeanProportionOfInfectedIsolations...
    ,MeanTotalAbsentDays,MeanNumberOfIsolatingRecovered...
    ,MeanNumberOfIsolatingHealthy] = General_Delay_LFT(Para)


%We now assume that simulations start on a sunday and tests only start with
%and initial initial cases could be either symp to asymp and with normal
%infect,detect,recover timers.

%We assume a symptomatic shows signs of infection after tInf days, the same
%time the student becomes infectious

%Isolate the whole class on infection: SubgroupSize = Ntot
%Isolate the individual on infection: SubgroupSize = 1
%Isolate the group on infection: SubgroupSize in [2,Ntot-1] and divisor

%ActivePop key: 0=susceptible, 1=infected (asympt), 2=recovered, NaN=isolating,
%3=isolating but healthy, 4 = infected (sympt), 5 = isolating but recovered
%% Simulation parameters

Ntot = Para(1); %total pop. size
SubgroupSize = Para(2); %number of individuals in each group
Irb = Para(3);%Background prevalence
SocialDistancingModifier = Para(17);%for FE sims
R = Para(4)/SocialDistancingModifier; %R number
Cf = Para(5); %far R ratio
Cc = Para(6); %close R ratio
PerFalseNeg = Para(7); %percentage of false negatives from tests
CompIso = Para(8); %compliance with isolation for any individual
MixInd = Para(9);%see mixing days
TestInd = Para(10); %this now goes into the changeTestingSchedule fucntionNg
PrevBool = Para(11);
RecoveryBool = Para(12); %bool
PerSymptomatic = Para(13); %percentage of symptomatic
RemoveSympt = Para(14); %Boolean to remove those that show symptoms
PerFalsePos = Para(15); %percentage of false positives
InitWithTesting = Para(16); %start simulations with testing days active - false means to wait for sympto

%Local infection prevalence of symptomatic (Percentage)
Ir =  Irb;

%if the subgroups are individuals, only infect others
if SubgroupSize == 1
    Cc = 0;
elseif SubgroupSize == Ntot
    Cf = 0;
end

%heterogeneous transmission
Rinter =R*(Cf/(Cc+Cf));
Rintra = R*(Cc/(Cc+Cf));

%infectious time (days) - should probably set this as a distribution for each
%agent at some point (also get data for this)
tInf = 3;

%detectable time (days)
tDet = 5;

%Recovery time (days)
tRec = 10;

%number of simulations
NoS = 1e3;

%intiail testing days
TestingDays = zeros(1,28);

if MixInd == 1
    %school days starting on a sunday
    MixingDays =  [0,1,1,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,1,1,0];
elseif MixInd == 2
    %HE - continuous mixing 
    MixingDays =  ones(1,28);
elseif MixInd == 3
    %HE - Arrival scenario
    MixingDays =  [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
end

if InitWithTesting == 1
    TestingDays = ChangeTestingSchedule(InitWithTesting,1,TestInd);
    SymptFlag = 1;
end

%number of uniform groups in the pop.
Groups = Ntot/SubgroupSize;

%checking for errors
if length(TestingDays) ~= length(MixingDays)
    error("EXPECTION: Testing days and agent mixing days must be the same length")
end

if floor(Ntot/SubgroupSize) ~= Ntot/SubgroupSize
    error('EXCEPTION: subgroup size must be a divisor of pop size')
end


%% simulation

%ActivePop key: 0=susceptible, 1=infected (asympt), 2=recovered, NaN=isolating,
%3=isolating but healthy, 4 = infected (sympt), 5 = isolating but recovered

Ir = Ir/100;
Irb = Irb/100;
PerFalseNeg = PerFalseNeg/100;
CompIso = CompIso/100;
PerSymptomatic = PerSymptomatic/100;
PerFalsePos = PerFalsePos/100;

NumberOfInfections = zeros(NoS,numel(MixingDays));
MeanNumOfInfections = zeros(1,numel(MixingDays));
stdNumOfInfections = zeros(1,numel(MixingDays));

HealthyPeople = zeros(NoS,numel(MixingDays));
MeanNumOfHealthy = zeros(1,numel(MixingDays));
stdNumOfHealthy = zeros(1,numel(MixingDays));

NumberOfSympt = zeros(NoS,numel(MixingDays));
MeanNumberOfSympt = zeros(1,numel(MixingDays));
stdNumberOfSympt = zeros(1,numel(MixingDays));

NumberOfAsympt = zeros(NoS,numel(MixingDays));
MeanNumberOfAsympt = zeros(1,numel(MixingDays));
stdNumberOfAsympt = zeros(1,numel(MixingDays));

NumberOfRecoveries = zeros(NoS,numel(MixingDays));
MeanNumberOfRecoveries = zeros(1,numel(MixingDays));
stdNumberOfRecoveries = zeros(1,numel(MixingDays));

NumberOfCasesDetected = zeros(NoS,numel(MixingDays));
MeanNumberOfCasesDetected = zeros(1,numel(MixingDays));
stdNumberOfCasesDetected = zeros(1,numel(MixingDays));

TotalOutbreakCounter = zeros(1,NoS);
TotalAbsentDays = zeros(1,NoS);

initTestingDay = zeros(1,NoS);

for i = 1:NoS
    %% Initialisation
    if InitWithTesting == 0
        SymptFlag = 0;
    end
    %initialise the number of symptomatic
    Infectedsym = rand(1,Ntot)<Ir;
    Ninfsym = max(sum(Infectedsym),1);
    
    %vectors for timers
    CountDownDetection = Inf(1,Ntot);
    CountDownInfectious = Inf(1,Ntot);
    CountDownRecovery = Inf(1,Ntot);
    
    %vector for population dynamics
    ActivePop = zeros(1,Ntot);
    initInfectedIndicies = randi([1 Ntot],1,Ninfsym);
    
    %check the each index is different
    while length(initInfectedIndicies) ~= length(unique(initInfectedIndicies))
        initInfectedIndicies = randi([1 Ntot],1,Ninfsym);
    end
    
    %randomly select an agent for infection
    for z = 1:Ninfsym
        if rand(1,1) < PerSymptomatic
            ActivePop(initInfectedIndicies(z)) = 4;
        else
            ActivePop(initInfectedIndicies(z)) = 1;
        end
        CountDownDetection(initInfectedIndicies(z)) = tDet;
        CountDownInfectious(initInfectedIndicies(z)) = tInf;
        CountDownRecovery(initInfectedIndicies(z)) = tRec;
    end
    
    
    %vector of inter and intra infections for each infected agent
    if MixingDays(1) == 1
        
        infected_indices = find(CountDownInfectious<=0);
        NumOfInfectious = sum(CountDownInfectious<=0);
        
        %split the total pop into subgroups
        Subgroups = mat2cell(ActivePop,1, SubgroupSize.*ones(1,Groups));
        SubDetectTimer = mat2cell(CountDownDetection,1, SubgroupSize.*ones(1,Groups));
        SubInfectTimer = mat2cell(CountDownInfectious,1, SubgroupSize.*ones(1,Groups));
        SubRecoveryTimer = mat2cell(CountDownRecovery,1, SubgroupSize.*ones(1,Groups));
        
        %intra-group infections
        SiReIndicesIntra = poissrnd(Rintra, NumOfInfectious,1);
        
        %inter-group infections
        SiReIndicesInter = poissrnd(Rinter, NumOfInfectious,1);
        
        %check that infections cannot be larger than 0's
        for a = 1:length(infected_indices)
            b = floor((infected_indices(a)-1)/SubgroupSize)+1;
            SiReIndicesIntra(a) = min(SiReIndicesIntra(a),sum(Subgroups{b}==0));
        end
        
        for a = 1:length(infected_indices)
            b = floor((infected_indices(a)-1)/SubgroupSize)+1;
            availGroups = 1:Groups;
            availGroups(b) = [];
            numOfSus = 0;
            for c  = 1:length(availGroups)
                numOfSus = numOfSus + sum(Subgroups{c}==0);
            end
            SiReIndicesInter(a) = min(SiReIndicesInter(a),numOfSus);
        end
        
        %final check for over inter infections
        while sum(SiReIndicesInter) > sum(ActivePop==0)
            SiReIndicesInter = SiReIndicesInter-1;
            for z = 1:length(initSiReIndicesInter)
                if SiReIndicesInter(z)<0
                    SiReIndicesInter(z) = 0;
                end
            end
        end
        
        
        %allocate the intra infections
        for v = 1:length(SiReIndicesIntra)
            b = floor((infected_indices(v)-1)/SubgroupSize)+1;
            for c = 1:SiReIndicesIntra(v)
                for d = 1:length(Subgroups{b})
                    if Subgroups{b}(d) == 0
                        if rand(1,1) < PerSymptomatic
                            Subgroups{b}(d) = 4;
                        else
                            Subgroups{b}(d) = 1;
                        end
                        SubDetectTimer{b}(d) = tDet;
                        SubInfectTimer{b}(d) = tInf;
                        SubRecoveryTimer{b}(d) = tRec;
                        break
                    end
                end
            end
        end
        
        
        %randomly allocate the inter infections
        for v = 1:length(SiReIndicesInter)
            b = floor((infected_indices(v)-1)/SubgroupSize)+1;
            availGroups = 1:Groups;
            availGroups(b) = [];
            availGroups = repmat(availGroups,1,SubgroupSize);
            availGroups = availGroups(randperm(length(availGroups)));
            availIndices = repmat(1:SubgroupSize,1,Groups-1);
            availIndices = availIndices(randperm(length(availIndices)));
            
            randomAvailIndices = [availGroups',availIndices'];
            
            for z = 1:SiReIndicesInter(v)
                for c = 1:((Groups-1)*SubgroupSize)
                    groupAndIndex = randomAvailIndices(c,:);
                    if Subgroups{groupAndIndex(1)}(groupAndIndex(2)) == 0
                        if rand(1,1) < PerSymptomatic
                            Subgroups{groupAndIndex(1)}(groupAndIndex(2)) = 4;
                        else
                            Subgroups{groupAndIndex(1)}(groupAndIndex(2)) = 1;
                        end
                        SubDetectTimer{groupAndIndex(1)}(groupAndIndex(2)) = tDet;
                        SubInfectTimer{groupAndIndex(1)}(groupAndIndex(2)) = tInf;
                        SubRecoveryTimer{groupAndIndex(1)}(groupAndIndex(2)) = tRec;
                        break
                    end
                end
            end
        end
        
        
        %update the counters and population
        ActivePop = [];
        CountDownDetection = [];
        CountDownInfectious = [];
        CountDownRecovery = [];
        
        for a = 1:Groups
            ActivePop = [ActivePop,Subgroups{a}];
            CountDownDetection = [CountDownDetection,SubDetectTimer{a}];
            CountDownInfectious = [CountDownInfectious,SubInfectTimer{a}];
            CountDownRecovery = [CountDownRecovery,SubRecoveryTimer{a}];
        end
        
    end
    
    HealthyPeople_plus_recovered(i,1)= sum(ActivePop==0) + sum(ActivePop==3)-sum((ActivePop==2)) ;
    NumberOfInfections(i,1) = sum(ActivePop==1) + sum(ActivePop==4);
    NumberOfSympt(i,1) = sum(ActivePop==4);
    NumberOfAsympt(i,1) = sum(ActivePop==1);
    
    %Day is over: remove a day from each timer
    CountDownDetection = CountDownDetection - 1;
    CountDownInfectious = CountDownInfectious - 1;
    CountDownRecovery = CountDownRecovery - 1;
    
    %if on - prevalence based infections after each day
    if PrevBool == 1
        for x = 1:length(ActivePop)
            %if agent is not infected then there is a possibility of
            %becoming infected from background
            if ActivePop(x) == 0
                if rand(1,1) < Irb
                    if rand(1,1) < PerSymptomatic
                        ActivePop(x) = 4;
                    else
                        ActivePop(x) = 1;
                    end
                    CountDownDetection(x) = tDet;
                    CountDownInfectious(x) = tInf;
                    CountDownRecovery(x) = tRec;
                end
            end
            
        end
    end
    
    %check for sympt flag
    if SymptFlag == 0
        for x = 1:length(ActivePop)
            if CountDownInfectious(x) <= 0 && ActivePop(x) == 4
                SymptFlag = 1;
                initTestingDay(i) = 1;
            end
        end
        TestingDays = ChangeTestingSchedule(SymptFlag,1,TestInd);
    end
    
    %initialisation over
    
    %% Loop over each day post initalisation day
    
    %each loop accounts for a single day
    HealthyPeople(i,1)=sum(ActivePop==0);
    for j = 2:numel(MixingDays)
        
        %if on - loop to return all recovered pop back.
        if RecoveryBool == 1
            for l = 1:length(ActivePop)
                if isnan(ActivePop(l)) && CountDownRecovery(l) <=0
                    ActivePop(l) = 2;
                    CountDownDetection(l) = Inf;
                    CountDownInfectious(l) = Inf;
                    CountDownRecovery(l) = Inf;
                    
                elseif ActivePop(l)== 1 && CountDownRecovery(l) <=0
                    ActivePop(l) = 2;
                    CountDownDetection(l) = Inf;
                    CountDownInfectious(l) = Inf;
                    CountDownRecovery(l) = Inf;
                    
                elseif ActivePop(l)== 3 && CountDownRecovery(l) <=0
                    ActivePop(l) = 0;
                    CountDownDetection(l) = Inf;
                    CountDownInfectious(l) = Inf;
                    CountDownRecovery(l) = Inf;
                    
                elseif ActivePop(l) == 4 && CountDownRecovery(l) <=0
                    ActivePop(l) = 2;
                    CountDownDetection(l) = Inf;
                    CountDownInfectious(l) = Inf;
                    CountDownRecovery(l) = Inf;
                
                elseif ActivePop(l) == 5 && CountDownRecovery(l) <=0
                    ActivePop(l) = 2;
                    CountDownDetection(l) = Inf;
                    CountDownInfectious(l) = Inf;
                    CountDownRecovery(l) = Inf;
                end
            end
        end
          
        
        %Update LFT based on agent timer
        if TestingDays(j) == 1
            posTestIndices = [];
            for z = 1:length(ActivePop)
                %if the agent is detectable then run a LFT AND infected
                if (CountDownDetection(z) <= 0 && ActivePop(z) == 1) ...
                        || (CountDownDetection(z) <= 0 && ActivePop(z) == 4)
                    if rand(1,1)>PerFalseNeg
                        if rand(1,1) < CompIso
                            %successful LFT: 'remove' agent from pop. and
                            ActivePop(z) = NaN;
                            CountDownDetection(z) = Inf;
                            CountDownInfectious(z) = Inf;
                            CountDownRecovery(z) = tRec;
                            posTestIndices = [posTestIndices ,z];
                        end
                    end
                end
                
                %Produce False positives
                if (ActivePop(z) == 0)
                    if rand(1,1)<PerFalsePos
                        if rand(1,1) < CompIso
                            %false LFT: 'remove' agent from pop. and
                            ActivePop(z) = 3;
                            CountDownDetection(z) = Inf;
                            CountDownInfectious(z) = Inf;
                            CountDownRecovery(z) = tRec;
                            posTestIndices = [posTestIndices ,z];
                        end
                    end
                end
                
                 %Produce False positives
                if (ActivePop(z) == 2)
                    if rand(1,1)<PerFalsePos
                        if rand(1,1) < CompIso
                            %false LFT: 'remove' agent from pop. and
                            ActivePop(z) = 5;
                            CountDownDetection(z) = Inf;
                            CountDownInfectious(z) = Inf;
                            CountDownRecovery(z) = tRec;
                            posTestIndices = [posTestIndices ,z];
                        end
                    end
                end
                
            end
         
            
            %split the total pop into subgroups
            Subgroups = mat2cell(ActivePop,1, SubgroupSize.*ones(1,Groups));
            SubDetectTimer = mat2cell(CountDownDetection,1, SubgroupSize.*ones(1,Groups));
            SubInfectTimer = mat2cell(CountDownInfectious,1, SubgroupSize.*ones(1,Groups));
            SubRecoveryTimer = mat2cell(CountDownRecovery,1, SubgroupSize.*ones(1,Groups));
            
            
            for a = 1:length(posTestIndices)
                b = floor((posTestIndices(a)-1)/SubgroupSize)+1;
                for c = 1:length(Subgroups{b})
                    if Subgroups{b}(c) == 1 || isnan(Subgroups{b}(c)) ...
                            || Subgroups{b}(c) == 4
                        if rand(1,1) < CompIso
                            Subgroups{b}(c) = NaN;
                            SubDetectTimer{b}(c) = Inf;
                            SubInfectTimer{b}(c) = Inf;
                            SubRecoveryTimer{b}(c) = tRec;
                        end
                    elseif (Subgroups{b}(c) == 2 || Subgroups{b}(c) == 5)
                        if rand(1,1) < CompIso
                            Subgroups{b}(c) = 5;
                            SubDetectTimer{b}(c) = Inf;
                            SubInfectTimer{b}(c) = Inf;
                            SubRecoveryTimer{b}(c) = tRec;
                        end
                    elseif (Subgroups{b}(c) == 0 || Subgroups{b}(c) == 3)
                        if rand(1,1) < CompIso
                            Subgroups{b}(c) = 3;
                            SubDetectTimer{b}(c) = Inf;
                            SubInfectTimer{b}(c) = Inf;
                            SubRecoveryTimer{b}(c) = tRec;
                        end
                    end
                end
            end
            
            %update the counters and population
            ActivePop = [];
            CountDownDetection = [];
            CountDownInfectious = [];
            CountDownRecovery = [];
            
            for a = 1:Groups
                ActivePop = [ActivePop,Subgroups{a}];
                CountDownDetection = [CountDownDetection,SubDetectTimer{a}];
                CountDownInfectious = [CountDownInfectious,SubInfectTimer{a}];
                CountDownRecovery = [CountDownRecovery,SubRecoveryTimer{a}];
            end
        end
        
        
        
       %update secondary infections
        if MixingDays(j) == 1
            
            infected_indices = find(CountDownInfectious<=0);
            NumOfInfectious = sum(CountDownInfectious<=0);
            
            %split the total pop into subgroups
            Subgroups = mat2cell(ActivePop,1, SubgroupSize.*ones(1,Groups));
            SubDetectTimer = mat2cell(CountDownDetection,1, SubgroupSize.*ones(1,Groups));
            SubInfectTimer = mat2cell(CountDownInfectious,1, SubgroupSize.*ones(1,Groups));
            SubRecoveryTimer = mat2cell(CountDownRecovery,1, SubgroupSize.*ones(1,Groups));
            
            %intra-group infections
            SiReIndicesIntra = poissrnd(Rintra, NumOfInfectious,1);
            
            %inter-group infections
            SiReIndicesInter = poissrnd(Rinter, NumOfInfectious,1);
            
            %check that infections cannot be larger than the numnber of
            %sus.
            for a = 1:length(infected_indices)
                b = floor((infected_indices(a)-1)/SubgroupSize)+1;
                SiReIndicesIntra(a) = min(SiReIndicesIntra(a),sum(Subgroups{b}==0));
            end
            
            for a = 1:length(infected_indices)
                b = floor((infected_indices(a)-1)/SubgroupSize)+1;
                availGroups = 1:Groups;
                availGroups(b) = [];
                numOfSus = 0;
                for c  = 1:length(availGroups)
                    numOfSus = numOfSus + sum(Subgroups{c}==0);
                end
                SiReIndicesInter(a) = min(SiReIndicesInter(a),numOfSus);
            end
            
            %final check for over inter infections
            while sum(SiReIndicesInter) > sum(ActivePop==0)
                SiReIndicesInter = SiReIndicesInter-1;
                for z = 1:length(SiReIndicesInter)
                    if SiReIndicesInter(z)<0
                        SiReIndicesInter(z) = 0;
                    end
                end
            end
            
            
            %allocate the intra infections
            for v = 1:length(SiReIndicesIntra)
                b = floor((infected_indices(v)-1)/SubgroupSize)+1;
                for c = 1:SiReIndicesIntra(v)
                    for d = 1:length(Subgroups{b})
                        if Subgroups{b}(d) == 0
                            if rand(1,1) < PerSymptomatic
                                Subgroups{b}(d) = 4;
                            else
                                Subgroups{b}(d) = 1;
                            end
                            SubDetectTimer{b}(d) = tDet;
                            SubInfectTimer{b}(d) = tInf;
                            SubRecoveryTimer{b}(d) = tRec;
                            break
                        end
                    end
                end
            end
            
            %randomly allocate the inter infections
            for v = 1:length(SiReIndicesInter)
                b = floor((infected_indices(v)-1)/SubgroupSize)+1;
                availGroups = 1:Groups;
                availGroups(b) = [];
                availGroups = repmat(availGroups,1,SubgroupSize);
                availGroups = availGroups(randperm(length(availGroups)));
                availIndices = repmat(1:SubgroupSize,1,Groups-1);
                availIndices = availIndices(randperm(length(availIndices)));
                
                randomAvailIndices = [availGroups',availIndices'];
                
                for z = 1:SiReIndicesInter(v)
                    for c = 1:((Groups-1)*SubgroupSize)
                        groupAndIndex = randomAvailIndices(c,:);
                        
                        if Subgroups{groupAndIndex(1)}(groupAndIndex(2)) == 0
                            if rand(1,1) < PerSymptomatic
                                Subgroups{groupAndIndex(1)}(groupAndIndex(2)) = 4;
                            else
                                Subgroups{groupAndIndex(1)}(groupAndIndex(2)) = 1;
                            end
                            
                            SubDetectTimer{groupAndIndex(1)}(groupAndIndex(2)) = tDet;
                            SubInfectTimer{groupAndIndex(1)}(groupAndIndex(2)) = tInf;
                            SubRecoveryTimer{groupAndIndex(1)}(groupAndIndex(2)) = tRec;
                            break
                        end
                    end
                end
            end
            
            
            %update the counters and population
            ActivePop = [];
            CountDownDetection = [];
            CountDownInfectious = [];
            CountDownRecovery = [];
            
            for a = 1:Groups
                ActivePop = [ActivePop,Subgroups{a}];
                CountDownDetection = [CountDownDetection,SubDetectTimer{a}];
                CountDownInfectious = [CountDownInfectious,SubInfectTimer{a}];
                CountDownRecovery = [CountDownRecovery,SubRecoveryTimer{a}];
            end
            
        end
        
        %check for sympt flag - to init testing days
        if SymptFlag == 0
            for x = 1:length(ActivePop)
                if CountDownInfectious(x) <= 0 && ActivePop(x) == 4
                    SymptFlag = 1;
                    initTestingDay(i) = j;
                end
            end
            TestingDays = ChangeTestingSchedule(SymptFlag,j,TestInd);
        end
        
        
        %remove the symptomatics and group showing showing infection
        if RemoveSympt == 1
            posTestIndices = [];
            for z = 1:length(ActivePop)
                if ActivePop(z) == 4 && CountDownInfectious(z) <=0
                    if rand(1,1) < CompIso
                        %successful LFT: 'remove' agent from pop. and
                        ActivePop(z) = NaN;
                        CountDownDetection(z) = Inf;
                        CountDownInfectious(z) = Inf;
                        CountDownRecovery(z) = tRec;
                        posTestIndices = [posTestIndices ,z];
                    end
                end
            end
            
            %split the total pop into subgroups
            Subgroups = mat2cell(ActivePop,1, SubgroupSize.*ones(1,Groups));
            SubDetectTimer = mat2cell(CountDownDetection,1, SubgroupSize.*ones(1,Groups));
            SubInfectTimer = mat2cell(CountDownInfectious,1, SubgroupSize.*ones(1,Groups));
            SubRecoveryTimer = mat2cell(CountDownRecovery,1, SubgroupSize.*ones(1,Groups));
            
            
            for a = 1:length(posTestIndices)
                b = floor((posTestIndices(a)-1)/SubgroupSize)+1;
                for c = 1:length(Subgroups{b})
                    if Subgroups{b}(c) == 1 || isnan(Subgroups{b}(c)) ...
                            || Subgroups{b}(c) == 4
                        if rand(1,1) < CompIso
                            Subgroups{b}(c) = NaN;
                            SubDetectTimer{b}(c) = Inf;
                            SubInfectTimer{b}(c) = Inf;
                            SubRecoveryTimer{b}(c) = tRec;
                        end
                    elseif (Subgroups{b}(c) == 2 || Subgroups{b}(c) == 5)
                        if rand(1,1) < CompIso
                            Subgroups{b}(c) = 5;
                            SubDetectTimer{b}(c) = Inf;
                            SubInfectTimer{b}(c) = Inf;
                            SubRecoveryTimer{b}(c) = tRec;
                        end
                    elseif (Subgroups{b}(c) == 0 || Subgroups{b}(c) == 3)
                        if rand(1,1) < CompIso
                            Subgroups{b}(c) = 3;
                            SubDetectTimer{b}(c) = Inf;
                            SubInfectTimer{b}(c) = Inf;
                            SubRecoveryTimer{b}(c) = tRec;
                        end
                    end
                end
            end
            
            %update the counters and population
            ActivePop = [];
            CountDownDetection = [];
            CountDownInfectious = [];
            CountDownRecovery = [];
            
            for a = 1:Groups
                ActivePop = [ActivePop,Subgroups{a}];
                CountDownDetection = [CountDownDetection,SubDetectTimer{a}];
                CountDownInfectious = [CountDownInfectious,SubInfectTimer{a}];
                CountDownRecovery = [CountDownRecovery,SubRecoveryTimer{a}];
            end
        end
        
        
        %record infections
        NumberOfInfections(i,j) = sum((ActivePop==1)) + sum((ActivePop==4)) ;
        NumberOfRecoveries(i,j) = sum((ActivePop==2)) +sum(ActivePop==5)  ;
        NumberOfCasesDetected(i,j) = sum(isnan(ActivePop))+sum(ActivePop==3) +sum(ActivePop==5) ;
        HealthyPeople(i,j)= sum(ActivePop==0) + sum(ActivePop==3) ;
        NumberOfIsolatingRecovered(i,j) = sum(ActivePop==5);
        NumberOfIsolatingHealthy(i,j) = sum(ActivePop==3);
      

        NumberOfSympt(i,j) = sum(ActivePop==4);
        NumberOfAsympt(i,j) = sum(ActivePop==1);
        
        IsolatingInfections(i,j) = sum(isnan(ActivePop));
        IsolatingHealthy(i,j) = sum(ActivePop==3);
        
        if mod(j,7) ~=1 || mod(j,7) ~= 0
            TotalAbsentDays(i) = TotalAbsentDays(i) + NumberOfCasesDetected(i,j);
        end
        
        %Day is over: remove a day from each timer
        CountDownDetection = CountDownDetection - 1;
        CountDownInfectious = CountDownInfectious - 1;
        CountDownRecovery = CountDownRecovery - 1;
        
        
        %if on - prevalence based infections after each day
        if PrevBool == 1
            for x = 1:length(ActivePop)
                %if agent is not infected then there is a possibility of
                %becoming infected from background
                if ActivePop(x) == 0
                    if rand(1,1) < Irb
                        if rand(1,1) < PerSymptomatic
                            ActivePop(x) = 4;
                        else
                            ActivePop(x) = 1;
                        end
                        CountDownDetection(x) = tDet;
                        CountDownInfectious(x) = tInf;
                        CountDownRecovery(x) = tRec;
                    end
                end
            end
        end
        
    end
    
    if sum(ActivePop==0) == 0
        TotalOutbreakCounter(i) = 1;
    end
    
    % error checking for total infections
    for j = 2:length(MixingDays)
        if  HealthyPeople(i,j-1) - HealthyPeople(i,j) <0
            error("healthy switching occured ")
        end
    end
    
    
end


%% Collect Data

for i = 1:numel(MixingDays)
    
    MeanNumOfInfections(i) = mean(NumberOfInfections(:,i),'omitnan');
    stdNumOfInfections(i) = std(NumberOfInfections(:,i),'omitnan');
    
    MeanNumOfRecoveries(i) = mean(NumberOfRecoveries(:,i),'omitnan');
    stdNumOfRecoveries(i) = std(NumberOfRecoveries(:,i),'omitnan');
    
    MeanNumberOfSympt(i) = mean(NumberOfSympt(:,i),'omitnan');
    stdNumberOfSympt(i) = std(NumberOfSympt(:,i),'omitnan');
    
    MeanNumberOfAsympt(i) = mean(NumberOfAsympt(:,i),'omitnan');
    stdNumberOfAsympt(i) = std(NumberOfAsympt(:,i),'omitnan');
    
    MeanNumOfHealthy(i) = mean(HealthyPeople(:,i),'omitnan');
    stdNumOfHealthy(i) = std(HealthyPeople(:,i),'omitnan');
    
    MeanNumberOfCasesDetected(i) = mean(NumberOfCasesDetected(:,i));
    stdNumberOfCasesDetected(i) = std(NumberOfCasesDetected(:,i));
    
    MeanNumberOfIsolatingRecovered(i) = mean(NumberOfIsolatingRecovered(:,i));
    MeanNumberOfIsolatingHealthy(i) = mean(NumberOfIsolatingHealthy(:,i));
    
end


ConLower=min(MeanNumOfInfections,1.96.*stdNumOfInfections);
ConUpper=min(Ntot-MeanNumOfInfections,1.96.*stdNumOfInfections);

ConLowerRec=min(MeanNumOfRecoveries,1.96.*stdNumOfRecoveries);
ConUpperRec=min(Ntot-MeanNumOfRecoveries,1.96.*stdNumOfRecoveries);

ConLowerDet=max(MeanNumberOfCasesDetected,1.96*stdNumberOfCasesDetected);
ConUpperDet=min(Ntot-MeanNumberOfCasesDetected,1.96*stdNumberOfCasesDetected);

PercentageOfSimOutbreaks = (sum(TotalOutbreakCounter)/NoS)*100;

MeanTotalInfections=mean(Ntot-HealthyPeople);
stdTotalInfections=std(Ntot-HealthyPeople);

TotalInfectionsLower=min(MeanTotalInfections,1.96*stdTotalInfections);
TotalInfectionsUpper=min(Ntot-MeanTotalInfections,1.96*stdTotalInfections);

%collect data for summary cost benefit plot
meanTotalNumberOfInfections = MeanTotalInfections(end);
MeanTotalAbsentDays = mean(TotalAbsentDays);

MeanProportionOfInfectedIsolations=mean(sum(IsolatingInfections,2)./(sum(IsolatingInfections,2)+sum(IsolatingHealthy,2)),'omitnan');


end
