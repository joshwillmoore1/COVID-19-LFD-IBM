
clc
clear 

%testing only starts at the sign of a symptomatic in the class

Ntot = 30; %total pop. size
SubgroupSize = 6; %number of individuals in each group
Irb = 2; %Background prevalence
R = 1.7; %R number
Cf = 1; %far R ratio
Cc = 1; %close R ratio
PerFalseNeg = 30; %percentage of false negatives from tests
CompIso = 81; %compliance with isolation for any individual
MixInd = 1;%see mixing days in genral function
TestInd = 2; %this now goes into the changeTestingSchedule function
PrevBool = false;
RecoveryBool = true; %bool  for infected recovery
PerSymptomatic = 40; %percentage of symptomatic
RemoveSympt = true; %Boolean to remove those that show symptoms
PerFalsePos = 0.3; %percentage of false positives
InitWithTesting = true; %start simulations with testing days active - false means to wait for sympto

Para = [Ntot,SubgroupSize,Irb,R,Cf,Cc,PerFalseNeg,CompIso,MixInd,TestInd...
    ,PrevBool,RecoveryBool,PerSymptomatic,RemoveSympt,PerFalsePos,InitWithTesting,1];

tic
[MeanNumOfHealthy,MeanTotalInfections,MeanNumberOfSympt...
    ,MeanNumberOfAsympt,initTestingDay,MeanNumOfRecoveries,...
    MeanNumberOfCasesDetected,MeanProportionOfInfectedIsolations] = General_Delay_LFT(Para);


%% plot the results
close all
subplot(1,4,1)
bar(1:28,MeanNumberOfAsympt+MeanNumberOfSympt,'FaceColor',[1 0 0])
alpha(0.5)
ylabel("Mean number of active infections")
xlabel("Number of days")
ylim([0,Ntot])

subplot(1,4,2)
plot(1:28,MeanTotalInfections,'r')
ylabel("Mean number of total infections")
xlabel("Number of days")
ylim([0,Ntot])

subplot(1,4,3)
bar(1:28,MeanNumberOfCasesDetected,'FaceColor',[0 0 1])
alpha(0.5)
ylabel("Mean number of isolating students")
xlabel("Number of days")
ylim([0,Ntot])

subplot(1,4,4)
bar(1:28,MeanNumOfRecoveries,'FaceColor',[0 .7 0])
alpha(0.5)
ylabel("Mean number of recoveries")
xlabel("Number of days")
ylim([0,Ntot])
