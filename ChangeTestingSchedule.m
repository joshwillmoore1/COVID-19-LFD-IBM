function TestingSchedule=ChangeTestingSchedule(Flag,Day,TestingScenario)
% Initially there are no tests running. Once a symptomatic is found the
% Testing Schedule is updated to reflect a number of scenarios. In all
% scenarios a test will occur the next day after a symptomatic is
% found (or Monday if they are found on a Friday).

% Flag = 0 if there are no Symptomatics
% Flag = 1 if there are symptomatics
% Flag can be changed to whatever is needed

% Day is the day number on which the Symptomatic is found
if Day<0||Day>28
    error("EXPECTION: Day should be and integer between 1 and 28, inclusive.")
end
% TestingScenario flags the recurrence of a testing:
% TestingScenario = 0 -> no testing
% TestingScenario = 1 -> once a week
% TestingScenario = 2 -> twice a wekk
% TestingScenario = 5 -> daily
% TestingScenario = 6 -> (HE arrival) once a week with a prior test
% TestingScenario = 6 -> (HE arrival) once a week


if ~ismember(TestingScenario,[0 1 2 5 6 7])
    error("EXPECTION: TestingScenario should be 0, 1, 2, or 5")
end

if Flag==1
    % If you're checking that there is a symptomatic before this you can
    % delete this check
    
    if TestingScenario==0
        % No Testing
        %                 [S,M,T,W,T,F,S,S,M,T,W,T,F,S,S,M,T,W,T,F,S,S,M,T,W,T,F,S]
        TestingSchedule = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    elseif TestingScenario==1
        % Once a week, Monday
        TestingSchedule = [0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0];
    elseif TestingScenario==2
        % Twice a week, Monday - Friday
        TestingSchedule = [0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0];
    elseif TestingScenario==5
        % 5 days a week
        TestingSchedule = [0,1,1,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,1,1,0];
    elseif TestingScenario==6
        % Prearrival Test - HE - one pre-arriv Fri and Mon weekly
        TestingSchedule = [0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0];
    elseif TestingScenario==7
        % Prearrival Test - HE - one weekly on Mon
        TestingSchedule = [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0];
    end
    
    %Zero all testing days on Day and before
    TestingSchedule(1:Day)=0;
    
    %Make the next applicable day a testing day. Since we always test on a
    %Monday we do not technically need to define the first two cases,
    %because they are simply putting a Monday as a testing day, which is
    %already done. However, I am keeping it in just in case we want to vary
    %the weekly days, e.g. Wednesday and Friday.
    
    if Day<27 %If Day is 27 or greater then we are on the last Friday, so no more testing can be added.
        if mod(Day,7)==0 %Check if symptomatic is found on a Saturday.
            TestingSchedule(Day+2)=1;
        elseif mod(Day,7)==6 %Check if symptomatic is found on a Friday.
            TestingSchedule(Day+3)=1;
        else %Make next day a testing day.
            TestingSchedule(Day+1)=1;
        end
    end
else
    TestingSchedule=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
end

end

