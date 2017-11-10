
function SimParams = digitalBeamformerDesign(SimParams)

sBeamM = SimParams.sBeamM;
precType = strsplit(SimParams.innerPrecoder,'_');

switch precType{1}
    
    case 'ZF'
        
        %         for iRun = 1:SimParams.nMontRuns
        %             B = zeros(SimParams.nTransmit,size(SimParams.groupInfo(1).statBeams,2),SimParams.nGroups);
        %             H = zeros(SimParams.nTransmit,SimParams.usersPerGroup(1),SimParams.nGroups);
        %             for iGroup = 1:SimParams.nGroups
        %                 for iUser = 1:SimParams.groupInfo(iGroup).nUsers
        %                     H(:,iUser,iGroup) = getChannel(SimParams,iGroup,iUser);
        %                 end
        %                 B(:,:,iGroup) = SimParams.groupInfo(iGroup).statBeams;
        %             end
        %
        %             GG = cell(SimParams.nGroups,1);
        %             for iGroup = 1:SimParams.nGroups
        %                 for jGroup = 1:SimParams.nGroups
        %                     GG{iGroup,1} = [GG{iGroup,1}, H(:,:,iGroup).' * B(:,:,jGroup)];
        %                 end
        %             end
        %
        
        %end
        
%         for iRun = 1:SimParams.nMontRuns
%             B = zeros(SimParams.nTransmit,size(SimParams.groupInfo(1).statBeams,2),SimParams.nGroups);
%             H = zeros(SimParams.nTransmit,SimParams.usersPerGroup(1),SimParams.nGroups);
%             for iGroup = 1:SimParams.nGroups
%                 for iUser = 1:SimParams.groupInfo(iGroup).nUsers
%                     H(:,iUser,iGroup) = getChannel(SimParams,iGroup,iUser);                    
%                 end
%                 B(:,:,iGroup) = SimParams.groupInfo(iGroup).statBeams;
%             end
%             
%             GG = cell(SimParams.nGroups,1);
%             for iGroup = 1:SimParams.nGroups
%                 for jGroup = 1:SimParams.nGroups
%                     GG{iGroup,1} = [GG{iGroup,1}, H(:,:,iGroup).' * B(:,:,jGroup)];                    
%                 end
%             end
%             
            
       % end
                 
        for iRun = 1:SimParams.nMontRuns
            zeroMatrix = zeros(SimParams.gStatBeams * SimParams.nGroups,SimParams.nUsers);
            H = zeros(SimParams.nUsers,SimParams.nTransmit);
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getChannel(SimParams,iGroup,iUser); %user specific channel
                    H(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser);
                    %zeroMatrix(:,cUser) = SimParams.groupInfo(iGroup).activeAntennas;
                end
            end
            
            effChannel = H * sBeamM; % * zeroMatrix;
            ZF = effChannel' * pinv(effChannel * effChannel');
            
            cvx_quiet('true');
            cvx_expert('true');
            
            cvx_begin
            
            variable absGain(SimParams.nUsers,1)
            variables userRate(SimParams.nUsers,1) userGamma(SimParams.nUsers,1)
            expressions xSignal totalIF
            
            maximize(sum(userRate))
            
            subject to
            
            for iGroup = 1:SimParams.nGroups
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    
                    cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                    
                    totalIF = SimParams.N0;
                    userRate(cUserIndex,1) <= log(1 + userGamma(cUserIndex,1));
                    
                    xSignal = abs(SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * ZF(:,cUserIndex))^2 * absGain(cUserIndex,1);
                    
                    %xSignal = abs(SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * zeroMatrix * ZF(:,cUserIndex))^2 * absGain(cUserIndex,1);
                    userGamma(cUserIndex,1) == xSignal / totalIF;
                end
                
            end
            
            %Total power constraint
            txPower = 0;
            for iUser = 1:SimParams.nUsers
                
                tempBeam = norm(sBeamM * ZF(:,iUser))^2 * absGain(iUser,1);
                
                %tempBeam = norm(sBeamM * zeroMatrix * ZF(:,iUser))^2 * absGain(iUser,1);
                txPower = txPower + tempBeam;
            end
            
            txPower <= SimParams.txPower;
            absGain >= 0;
            
            cvx_end
            
            if contains(cvx_status,'Solved')
                SimParams.sBeamM = sBeamM;
                ZF = ZF * sqrt(diag(absGain));
                SimParams.groupSumRate.srate(iRun,1) = cvx_optval;
                SimParams.groupSumRate.tsca(iRun,1) = 1;
            end
            
        end
        
    case 'CVX'
        
        % plotBeamPerformance(SimParams, sBeamM)
        
        cvx_quiet('true');
        cvx_expert('true');
        
        nSCAIterations = 10;
        SimParams.tStatBeams = size(sBeamM,2);
        
        for iRun = 1:SimParams.nMontRuns
            
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams = getChannel(SimParams,iGroup,iUser,'Reset');
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getChannel(SimParams,iGroup,iUser);
                end
            end
            
            SimParams.tStatBeams = size(sBeamM,2);
            cvxInnerMP = complex(randn(SimParams.tStatBeams,SimParams.nUsers),...
                randn(SimParams.tStatBeams,SimParams.nUsers)) / sqrt(2);
            userBetaP = rand(SimParams.nUsers,1) * 10;
            %             SimParams = evaluateUserRatesWithPerGroupPrecoders(SimParams,cvxInnerMP);
            %             userBetaP = SimParams.totUserBeta;
            
            
            
            for iSca = 1:nSCAIterations
                
                cvx_begin
                
                variable cvxInnerM(SimParams.tStatBeams,SimParams.nUsers) complex
                variables userRate(SimParams.nUsers,1) userGamma(SimParams.nUsers,1) userBeta(SimParams.nUsers,1)
                expressions xSignal totalIF
                
                maximize(sum(userRate))
                
                subject to
                
                for iGroup = 1:SimParams.nGroups
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        totalIF = sqrt(SimParams.N0);
                        
                        %Interference from my own group: Intragroup
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if jUser ~= iUser
                                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * cvxInnerM(:,xUserIndex)];
                            end
                        end
                        
                        %Interference from my neighbouring group: Intergroup
                        for jGroup = 1:SimParams.nGroups
                            if jGroup ~= iGroup
                                for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                    xUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                    totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * cvxInnerM(:,xUserIndex)];
                                end
                            end
                        end
                        
                        totalIF * totalIF' <= userBeta(cUserIndex,1);
                        userRate(cUserIndex,1) <= log(1 + userGamma(cUserIndex,1));
                        
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM;
                        xSignal = abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / userBetaP(cUserIndex,1) ...
                            + 2 * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) * (cvxInnerM(:,cUserIndex) - cvxInnerMP(:,cUserIndex)) / userBetaP(cUserIndex,1) ...
                            - (abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / (userBetaP(cUserIndex,1).^2)) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                        userGamma(cUserIndex,1) == xSignal;
                    end
                    
                    if SimParams.limitToGroupBeamsOnly
                        tempVector = cvxInnerM(~SimParams.groupInfo(iGroup).activeBeams,SimParams.groupInfo(iGroup).gUserIndices);
                        norm(tempVector(:)) <= 0;
                    end
                    
                end
                
                %Total power constraint
                txPower = 0;
                for iUser = 1:SimParams.nUsers
                    tempBeam = sBeamM * cvxInnerM(:,iUser);
                    txPower = txPower + tempBeam' * tempBeam;
                end
                
                txPower <= SimParams.txPower;
                userGamma >= 0;
                userBeta >= 0;
                
                cvx_end
                
                if contains(cvx_status,'Solved')
                    SimParams.sBeamM = sBeamM;
                    cvxInnerMP = cvxInnerM;
                    userBetaP = userBeta;
                    SimParams.groupSumRate.srate(iRun,iSca) = cvx_optval;
                    SimParams.groupSumRate.tsca(iRun,1) = iSca;
                    
                    for iGroup = 1:SimParams.nGroups
                        SimParams.groupInfo(iGroup).Perf(iRun,:) = userRate(SimParams.groupInfo(iGroup).gUserIndices);
                        SimParams.groupInfo(iGroup).InnerBF(:,:,iRun) = cvxInnerMP(:,SimParams.groupInfo(iGroup).gUserIndices);
                    end
                    
                    if ((std(SimParams.groupSumRate.srate(iRun,max(iSca - 4,1):iSca)) <= 1e-1) && (iSca > 10))
                        SimParams = getIntraGroupIF(SimParams,sBeamM,cvxInnerMP);
                        SimParams = getInterGroupIF(SimParams,sBeamM,cvxInnerMP);
                        fprintf('Completed with SR - %f \n',cvx_optval);
                        break;
                    else
                        fprintf('SR progress for SCA iteration - [%3d] %f \n',iSca,cvx_optval);
                    end
                else
                    fprintf('Resetting (failure !)');
                    userBetaP = rand(SimParams.nUsers,1) * sqrt(0.01) + userBetaP;
                    cvxInnerMP = cvxInnerMP + sqrt(0.01) * complex(randn(SimParams.tStatBeams,SimParams.nUsers),randn(SimParams.tStatBeams,SimParams.nUsers));
                end
                
            end
            
            
        end
        
    case 'CVXG'

        
       
        cvx_quiet('true');
        cvx_expert('true');
        cvx_solver('sdpt3');
        
        nSCAIterations = 10;

        SimParams.tStatBeams = size(sBeamM,2);
        
        for iRun = 1:SimParams.nMontRuns
            
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams = getChannel(SimParams,iGroup,iUser,'Reset');
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getChannel(SimParams,iGroup,iUser);
                end
            end
            
            cvxInnerMP = initializeSCAPoints(SimParams);
            SimParams = evaluateUserRatesWithPerGroupPrecoders(SimParams,cvxInnerMP);
            userBetaP = SimParams.totUserBeta;
            
            for iSca = 1:nSCAIterations
                
                cvx_begin
                
                variable cvxInnerM(SimParams.gStatBeams,SimParams.nUsers) complex
                variables userRate(SimParams.nUsers,1) userGamma(SimParams.nUsers,1) userBeta(SimParams.nUsers,1) groupBeta(SimParams.nUsers,SimParams.nGroups)
                expressions xSignal totalIF
                
                maximize(sum(userRate))
                
                subject to
                
                for iGroup = 1:SimParams.nGroups
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        totalIF = sqrt(SimParams.N0);
                        
                        %Interference from my own group: Intragroup
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if jUser ~= iUser
                                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,xUserIndex)];
                            end
                        end
                        %totalIF == real(totalIF);
                        
%                         Interference from my neighbouring group: Intergroup
%                         for jGroup = 1:SimParams.nGroups
%                             if jGroup ~= iGroup
%                                 xIF = [0];
%                                 for jUser = 1:SimParams.groupInfo(jGroup).nUsers
%                                     xUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
%                                     xIF = [xIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(jGroup).statBeams * cvxInnerM(:,xUserIndex)];
%                                 end
%                                 groupBeta(cUserIndex,jGroup) >= xIF * xIF';
%                             end
%                         end
%                         
%                         if strcmpi(precType{2},'Opt')
%                             totalIF * totalIF' + sum(groupBeta(cUserIndex,:)) <= userBeta(cUserIndex,1);
%                         else
%                             for jGroup = 1:SimParams.nGroups
%                                 if jGroup ~= iGroup
%                                     groupBeta(cUserIndex,jGroup) <= str2double(precType{2});
%                                 end
%                                 totalIF * totalIF' + str2double(precType{2}) * (SimParams.nGroups - 1) <= userBeta(cUserIndex,1);
%                             end
%                         end
                         
                        if strcmpi(precType{2},'Opt')
                            totalIF * totalIF' <= userBeta(cUserIndex,1);
                        else
                            totalIF * totalIF' + str2double(precType{2}) * (SimParams.nGroups - 1) <= userBeta(cUserIndex,1);
                        end
                            
                        
                        userRate(cUserIndex,1) <= log(1 + userGamma(cUserIndex,1));
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        
                        xSignal = abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / userBetaP(cUserIndex,1) ...
                            + 2 * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) * (cvxInnerM(:,cUserIndex) - cvxInnerMP(:,cUserIndex)) / userBetaP(cUserIndex,1) ...
                            - (abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / (userBetaP(cUserIndex,1).^2)) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                        userGamma(cUserIndex,1) == xSignal;
                        
                    end
                    
                end
                
                %Total power constraint
                txPower = 0;
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        tempBeam = SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,cUserIndex);
                        txPower = txPower + tempBeam' * tempBeam;
                    end
                end
                
                txPower <= SimParams.txPower;
                
                cvx_end
                
                if contains(cvx_status,'Solved')
                    SimParams.sBeamM = sBeamM;
                    cvxInnerMP = cvxInnerM;
                    userBetaP = userBeta;
                    SimParams.groupSumRate.srate(iRun,iSca) = cvx_optval;
                    SimParams.groupSumRate.tsca(iRun,1) = iSca;
                    
                    for iGroup = 1:SimParams.nGroups
                        SimParams.groupInfo(iGroup).Perf(iRun,:) = userRate(SimParams.groupInfo(iGroup).gUserIndices);
                        SimParams.groupInfo(iGroup).InnerBF(:,:,iRun) = cvxInnerMP(:,SimParams.groupInfo(iGroup).gUserIndices);
                    end
                    
                    if ((std(SimParams.groupSumRate.srate(iRun,max(iSca - 4,1):iSca)) <= 1e-2) && (iSca > 10))
                        fprintf('Completed with SR - %f \n',cvx_optval);
                        break;
                    else
                        fprintf('SR progress for SCA iteration - [%3d] %f \n',iSca,cvx_optval);
                    end
                else
                    fprintf('Resetting (failure !)');
                    userBetaP = rand(SimParams.nUsers,1) * sqrt(0.01) + userBetaP;
                    cvxInnerMP = complex(randn(SimParams.gStatBeams,SimParams.nUsers),randn(SimParams.gStatBeams,SimParams.nUsers));
                end
                
            end
            
            SimParams = evaluateUserRatesWithPerGroupPrecoders(SimParams, cvxInnerMP);
            SimParams.groupSumRate.srate(iRun,SimParams.groupSumRate.tsca(iRun,1)) = sum(SimParams.totUserRateE);
            
        end
        
        
end


end