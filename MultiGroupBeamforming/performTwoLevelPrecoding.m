

function gParams = performTwoLevelPrecoding(inParams)

gParams = cell(length(inParams.beamsPerGroup),1);

cBeamIndex = 0;
for iBeamPerGroup = inParams.beamsPerGroup
    
    %SimParams.rng = inParams.rng;
    
    rng('shuffle');
    cBeamIndex = cBeamIndex + 1;
        
    SimParams.frequency = inParams.frequency;    
    SimParams.legendName = inParams.simulationLegend;
    SimParams.usersPerGroup = inParams.usersPerGroup;
    
    SimParams.nReceive = inParams.nReceive;
    SimParams.nTransmit = inParams.nTransmit;
    
    SimParams.cFrequency = 3e9;
    SimParams.sFrequency = 1e6;
    SimParams.numScatterers = 20;
    SimParams.userRingRadius = 30;
    
    SimParams.lambda = physconst('LightSpeed') / SimParams.cFrequency;
    SimParams.elemSpacing = SimParams.lambda * 0.5;

    SimParams.txArray = phased.ULA('NumElements',SimParams.nTransmit,'ElementSpacing',SimParams.elemSpacing,'ArrayAxis','y');
    SimParams.stVector = phased.SteeringVector('SensorArray',SimParams.txArray);

    if SimParams.nReceive == 1
        SimParams.rxArray = phased.ConformalArray();
    else
        SimParams.rxArray = phased.ULA('NumElements',SimParams.nReceive,'ElementSpacing',SimParams.elemSpacing,'ArrayAxis','y');
    end
    
    SimParams.fixedGroups = 16;
    SimParams.nRealizations = 100;
    
    SimParams.N0 = 1;
    SimParams.txPower = 10^(inParams.txSNR/10);
    SimParams.nUsers = sum(SimParams.usersPerGroup);
    
    SimParams.nMontRuns = inParams.nDrops;
    SimParams.statBeamType = inParams.statBeamType;
    SimParams.innerPrecoder = inParams.innerPrecoderType;
    SimParams.nGroups = length(SimParams.usersPerGroup);
    SimParams.beamsPerGroup = inParams.beamsPerGroup;
    
    SimParams.totStatBeams = iBeamPerGroup;
    SimParams.gStatBeams = floor(iBeamPerGroup / SimParams.nGroups);
    SimParams.groupUserIndices = cell(SimParams.nGroups,1);
    
    SimParams.limitToGroupBeamsOnly = inParams.limitToGroupBeamsOnly;
    
    SimParams.fixedGroupLocs = linspace(-60,60,SimParams.fixedGroups + 1);
    SimParams.fixedGroupLocs = SimParams.fixedGroupLocs(2:2:SimParams.fixedGroups);
    userClustersPerGroup = ceil(length(SimParams.fixedGroupLocs) / SimParams.nGroups);
    
    SimParams.angSpread = inParams.uAngularSpread;
    SimParams.antennasPerGroup = SimParams.nTransmit / SimParams.nGroups; %This is not required
    
    userIndex = 0;
    userIndices = randperm(SimParams.nUsers);
    for iGroup = 1:SimParams.nGroups
        possibleBeams = SimParams.fixedGroupLocs((iGroup - 1) * userClustersPerGroup + 1 : min(iGroup * userClustersPerGroup,length(SimParams.fixedGroupLocs)));
        SimParams.userLocs{iGroup,1} = possibleBeams(randi([1, length(possibleBeams)],1,SimParams.usersPerGroup(1,iGroup)));
        SimParams.groupUserIndices{iGroup,1} = userIndices((userIndex + 1) : sum(SimParams.usersPerGroup(1:iGroup)));
        SimParams.groupUserBaseAngles{iGroup,1} = possibleBeams;
        userIndex = sum(SimParams.usersPerGroup(1:iGroup));
    end
    
    SimParams.groupInfo = struct();
    for iGroup = 1:SimParams.nGroups
        SimParams.groupInfo(iGroup).nUsers = SimParams.usersPerGroup(1,iGroup);
        SimParams.groupInfo(iGroup).userLocs = SimParams.userLocs{iGroup,1};
        SimParams.groupInfo(iGroup).baseTheta = SimParams.groupUserBaseAngles{iGroup,1};
        SimParams.groupInfo(iGroup).gUserIndices = SimParams.groupUserIndices{iGroup,1};
        SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
            SimParams.groupInfo(iGroup).mimoChannel{iUser,1} = [];
            SimParams = getChannel(SimParams,iGroup,iUser,'Initialize');
        end
    end
    
    SimParams = outerBeamformerDesign(SimParams);
    SimParams = digitalBeamformerDesign(SimParams);
    gParams{cBeamIndex,1} = SimParams;
    
end

end

