function ZF = initializeSCAPoints(SimParams)

    ZF = zeros(SimParams.gStatBeams,SimParams.nUsers);
    effChannel = zeros(SimParams.nUsers,SimParams.gStatBeams);
    for iGroup = 1:SimParams.nGroups
        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
            cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
            effChannel(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
            ZF(:,cUser) = effChannel(cUser,:)' / norm(effChannel(cUser,:));
        end
    end
    
    ZF = ZF * (sqrt(SimParams.txPower / sqrt(real(trace(ZF * ZF')))));
        
end
