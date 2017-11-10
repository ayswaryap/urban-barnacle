
function SimParams = outerBeamformerDesign(SimParams)

% * * Depending upon the choice of outer BF we design the outer BF * *


switch SimParams.statBeamType
    
    case {'Eigen'}
        
        for iGroup = 1:SimParams.nGroups
            covMatrix = zeros(SimParams.nTransmit);
            for iRealization = 1:SimParams.nRealizations
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams = getChannel(SimParams,iGroup,iUser,'Reset');
                    cH = getChannel(SimParams,iGroup,iUser);
                    covMatrix = cH' * cH + covMatrix;
                end
            end
            SimParams.groupInfo(iGroup).covMatrix = covMatrix / SimParams.nRealizations;
        end
        
        for iGroup = 1:SimParams.nGroups
            [U,~,~] = svd(SimParams.groupInfo(iGroup).covMatrix);
            SimParams.groupInfo(iGroup).statBeams = U(:,1:SimParams.gStatBeams);
            %plot = plotBeamPerformance(SimParams, SimParams.groupInfo(iGroup).statBeams);
        end
        
    case 'GDFT'
        
        U = zeros(SimParams.nTransmit,SimParams.gStatBeams);
        for iGroup = 1:SimParams.nGroups
            xDivision = length(SimParams.groupInfo(iGroup).baseTheta) * SimParams.nGroups;
            beamDirection = linspace(-(180 / xDivision),(180 / xDivision),(xDivision * 2) + 1);
            beamDirection = repmat(beamDirection(2:2:end),length(SimParams.groupInfo(iGroup).baseTheta),1);
            beamDirection = beamDirection + repmat(SimParams.groupInfo(iGroup).baseTheta',1,size(beamDirection,2));
            
            beamDirection = reshape(beamDirection,1,[]);
            elementLocs = linspace(0,(SimParams.nTransmit - 1),SimParams.nTransmit) * 0.5;
            for iBeam = 1:SimParams.gStatBeams
                U(:,iBeam) = steervec(elementLocs,wrapTo180(beamDirection(1,iBeam)));
            end
            SimParams.groupInfo(iGroup).statBeams = U;
            SimParams.groupInfo(iGroup).beamDirections = beamDirection;
        end
        
    case 'Eye'
        
        for iGroup = 1:SimParams.nGroups
            SimParams.limitToGroupBeamsOnly = 0;
            SimParams.groupInfo(iGroup).activeAntennas = ones(size(SimParams.nTransmit));
            SimParams.groupInfo(iGroup).statBeams = eye(SimParams.nTransmit);
        end
        SimParams.gStatBeams = SimParams.nTransmit;
        SimParams.innerPrecoder = 'CVX_Opt';
        
    case 'BEAMS'
        oversampFactor = 2^0;groupAngularSpread = 60 / SimParams.nGroups;
        U = zeros(SimParams.nTransmit,SimParams.gStatBeams * oversampFactor);
        couplingMatrix = kron(eye(SimParams.gStatBeams),ones(oversampFactor,1));
        for iGroup = 1:SimParams.nGroups
            hBeam = mean(SimParams.groupInfo(iGroup).baseTheta);
            thetaBegin = hBeam - groupAngularSpread;
            thetaEnd = hBeam + groupAngularSpread;
            beamDirections = wrapTo180(linspace(thetaBegin,thetaEnd,oversampFactor * SimParams.gStatBeams));
            for iBeam = 1:length(beamDirections)
                U(:,iBeam) = conj(step(SimParams.stVector,SimParams.cFrequency,beamDirections(1,iBeam)));
            end
            SimParams.groupInfo(iGroup).statBeams = U * couplingMatrix;
            SimParams.groupInfo(iGroup).beamDirections = beamDirections;
            SimParams.groupInfo(iGroup).statBeams = SimParams.groupInfo(iGroup).statBeams / (norm(SimParams.groupInfo(iGroup).statBeams,'fro') * sqrt(numel(U)));
            plotBeamPerformance(SimParams, sum(SimParams.groupInfo(iGroup).statBeams,2));
        end
        
        SimParams.stVector.release;       
        
end

%Depending on the choice of outer BF we design the B accordingly

if any(strcmpi(SimParams.statBeamType,{'Eigen','GDFT','BEAMS'}))
    SimParams.sBeamM = [];
    activeBeamsPerGroup = zeros(1,SimParams.totStatBeams);
    activeBeamsPerGroup(1:SimParams.gStatBeams) = 1;
    for iGroup = 1:SimParams.nGroups
        SimParams.groupInfo(iGroup).activeBeams = circshift(activeBeamsPerGroup, SimParams.gStatBeams * (iGroup - 1), 2);
        SimParams.sBeamM = [SimParams.sBeamM, SimParams.groupInfo(iGroup).statBeams];
    end
else
    SimParams.sBeamM = eye(SimParams.nTransmit);
end

end