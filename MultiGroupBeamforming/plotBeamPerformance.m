function plotBeamPerformance(SimParams,sBeamM)

nResols = 1024;
cTheta = linspace(-180,179,nResols);
xTheta = wrapTo180(cTheta + cTheta);

plotGains = zeros(nResols,2,size(sBeamM,2));
elementLocs = linspace(0,(SimParams.nTransmit - 1),SimParams.nTransmit) * 0.5;

for iTheta = 1:nResols
    for iBeam = 1:size(sBeamM,2)
        plotGains(iTheta,2,iBeam) = cbfweights(elementLocs,xTheta(1,iTheta))' * sBeamM(:,iBeam);
        plotGains(iTheta,1,iBeam) = xTheta(1,iTheta);
    end
end

for iBeam = 1:size(sBeamM,2)
    plot(plotGains(:,1,iBeam),abs(plotGains(:,2,iBeam)));hold all;
end

end

