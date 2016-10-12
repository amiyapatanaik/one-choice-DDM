% sample.mat has 7 days PVT data for two subjects. Each subject had 2 PVTs each day.
% RT data from each day is combined and DDM parameters are estimated from them
% As all five parameters of the DDM are not uniquely identifiable, drift and drift 
% variability parameters are normalized by the boundary parameter.
load('sample.mat');
N = size(RTs);
thetas = cell(N);
fits = zeros(N);
for i=1:N(1),
    for j=1:N(2),
        RT1 = RTs{i,j}{1};
        %RTs < 150 ms are fast guesses
        RT1(RT1<150) = [];

        RT2 = RTs{i,j}{2};
        RT2(RT2<150) = [];
        [param, fits(i,j)] = mixedEstimator([RT1;RT2]/1000,10);
        %Normalized drift and drift variability
        param(3) = param(3)/param(1);
        param(5) = param(5)/param(1);
        %boundary can be deleted 
        param(1) = [];
        thetas{i,j} = param;
    end
end

