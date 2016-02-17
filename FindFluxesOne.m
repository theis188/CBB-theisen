%FLUXES(GoodResults(:,PertStep,i), ...
%        EnsembleKvec(:,i), 1, U + (Uf1 - U)*(PertStep/200));
    
% FLUXES(X, K, 1, U)
clear all
load('30_pct_oPPP_rev_300_F');
%EnsembleSize = 100;
%StepsUp = 800;
BifurcPoints = zeros(size(EnzName,1),EnsembleSize,2); %(:,:,1) is up --- (:,:,2) is down
IsoButFluxes = zeros(1,EnsembleSize,1); %(:,:,1) is up --- (:,:,2) is down
U = ones(size(EnzName,1),1);
UUp = ones(size(EnzName,1),1);
UDown = ones(size(EnzName,1),1);

for a = 1:EnsembleSize; %finding point of bifurcation (or 10x up/down)
    for j = UniqueEnzymes;%1:size(EnzName,1);
        BifurcPoints(j,a,1) = max(find(isnan(ModelResults{a,1}(1,1:1 + StepsUp,j))==0));
        %BifurcPoints(j,a,2) = max(find(isnan(ModelResults{a,1}(1,202:402,j))==0));
        %-------------------------------------Up
        UUp = ones(size(EnzName,1),1);
        UUp(j) = PertUp; %FLUXES(X,Kvec,1,U)
        CurrFluxes = FLUXES(... 
            ModelResults{a,1}(:,BifurcPoints(j,a,1),j), ...
            EnsembleKvec(:,a), ...
            1, ...
            U + (UUp - U)*((BifurcPoints(j,a,1) - 1)/StepsUp));
        %CurrFluxes(39)
        AllFluxes(:,a,j) = CurrFluxes;%(find(strcmp('EXCH_lac-l(e)',EnzName)));
        
        %---------------------------------------Dwon
%         UDown = ones(size(EnzName,1),1);
%         UDown(j) = 0.1;
%         CurrFluxes = FLUXES(...
%             ModelResults{a,1}(:,BifurcPoints(j,a,2)+201,j), ...
%             EnsembleKvec(:,a), ...
%             1, ...
%             U + (UDown - U)*((BifurcPoints(j,a,2) - 1)/200));
%         IsoButFluxes(j,a,2) = CurrFluxes(39);
    end
end

RefFlux = Net.Vref;
stdFluxes = std(AllFluxes,1,2)
AvgFluxes = mean(AllFluxes,2);%/RefFlux(find(strcmp('EXCH_(e)',EnzName)));
%AvgFluxes = AvgFluxes(:,[2 1]);
%StabFrac = [mean(BifurcPoints(:,:,2)==201,2) ...
StabFrac = [mean(BifurcPoints(:,:,1)==StepsUp + 1,2)];

%Report = [StabFrac(1) AvgFluxes]

RefFlux(find(strcmp('RuBiSCO',EnzName)))
AvgFluxes(find(strcmp('RuBiSCO',EnzName)),1,find(strcmp('SBPase',EnzName)))
AvgFluxes(find(strcmp('RuBiSCO',EnzName)),1,find(strcmp('FBPase',EnzName)))
AvgFluxes(find(strcmp('RuBiSCO',EnzName)),1,find(strcmp('RuBiSCO',EnzName)))

A = permute(AvgFluxes(find(strcmp('RuBiSCO',EnzName)),1,[1:13]),[3 2 1]);
B = permute(stdFluxes(find(strcmp('RuBiSCO',EnzName)),1,[1:13]),[3 2 1]);
C = permute(AvgFluxes(find(strcmp('G6PDH',EnzName)),1,[1:13]),[3 2 1]);
D = permute(stdFluxes(find(strcmp('G6PDH',EnzName)),1,[1:13]),[3 2 1]);

% AvgFluxes(find(strcmp('EXCH_leu-l(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_asp(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_ala-l(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_ibut(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_etoh(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_mal(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_lac-l(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_for(e)',EnzName)))
% AvgFluxes(find(strcmp('EXCH_ac(e)',EnzName)))

