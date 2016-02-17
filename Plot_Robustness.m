clear

ResultName = {'noPPP_Revs_50_WPi'}
%{...
%    'ModelMattMC2Model_Lo_glc_flux Results';...
%    'ModelMattMC2Model_VLo_glc_flux Results';...
%    'ModelMattMC2Model_VVLo_glc_flux Results';}

load(strcat(ResultName{1,1},'.mat'));

TestedUs = [PertDown:(1-PertDown)/StepsDown:1 1+(PertUp-1)/StepsUp:(PertUp-1)/StepsUp:PertUp];
[DVDU2Us, IX] = sort(TestedUs([find(TestedUs==1):end find(TestedUs==1):-1:1 ]));

figure
hold on

% PlotOrder = [2 1 3:size(EnzName,1)-2];

for Count = [1:size(EnzName,1)];
    %PlotOrder = [2 3 4 5 6 7 9 10]; %PlotOrder(Count)
    Enzyme = Count
%Find fraction of models that are robust
for j = 1:EnsembleSize;
IsRobust(j,:) = (isnan(ModelResults{j,1}(1,IX,Enzyme))==0);
end
Robust_Frac(Enzyme,:) = mean(IsRobust);

Edge_Length = ceil(sqrt(size(EnzName,1)-2))
subplot(5,5,Count)
plot(log10(DVDU2Us), smooth(Robust_Frac(Enzyme,:)),'b','LineWidth',2)
%if Count == 7
    set(gca,'LineWidth',2, 'fontsize', 11, 'FontName', 'Cambria')
    set(gca,'xtick',[-1 0 1], 'xticklabel',{'0.1' '1' '10'}, 'fontsize', 11, 'fontname', 'Cambria')
    set(gca,'ytick',[0 0.5 1], 'yticklabel',{'0' '0.5' '1'}, 'fontsize', 11, 'fontname', 'Cambria')
%else
%    set(gca,'LineWidth',2, 'fontsize', 11, 'FontName', 'Cambria')
%    set(gca,'xtick',[-1 0 1],'xticklabel',{' ', ' ', ' '}, 'fontsize', 11, 'FontName', 'Cambria')
%    set(gca,'ytick',[0 0.5 1],'yticklabel',{' ', ' ', ' '}, 'fontsize', 11, 'FontName', 'Cambria') 
%end

annotation('textarrow',[0.09 0.5],[0.5 0.5],'string','Fraction of Stable Models', ...
      'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,...
      'fontname', 'Cambria', 'fontsize', 16);

annotation('textarrow',[0.5 0.5],[0.06 0.5],'string','Fold Change in Enzyme Level', ...
      'HeadStyle','none','LineStyle', 'none', 'TextRotation',0,...
      'fontname', 'Cambria', 'fontsize', 16)
 %for j=1:5;

%plot(log10(DVDU2Us), Average_Fluxes(2,:),'g')
%plot(log10(DVDU2Us), Average_Fluxes(3,:),'r')

%legend('Glc Flux = 0.2','Glc Flux = 0.01','Glc Flux = 1e-5','Location','NorthWest')
%xlabel(strcat('Log_1_0 of ',EnzName{Enzyme,1},' Activity Fold Change'))
xlabel(EnzName{Enzyme,1}(1:min(6,numel(EnzName{Enzyme,1}))), 'fontname', 'Cambria')
ylim([0 1.2])
xlabh = get(gca,'XLabel');
%set(xlabh,'Position',[0 -0.1+(Count==7)*0.4 0])
set(gcf,'color','white')
end

%{

for i=1:size(EnzName,1)
    
    subplot(S,S,i)
        hold on
        plot(log10(DVDU2Us), ModelResults{j,1}(Metabolite,IX,i))
    end
    xlabel(EnzName{i});
    ylabel(Net.MetabName{Metabolite+2});
    
end
%}

