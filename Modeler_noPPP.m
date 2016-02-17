clear
clc

load('CBB_Model.mat')

%% AddRxns

% RxnsToAdd = {...
%     'LipSyn',{'AcCoA'},[-1]},...
% };
% 
% for i = 1:size(Rxns
% EnzName{end + 1} = RxnsToAdd
% 

%% Fix S

S(find(strcmp('E4P',MetabName)),find(strcmp('Rpe',EnzName))) = 0;
S(find(strcmp('X5P',MetabName)),find(strcmp('Rpe',EnzName))) = -1;
S(find(strcmp('X5P',MetabName)),find(strcmp('Tkt2',EnzName))) = 1;
S(find(strcmp('R5P',MetabName)),find(strcmp('Tkt2',EnzName))) = 1;

%% DeleteMetabs
DeleteMetabs = {'2PG';'PEP';'PYR';'AcCoA';'6PGL';'6PGT';'SBA';'AMP'};

for i = 1:size(DeleteMetabs,1);
    CurrMetabIndex = find(strcmp(DeleteMetabs{i,1},MetabName));
    S(CurrMetabIndex,:) = [];
    MetabName(CurrMetabIndex,:) = [];
end

%% DelRxns
DeleteRxns = {'G6PDH';'GLNase';'6PGLDH'};

for i = 1:size(DeleteRxns,1);
    CurrRxnIndex = find(strcmp(DeleteRxns{i,1},EnzName));
    S(:,CurrRxnIndex) = [];
    EnzName(CurrRxnIndex,:) = [];
    Reversibilities(CurrRxnIndex,:) = [];
end

%% Add Rxns
AddRxns = {
{'SSr'},{'ADP','ADPG'},[-1 1],0;
{'G3Pr'},{'G3P'},[1],0;
};


for j = 1:size(AddRxns,1)
    S = [S zeros(size(S,1),1)]
    for i = 1:size(AddRxns{j,2},2)
        MetabIndex = find(strcmp(AddRxns{j,2}{i},MetabName));
        S(MetabIndex,end) = AddRxns{j,3}(1,i);
    end
    Reversibilities(end+1,1) = AddRxns{j,4};
    EnzName{end+1,1} = AddRxns{j,1}{1};
end


% Add Metabs = 
AddMetabs = {
{'Pi'},[0;0;1;0;0;1;0;0;1;0;0;0;0;1;0;0;2;0;0;-2;0;-1;]';
}

for i = 1:size(AddMetabs,1);
    S(end + 1,:) = AddMetabs{i,2}
    MetabName{end + 1,1} = AddMetabs{i,1}{1}
end


%% SpecifyRatios

Seq = S;

% % specifyRatios = [
% % {'SS','G3P_Out'},-2 ,1;
% % ];
% % 
% % for i=1:size(specifyRatios,1);
% %     rxn1Index = find(strcmp(specifyRatios{i,1},EnzName));
% %     rxn2Index = find(strcmp(specifyRatios{i,2},EnzName));
% %     newRow = zeros(1,size(S,2));
% %     newRow(1,rxn1Index) = specifyRatios{i,3};
% %     newRow(1,rxn2Index) = specifyRatios{i,4};
% %     Seq = [Seq; newRow];
% % end

%% Specify Reversibilities

Reversibilities(find(strcmp('G3P_Out',EnzName))) = 1;
Reversibilities(find(strcmp('SS',EnzName))) = 1;


%% Fixing everything

beq = zeros(size(Seq,1),1);
lb = Reversibilities*-1000;
ub = 1000*ones(size(lb));

specifyBounds = [
{'RuBiSCO'},10;
{'SS'},3;
{'SSr'},2;
{'G3Pr'},3;
];

for i = 1:size(specifyBounds, 1)
    lb(find(strcmp(specifyBounds{i,1},EnzName))) = specifyBounds{i,2}*0.99;
    ub(find(strcmp(specifyBounds{i,1},EnzName))) = specifyBounds{i,2}*1.01;
end

f = zeros(length(EnzName),1);
f(find(strcmp('SS',EnzName))) = 1;

Vref = linprog(f,[],[],Seq,beq,lb,ub);

Net.S = S;
Net.MetabName = MetabName;
Net.EnzName = EnzName;
Net.Sreg = zeros(size(S));
Net.Vref  = Vref;
Net.Reversibilities = Reversibilities;
%%Fix Sreg

Net.Sreg(find(strcmp('NADPH',MetabName)),find(strcmp('G6PDH',EnzName))) = 1;
Net.Sreg(find(strcmp('FBP',MetabName)),find(strcmp('PGM',EnzName))) = 1;
Net.Sreg(find(strcmp('RuBP',MetabName)),find(strcmp('PGM',EnzName))) = 1;
%Net.Sreg(find(strcmp('NADP',MetabName)),find(strcmp('TPI',EnzName))) = 1;
%Net.Sreg(find(strcmp('NADPH',MetabName)),find(strcmp('TPI',EnzName))) = 2;


MainPerturb(Net,'noPPP_WPi_300')