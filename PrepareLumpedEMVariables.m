function [ParamInfo, rVnet,Sreg, EnzName,S, V, KVEC, X, SubstrateConc, ParameterRange, rV, U]=PrepareLumpedEMVariables(Net)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The code below sets up all the variables to be used later. There is no
%need to change this
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input the stoichiometric matrix along with required indices and fluxes

S=Net.S;
Sreg=Net.Sreg;     %%%%%%%% This script will not use regulation
EnzName=Net.EnzName;
%% Identify passive transport for making Net structure
Vout_index = [];
Vin_index = [];
for n=1:length(EnzName)
    if sum(S(:,n)>0) == 0,
        Vout_index = [Vout_index; n];
    end
    if sum(S(:,n)<0) == 0,
        Vin_index = [Vin_index; n];
    end
end
% Vcof_index = Net.Vcof_index;
% input the thermodynamic constraints
rVnet=Net.Vref;
Revs = Net.Reversibilities;

%% Create the system of ODEs
TotalParams = FindTotalParams(S, Sreg, Vin_index, Vout_index,Revs);
KVEC = sym('KVEC', [TotalParams 1]);
U = sym('U', [length(EnzName) 1]);
X = sym('X', [size(S,1) 1]);

GroupIndices = {}

%if ~isempty(findstr('Isots',ModelID));
%G = sym('G', [Net.Groups 1]);
%GroupIndices = {};
%GroupName = {};
%    for i = 1:Net.Groups;
%        GroupName{i,1} = Net.MetabName{i,1};
%        GroupIndices{i,1} = find(strncmpi(GroupName{i,1},Net.MetabName,numel(GroupName{i,1})));
%    end
%[V, ParamInfo, SubstrateConc, ParameterRange, rV] = SetUpRxnRates_Isotopic(S,Sreg, KVEC, X, Vin_index, Vout_index, U,rVnet,Revs,G,GroupIndices);
%else
[V, ParamInfo, SubstrateConc, ParameterRange, rV] = SetUpRxnRates(S,Sreg, KVEC, X, Vin_index, Vout_index, U,rVnet,Revs);
end

% for n=1:length(Vcof_index),
%     V = subs(V,X(Vcof_index(n)+1), 2-X(Vcof_index(n)));
% end
% S(Vcof_index+1, :) = [];
% X(Vcof_index+1, :) = [];


