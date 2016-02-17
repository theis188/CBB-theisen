%%   Prepare  Set of Reactions
%    load(strcat(ModelID,'.mat'));
    
  
    [ParamInfo, rVnet,Sreg,EnzName,S, V, KVEC, X, SubstrateConc, ParamRange, rV,U] = PrepareLumpedEMVariables(Net);

    %% Find conserved moeities and replace them with algebraic expressions
    Xini = ones(length(X),1);
    Uini = ones(length(U),1);
   
   %if isempty(findstr('TD',ModelID))
    NullMat = null(S', 'r');
    NullMat = NullMat>1e-5;
    ToRemove = [];
    
    for n=1:size(NullMat,2),
        %Replace every set of conserved moeties with an algebraic
        %expression
        CurrentVec = NullMat(:,n);
        OtherVec = NullMat;
        OtherVec(:,n) = [];
        tempRemov = find(CurrentVec~=0 & all(OtherVec==0,2));
        ToRemove = [ToRemove tempRemov(1)];
        XExpression = solve(CurrentVec'*X - CurrentVec'*Xini, X(tempRemov(1)))
        Net.MetabName(tempRemov);
        V = subs(V, X(tempRemov(1)), XExpression);

        %SYMK1S = subs(SYMK1S, X(tempRemov(1)), XExpression);
        %K1S = matlabFunction(SYMK1S, 'vars', [{X} {KVEC} {SubstrateConc} {SYMVNET} {U}]);
    end
    Net.MetabName(ToRemove, 1)
    %Net.MetabName(ToRemove, :) = [];
    S(ToRemove, :) = [];
    X(ToRemove) = [];
    Xini(ToRemove) = [];
    %end



    
       
    
    %% Prepare Solution related functions
    % Calculate Dx
    PREDX = S*V;
    %Calculate the jacobian matrix
    SYMJAC = jacobian(PREDX,X);
    save('SYMJACFile','SYMJAC');
    %Calculate the flux jacobians
    SYMDVDX = jacobian(V,X);
    SYMDVDU = jacobian(V,U);
    %Set up for calculting Vf after sampling the other parameters
    SYMVNET = sym('VNET', [length(V) 1]);
    SYMK1S = sym('K1S', [length(V) 1]);
    for m=1:length(V),
        SYMK1S(m) = solve(V(m)-SYMVNET(m), KVEC(ParamInfo(m,1)));
        disp(m)
    end
    % Convert everything into matlab functions for speed
    syms t;
    %if ~isempty(findstr('TD',ModelID));
    %DX = matlabFunction(PREDX, 'vars', [{X} {KVEC} {SubstrateConc} {U}]);% {G}]);
    %else
    DX = matlabFunction(PREDX, 'vars', [{X} {KVEC} {SubstrateConc} {U}]);
    JACOBIAN = matlabFunction(SYMJAC, 'vars', [{t} {X}, {KVEC} {SubstrateConc} {U}]);
    K1S = matlabFunction(SYMK1S, 'vars', [{X} {KVEC} {SubstrateConc} {SYMVNET} {U}]);
    FLUXES = matlabFunction(V, 'vars', [{X} {KVEC} {SubstrateConc} {U}]);
    DVDX = matlabFunction(SYMDVDX, 'vars', [{X}, {KVEC} {SubstrateConc} {U}]);
    DVDU = matlabFunction(SYMDVDU, 'vars', [{X}, {KVEC} {SubstrateConc} {U}]);
    %end
    B = null(S);
