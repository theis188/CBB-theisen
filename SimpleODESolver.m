function [uf,conc,bif] = SimpleODESolver(NoSteps,InitialState,U0,Uf, ...
    Kvec,DVDX,DVDU,S,JACOBIAN)
%SimpleODESolver: Fixed Step No Correction Solver
x = InitialState;
U = U0;
du = (Uf-U0)./NoSteps;
%%%
%Phosphok = 0;
%PERTUP = 0;
%Phosphok = max(find(du)==7);
%PERTUP = max(du>0);
%%%
conc = NaN(1+NoSteps,length(x));
conc(1,:) = x';
XJac = S*DVDX(x,Kvec,1,U);

for i = 2:1+NoSteps
    
    dx = -(XJac)\S*DVDU(x,Kvec,1,U)*du;
    x = x + dx;
    U = U + du;
    
    if min(x) < 0
        uf = U;
        bif = 1;
        return
    else
        XJac = S*DVDX(x,Kvec,1,U);  % Calculate new Jacobian
        if max(real(eig(XJac))) > -1e-6
            uf = U;
            bif = 1;
            %%%
            %if Phosphok == 1 && PERTUP == 1;
            %i/NoSteps
            %    if i/NoSteps > 0.3 && i/NoSteps <0.5;
            %        Kvec
            %        pause
            %    end
            %end
            %%%
            return
        end
        conc(i,:) = x';
    end
    
end

uf = U;
bif = 0;