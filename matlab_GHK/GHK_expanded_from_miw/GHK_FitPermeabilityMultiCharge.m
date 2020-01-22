
function [ PositiveIonPermeability, NegativeIonPermeability, Offset ] = GHK_FitPermeabilityMultiCharge( IV, ConcI, Conc0, Vm,z )
%Fit the permebaility to a curve assuming 2 ions one positive and
%one negative in equal concentration.
%
% NOTE THE Conc0 is the ground or reservoir. ConcI is the capillary
R_Const = 8.3144621;
F_Const = 9.64853399e4;
%Create a temperature const so can change temperature though out all
%equations easily and consistently
T_Const = 300;
e_Const = 1.6021766e-19;
Kb_Const = 1.38064852e-23;
epsilon_0_Const = 8.8541878176e-12;
N_A_Const = 6.022140857e23;    %Avogadros number;

R = 8.3144621;
F = 9.64853399e4;
T = 300;

if nargin<4
    %A voltage to remove Vm - need to fix in the center..
    Vm = -1 *  mean(IV(:,2)) * 1e-3;
else
    Vm = -1 * Vm * 1e-3;
end

%Can clean the IV curve to remove bad points.
IVc = IVClean(IV);

I = IVc(:,1) ;%* 1e-9;
V = IVc(:,2) * 1e-3 + Vm;

%Use a linear fit  - calculate x values using equation with P = 1
%z = [1 -1];
n = max(size(V));

GHK = [GHK_CurrentByMultiIon_S( z(1), V , 1, ConcI*(abs(z(2)/z(1))), Conc0*(abs(z(2)/z(1))) ) GHK_CurrentByMultiIon_S( z(2), V , 1, ConcI*(abs(z(1)/z(2))), Conc0*(abs(z(1)/z(2))) ) ones(n,1)];

Params0 = [1e-7 1e-7 0];
%Does the fitting with V and nA  - so need to correct later
P = fminsearch(@(Params0) Cost( I, GHK, Params0  ), Params0);
P(1)
P(2)
CostError = (Cost( I, GHK, P ));
m = mean(IVc(:,1));
sd = sqrt(sum((IVc(:,1) - m).^2 ));
%e = CostError/(sd)
%e = CostError/(max(IVc(:,1))- min(IVc(:,1)))
e = CostError;


%Do any necessary Scaling to get into SI.
PositiveIonPermeability = abs(P(1)) * 1e-9;
NegativeIonPermeability = abs(P(2)) * 1e-9;
Offset = P(3)*1e-9;

%Plot the result against the raw data  - in SI Units
[It , Ic] = GHK_TotalCurrent_Multi(z,V, [PositiveIonPermeability NegativeIonPermeability] , ConcI, Conc0);

%Convert  - into nA and mV for plotting.
It = It * 1e9;
Ic = Ic * 1e9;
offsetnA = Offset * 1e9;

%Plot it in nA and mV
%figure(20);
hold off;
plot((V-Vm)*1000,Ic + offsetnA,' - ');
hold all
plot((V-Vm)*1000,It + offsetnA);
hold on
xlabel({['Voltage mV'],['ConcI: ' num2str(ConcI) ', Conc0: ' num2str(Conc0)]});
ylabel('Current nA')
title({['GHK Plot - Positive P: ' num2str(PositiveIonPermeability) ', Negative P: ' num2str(NegativeIonPermeability) ],[ ', Offset ' num2str(Offset*1e9) 'nA, error: ' num2str(e)],['Ratio: ' num2str(PositiveIonPermeability/NegativeIonPermeability)]});
%Plot the raw data - on top of the plot created by GHK_Total Current
plot(IV(:,2),IV(:,1),'o');
plot((V-Vm)*1e3,I(:,1),'+');
hold off

disp(['Positive P: ' num2str(PositiveIonPermeability) ', Negative P:' num2str(NegativeIonPermeability) ', Offset ' num2str(Offset*1e9) 'nA, Error: ' num2str(e)]);
disp(['Ratio ' num2str(PositiveIonPermeability/NegativeIonPermeability)]);
disp(['ConcI: ' num2str(ConcI) ', Conc0: ' num2str(Conc0)]);

    function [ CurrentS ] = GHK_CurrentByMultiIon_S( z, E , P, Cres, Ccap )
        an_expTerm = exp(-z.*F.*E.*(1./R).*(1./T));
        CurrentS = (P.*z.*z.*E.*F.*F).*(Ccap-(Cres.*an_expTerm)).*(1./(R.*T.*(1-an_expTerm)));
    end

    function [ I_Total, I_Components ] = GHK_TotalCurrent_Multi( z, E , P, ConcI, Conc0)
        
        n_1 = max(size(z));
        I_Components = zeros(max(size(E)),n_1);
        for i = 1:n_1
            if(i==1)
                n_z = 2;
            else
                n_z = 1;
            end
            I_Components(:,i) = GHK_CurrentByMultiIon_S( z(i), E , P(i), ConcI*abs(z(n_z)), Conc0*abs(z(n_z)) );
        end
        I_Total = sum(I_Components,2);
        
        DoPlot = 0;
        if(DoPlot >0)
            ORG = Matlab2OriginPlot();
            ORG.HoldOff();
            ORG.PlotScatter(E',I_Components(:,1)','Hi2', ORG.ColourPicker());
            ORG.PlotScatter(E',I_Components(:,2)','Hi2', ORG.ColourPicker());
            ORG.ylabel('Current','nA');
            ORG.xlabel('Voltage','mV');
            ORG.HideActiveWkBk()
            %ORG.xComment('');
            ORG.HoldOn()
            ORG.PlotScatter(E',I_Total','Hi2', ORG.ColourPicker());
            
            hold off;
            plot(E,I_Components,' - ');
            hold all;
            plot(E,I_Total);
            
            xlabel('Voltage V');
            ylabel('Current A')
            title(['GHK Plot']);
        end
        
    end

end

