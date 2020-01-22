function [total_ratio,anion_cation,conc] = sauvage_GHKpermeability(Caps)

DB = DBConnection;
E = Experiments(DB);
total_ratio = [];
conc = [];
AllowableSuppressionCodes = [0 16];
anion_cation = [];

for i = 1:length(Caps)
    cID = Caps(i);
    disp(['Capillary: ' num2str(cID) '. '])
    E.SELECT(['Capillary = ' num2str(cID) ' AND Suppressed = 0 AND Sealed > 0']);
    cap_conc = E.getCapillaryConc;
    IVs = 1;
    while IVs
        anion = E.getnPerm;
        cation = E.getpPerm;
        conc = [conc, E.getReservoirConc];
        if(anion == 0 & cation == 0)
        else
            ratio = cation/anion;
            total_ratio = [total_ratio ; ratio];% r diameter];
            anion_cation = [anion_cation; anion,cation];
        end
        IVs = E.NextResult;
    end
end

%%Altered to remove like for like values
idx = find(abs(conc - cap_conc)<=0.00001);

total_ratio(idx) = [];
anion_cation(idx,:) =[];
conc(idx)=[];

%

E.CloseConnection();

% Range = [-13 14];
% edges = logspace(-13,14,100);
% edgeLables = linspace(-13,14,100)';
%
%
% g = histc(total_ratio(:,1),edges);
%
% G = g/sum(g);
%
% figure(1);hold off;
% Data = G';
% subplot(1,2,1)
% logXHistogram(total_ratio(:,1),100,'r',Range);
%
% subplot(1,2,2); hold off;
% logXHistogram(total_ratio(:,1),60,'r');
%
% figure(2); hold off;
% Gc = cumsum(G);
% hold all;
% plot(edgeLables ,Gc,'r');
%
% figure(3)
% subplot(1,2,1); hold off;
% logXHistogram(anion_cation(:,1),50,'r');
%
% subplot(1,2,2); hold off;
% logXHistogram(anion_cation(:,2),50,'r');

end