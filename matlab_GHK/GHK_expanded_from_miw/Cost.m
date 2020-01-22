function [ CostFunc  ] = Cost(I, GHK, P )

CostFunc  =  sum((I - GHK * [abs(P(1)) ; abs(P(2)); P(3)] ).^2 );

CostFunc = sqrt(CostFunc) / (max(I)- min(I));

end

