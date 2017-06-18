 function [papr] = Fitness(Spider,a1)
% the number of possible phase factor combinations
p=[1 -1 j -j];
Spider = round(Spider);
minimum = 1;
maximum = 4;
for k = 1:length(Spider)
if Spider(1,k)< minimum
    Spider(1,k) = minimum;
elseif Spider(1,k)> maximum
    Spider(1,k) = maximum; 
end
end
B = p(Spider);
%%%%%%%%%
K=128;
aa1 = sum(a1.*repmat((B).',1,K));
Signal_Powersso = abs(aa1.^2); 
        Peak_Powersso= max(Signal_Powersso,[],2); 
        Mean_Powersso = mean(Signal_Powersso,2);
        papr = 10*log10(Peak_Powersso./Mean_Powersso);%%%%%%

end
      