function z=costfunction(u,a1)
K=128;
    sum_u = sum(a1.*(repmat((u).',1,K)));
     u_Power = abs(sum_u.^2); 
        u_Peak_Power = max(u_Power,[],2); 
        u_Mean_Power = mean(u_Power,2);
        papr_u = 10*log10(u_Peak_Power./u_Mean_Power);
        z=papr_u;