function [fsp,msp,fefit,mafit] = Survive(fsp,msp,ofspr,fefit,mafit,spfit,fn,a1)
%SURVIVE Summary of this function goes here
%   Detailed explanation goes here
    [n1, sp1] = size(ofspr);
    %Evalute the offspring
    for q1=1:n1
            offit(q1)=Fitness(ofspr(q1,:),a1);
    end
    for i=1:n1
        %Calculate the worst spider
        [w1, w2]=max(spfit);
        %If the offspring is better than the worst spider
        if offit(i)<w1
            %Check if is male or female
            if w2>fn
                %Male
                msp(w2-fn,:)=ofspr(i,:);
                mafit(w2-fn)=offit(i);
            else
                %Female
                fsp(w2,:)=ofspr(i,:);
                fefit(w2)=offit(i);
            end
            spfit(w2)=offit(i);
        end  
    end
end