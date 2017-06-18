function i=roulettewheelSelection(prob)

    r=rand;
    
    C=cumsum(prob);
    
    i=find(r<=C,1,'first');

end