function [ofsp] = Mating(femass,mamass,fsp,msp,dims)
%CROSSOVER Summary of this function goes here
%   Detailed explanation goes here
    % Generate the offsprings
    ofsp = [];
    cont = 1;
    % Check whether a spider is good or not (above median)
    [Indf,pvkd3] = find(femass);                % Female spiders
    [Indm,pvkd4] = find(mamass>median(mamass)); % Male spiders above median
    fespid=fsp(Indf,:);             % Female spiders
    maspid=msp(Indm,:);             % Only the Male spiders above median
    sp2mate = [];
    %% Calculate the radio
	rad = zeros(1,dims);
    spid = [fsp' msp']';
    for i=1:dims
        rad(i) = max(spid(:,i))-min(spid(:,i));
    end
    r=(sum(rad)/2)/(dims);
    %% Start looking if there's a good female near
    [sz,pvkd5] = size(Indf);
    dist=zeros(1,sz);
    for i=1:size(Indm)
        iaux = 1;       % Aux to form the elements to mate
        for q1=1:size(Indf)
            dist(q1)=norm(msp(Indm(i),:)-fsp(Indf(q1),:));
        end
        for k=1:size(Indf)
            if dist(k)<r
                mate(iaux,:) = fsp(Indf(k),:);
                mass(iaux) = femass(Indf(k));
                iaux = iaux+1;
            % Form the matrix with elements to mate
            sp2mate = [msp(Indm(i),:)' mate']';
            masmate = [mamass(Indm(i)) mass];                
            end
        end
        % Realizo el mate
        if isempty(sp2mate)
            % do nothing
        else
            [num2,n] = size(sp2mate);
            for k=1:num2
            for q1=1:n
                accumulation = cumsum(masmate);
                p = rand() * accumulation(end);
                chosen_index = -1;
                for index = 1 : length(accumulation)
                    if (accumulation(index) > p)
                        chosen_index = index;
                        break;
                    end
                end
                choice = chosen_index;
                % Forma the new element
                ofsp(k,q1)=sp2mate(choice,q1);
            end
            end
            cont = cont+1;
        end
    end
end