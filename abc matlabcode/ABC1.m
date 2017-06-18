function [Best_Weight]=ABC1(Pt1,Pt2,Pt3,Pt4)


minimum = 1;
maximum = 4;
nd=4;
Pop_size= 10; % Population size
bee = CreatSolution(Pop_size,nd);
Lb=ones(1,nd);
% Upper bounds
Ub=4*ones(1,nd);

%% Fitness

for i = 1 : size(bee,1)
    [Fit(i,1)] = Fitness(bee(i,:),Pt1,Pt2,Pt3,Pt4);
end

[rsize,csize] = size(bee);

cycle = 1;
trial = 0;
best = [];
scout_limit = 3;
Iteration = 2;
emp_bee = bee;
while cycle < Iteration
    
    %% Employ bee updation
    
    for i = 1:rsize
        k =randint(1,1,[1 rsize]);
        while  i == k
           % k =  randi([1 rsize],1,1);
           k= randint(1,1,[1 rsize]);
        end
        for j =  1 : csize
            pi = (-1 + (1-(-1)))*rand(1);
            emp_bee(i,j) = bee(i,j) + (pi * (bee(i,j)-bee(k,j)))*0.01;
            
            if emp_bee(i,j)< minimum
                emp_bee(i,j) = minimum;
            elseif emp_bee(i,j)> maximum
                emp_bee(i,j) = maximum;
            end
        end
    end
    
    % Fitness
    for i = 1 : size(emp_bee,1)
        [emp_fit(i,1)] = Fitness(emp_bee(i,:),Pt1,Pt2,Pt3,Pt4);
    end
    
    % Greedy Selection Employed Bee
    
    for i = 1:rsize
        if Fit(i)> emp_fit(i)
            emp_bee(i,:) = bee(i,:);
            emp_fit(i) = Fit(i);
        end
    end
    
    
    %% Probability Calculation
    
    
    for i = 1:rsize
        prob(i) = emp_fit(i)/sum(emp_fit);
    end
    
    [p_val,p_ind] = sort(prob,'ascend');
    
    
    %% Onlooker Bee Phase
    
    % select solution depending on probability
    
    for i = 1:rsize
        onlooker_bee(i,:) = emp_bee(p_ind(i),:);
    end
    
    % Produce new solution
    
    for i = 1:rsize
       % k =  randi([1 rsize],1,1);
        k=randint(1,1,[1 rsize]);
        while  i == k
           % k =  randi([1 rsize],1,1);
            k=randint(1,1,[1 rsize]);
        end
        
        for j =  1 : csize
            pi = (-1 + (1-(-1)))*rand(1);
            
            onlooker_bee(i,j) = emp_bee(i,j) + (pi * (emp_bee(i,j)-emp_bee(k,j)))*0.01;
            
            if onlooker_bee(i,j)< minimum
                onlooker_bee(i,j) = minimum;
            elseif onlooker_bee(i,j)> maximum
                onlooker_bee(i,j) = maximum;
            end
        end
        
    end
    
    % fitness
    for i = 1 : size(emp_bee,1)
        [onlooker_fit(i,1) ] = Fitness(onlooker_bee(i,:),Pt1,Pt2,Pt3,Pt4);
    end
    
    % Apply greedy selection
    
    for i = 1:rsize
        if onlooker_fit(i)> emp_fit(i)
            onlooker_bee(i,:) = emp_bee(i,:);
            onlooker_fit(i)= emp_fit(i);
        end
    end
    
    %% Memorizing Best
    
    [val ind] = sort(onlooker_fit,'descend');
    if cycle == 1
        Gval = val(1,1);
        bst_bee = onlooker_bee(ind(1,1),:);
        b(cycle,:) = [bst_bee Gval];
    else
        if val(1,1) <   Gval
            Gval = val(1,1);
            bst_bee = onlooker_bee(ind(1,1),:);
            b(cycle,:) = [bst_bee Gval];
        else
            b(cycle,:) = b(cycle-1) ;
        end
    end
    
    
    Best_Weight = b(cycle,1:end-1)/2;
    
    
    %% Trial Counter Check
    
    for i = 1:rsize
        if emp_fit(i) == onlooker_fit(i)
            trial = trial + 1;
            if trial > scout_limit
                onlooker_bee(i,:) = (minimum + (maximum-(minimum)))*rand(1);
            end
        else
            trial = 0;
        end
    end
    
    bee = onlooker_bee;
    Fit = onlooker_fit;
    cycle = cycle +1;
end
end
function [papr] = Fitness(Bee,Pt1,Pt2,Pt3,Pt4)
% the number of possible phase factor combinations
p=[1 -1 j -j];
Bee = round(Bee);
minimum = 1;
maximum = 4;
for k = 1:length(Bee)
if Bee(1,k)< minimum
    Bee(1,k) = minimum;
elseif Bee(1,k)> maximum
    Bee(1,k) = maximum;
end
end
B = p(Bee);

k = 1;
final_signal = B(k,1)*Pt1+B(k,2)*Pt2+B(k,3)*Pt3+B(k,4)*Pt4;
meank = mean(abs(final_signal).^2);
peak = max(abs(final_signal).^2);
papr = 10*log10(peak/meank);

end

function [Bee] = CreatSolution(p,c)
Bee = randint(p,c,[1 c/2]);
end
