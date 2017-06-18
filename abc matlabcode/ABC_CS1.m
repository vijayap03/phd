function [Best_Weight]=ABC_CS1(Pt1,Pt2,Pt3,Pt4)


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
scout_limit = 5;
Iteration = 10;
emp_bee = bee;
while cycle < Iteration
    
    %% Employ bee updation
    
    for i = 1:rsize
       % k =  randi([1 rsize],1,1);
        k=randint(1,1,[1 rsize])
        while  i == k
            %k =  randi([1 rsize],1,1);
            k=randint(1,1,[1 rsize])
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
        if Fit(i)< emp_fit(i)
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
        %k =  randi([1 rsize],1,1);
        k=randint(1,1,[1 rsize])
        while  i == k
            %k =  randi([1 rsize],1,1);
            k=randint(1,1,[1 rsize])
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
        if onlooker_fit(i)< emp_fit(i)
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
    
    
    Best_Weight = b(cycle,1:end-1);
    
    
    %% Trial Counter Check
    
    for i = 1:rsize
        if emp_fit(i) == onlooker_fit(i)
            trial = trial + 1;
            if trial > scout_limit
                onlooker_bee(i,:) = (minimum + (maximum-(minimum)))*rand(1);
                onlooker_bee(i,:) = get_cuckoos(onlooker_bee(i,:),Best_Weight,Lb,Ub);
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
%Bee = randi([1,c/2],[p,c]);
Bee=randint(p,c,[1 c/2])
end

function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n,
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    
    % In the next equation, the difference factor (s-best) means that
    % when the solution is the best solution, it remains unchanged.
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale;
    % otherwise, Levy flights may become too aggresive/efficient,
    % which makes new solutions (even) jump out side of the design domain
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
    % Apply simple bounds/limits
    nest(j,:)=simplebounds(s,Lb,Ub);
end
end
% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
% Apply the lower bound
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I)+1;

% Apply the upper bounds
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
% Update this new move
s=ns_tmp;
end