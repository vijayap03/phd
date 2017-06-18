function BS1=BestSpider(a1,B1)
% comment added by sp 
spidn=10;
itern=2;
Pop_size = 10;
nd = 4;
L = 10;
f=@Griewank;
xd=-1;
xu=4;
lb = ones(Pop_size,nd);
ub = 4*ones(Pop_size,nd);
%% Initial Parameters
rand('state',0');  % Reset the random generator
% Define the poblation of females and males
fpl = 0.65;     % Lower Female Percent
fpu = 0.9;      % Upper Female Percent
%fp = fpl+(fpu-fpl)*rand;
fp = fpu-(fpu-fpl)*rand;	% Aleatory Percent
fn = round(spidn*fp);   % Number of females
mn = spidn-fn;          % Number of males
dims = 4;
%Probabilities of attraction or repulsion
% Proper tunning for better results
pm = exp(-(0.1:(3-0.1)/(itern-1):3));
% Initialization of vectors
fsp = zeros(fn,dims);   % Initlize females
msp = zeros(mn,dims);   % Initlize males
fefit = zeros(fn,1);    % Initlize fitness females
mafit = zeros(mn,1);    % Initlize fitness males
spwei = zeros(spidn,1); % Initlize weigth spiders
fewei = zeros(fn,1); % Initlize weigth spiders
mawei = zeros(mn,1); % Initlize weigth spiders
%% Population Initialization
% Generate Females
%fsp = CreatSolution(Pop_size,nd);
fsp = CreatSolution(fn,nd);
% Generate Males
msp = CreatSolution(mn,nd);
%% **** Evaluations *****
% Evaluation of function for females
for i2=1:size(fsp,1)
    fefit(i2,1)= Fitness(fsp(i2,:),a1);
end

% Evaluation of function for males

for i3=1:size(msp,1)
    mafit(i3,1)= Fitness(msp(i3,:),a1);
end

%% ***** Assign weigth or sort ***********
% Obtain weight for every spider

spfit = [fefit' mafit']';   % Mix Females and Males
bfitw = min(spfit);          % best fitness
wfit = max(spfit);          % worst fitness
for i4=1:spidn
    spwei(i4) = 0.001+((spfit(i4)-wfit)/(bfitw-wfit));
end
fewei = spwei(1:fn);      % Separate the female mass
mawei = spwei(fn+1:spidn);% Separate the male mass
%% Memory of the best
% Check the best position
[maxval,Ibe] = max(spwei);
% Check if female or male
if Ibe > fn
    % Is Male
    spbest=msp(Ibe-fn,:);   % Asign best position to spbest
    bfit = mafit(Ibe-fn);      % Get best fitness for memory
else
    % Is Female
    spbest=fsp(Ibe,:);      % Asign best position to spbest
    bfit = fefit(Ibe);      % Get best fitness for memory
end

%% Start the iterations

for i5=1:itern
    %% ***** Movement of spiders *****
    % Move Females
    [fsp] = FeMove(spidn,fn,fsp,msp,spbest,Ibe,spwei,dims,lb,ub,pm(i5));
    % Move Males
    [msp] = MaMove(fn,mn,fsp,msp,fewei,mawei,dims,lb,ub,pm(i5));
    %% **** Evaluations *****
    % Evaluation of function for females
    for j1=1:fn
       % fefit(j1)=f(fsp(j1,:),dims);
       fefit(j1)=Fitness(fsp(j1,:),a1);
    end
    % Evaluation of function for males
    for j2=1:mn
        %mafit(j2)=f(msp(j2,:),dims);
        mafit(j2)=Fitness(msp(j2,:),a1);
    end
    %% ***** Assign weigth or sort ***********
    spfit = [fefit' mafit']';   % Mix Females and Males
    bfitw = min(spfit);          % best fitness
    wfit = max(spfit);          % worst fitness
    % Obtain weight for every spider
    for j3=1:spidn
        spwei(j3) = 0.001+((spfit(j3)-wfit)/(bfitw-wfit));
    end
    fewei = spwei(1:fn);      % Separate the female mass
    mawei = spwei(fn+1:spidn);% Separate the male mass
    %% Mating Operator
    [ofspr] = Mating(fewei,mawei,fsp,msp,dims);
    %% Selection of the Mating
    if isempty(ofspr)
        %             % Do nothing
    else
        %[fsp,msp,fefit,mafit] = Survive(fsp,msp,ofspr,fefit,mafit,spfit,f,fn,dims);
        [fsp,msp,fefit,mafit] = Survive(fsp,msp,ofspr,fefit,mafit,spfit,fn,a1)
        % ***** Recalculate the weigth or sort ***********
        spfit = [fefit' mafit']';   % Mix Females and Males
        bfitw = min(spfit);          % best fitness
        wfit = max(spfit);          % worst fitness
        % Obtain weight for every spider
        for j4=1:spidn
            spwei(j4) = 0.001+((spfit(j4)-wfit)/(bfitw-wfit));
        end
        fewei = spwei(1:fn);      % Separate the female mass
        mawei = spwei(fn+1:spidn);% Separate the male mass
    end
    %% Memory of the best
    % Check if best position belongs to male or female
    [sp, Ibe2] = max(spwei);
    if Ibe2 > fn
        % Is Male
        spbest2=msp(Ibe2-fn,:);      % Asign best position to spbest
        bfit2 = mafit(Ibe2-fn);      % Get best fitness for memory
    else
        % Is Female
        spbest2 = fsp(Ibe2,:);  % Asign best position to spbest
        bfit2 = fefit(Ibe2);    % Get best fitness for memory
    end
    %% Global Memory
    if bfit<=bfit2
        bfit = bfit;
        spbest = spbest;      % Asign best position to spbest
        befit(i5) = bfit;
    else
        bfit = bfit2;
        spbest = spbest2;      % Asign best position to spbest
        befit(i5) = bfit;
    end
    spbesth(i5,:)=spbest;
end
BS1 = spbest;