clc
clear all
close all
warning off

%% Initialize Data

 M=4;
NN = 100;       % Number of iterations
N = 16;         % number of subcarriers
L = 4;
r=randint(NN,N,[0 M-1]);
dm=pskmod(r,M);
% oversampling factor

%% Initialize the process
p=[1 -1 i -i]; % phase factor possible values
B=[];
U=100;
B=randsrc(U,N,p);
for i1=1:NN        
     x1f=ifft([dm(i1,1:N/2) zeros(1,(L-1)*N) dm(i1,N/2+1:N)]);    
     mean1=mean(abs(x1f).^2);
     peak1=max(abs(x1f).^2);
     paprofdm(i1)=10*log10(peak1/mean1);

PP1=[];
U1=16;
 for k=1:U
        PP1(k,:)=(B(k,:).*dm(i1,:));  
        pt1(k,:)=abs(ifft([PP1(k,1:N/2) zeros(1,(L-1)*N) PP1(k,N/2+1:N)]));
        ppr1(k)=10*log10(max(abs(pt1(k,:)).^2)/mean(abs(pt1(k,:)).^2));
 end     
    paprslm(i1)=min(ppr1);   

%%%%%%%abc%%%%%%%%%%%%%%
nvar=16;
varsize=[1 nvar];
for i2=1:N
varmin(1,i2)=-j;
varmax(1,i2)=j;
end
MaxIt=2;%%%%%%%%%generations=MaxIt%%%%%in this program%%%%%%
npop=8;
nonlooker=npop;
L1=round(0.6*nvar*npop);
a_A=1;
empty_bee.position=[];
empty_bee.cost=[];
pop_A=repmat(empty_bee,npop,1);
bestsol.cost=paprofdm(i1);
for i=1:npop;
    index=(randint(1,1,[1,100]));
     pop_A(i).position=B(index,:);
    %pop_A(i).position=2*(randn(1,nvar,[varmin,varmax]))-1;
    pop_A(i).cost=costfunctionslm(pop_A(i).position,dm(i1,:));
    if pop_A(i).cost<bestsol.cost
        bestsol=pop_A(i);
    end
end
C=zeros(npop,1);
bestcost=zeros(MaxIt,1);
%%%%%mainloop%%%%%%
for it=1:MaxIt
    for i=1:npop;
        k=[1:i-1 i+1:npop];
        k=randint(1,1,[1,numel(k)]);
        phi=a_A*unifrnd(-1,+1,[1 16]);
        %newbee.position=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position);
        new=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position);
        newbee.position=simplebounds(new,varmin,varmax);
             newbee.cost=costfunctionslm(newbee.position,dm(i1,:));
        if newbee.cost<=pop_A(i).cost
            pop_A(i)=newbee;
        else
            C(i)=C(i)+1;
        end
    end
    F=zeros(npop,1);
    meancost=mean([pop_A.cost]);
    for i=1:npop;
        F(i)=exp(-pop_A(i).cost/meancost);
    end
    prob=F/sum(F);
    for m=1:nonlooker
        i=roulettewheelselection(prob);
        k=[1:i-1 i+1:npop];
        k=randint(1,1,[1,numel(k)]);
        phi=a_A*unifrnd(-1,+1,[1 16]);
        %newbee.position=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position)
        new=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position);
        newbee.position=simplebounds(new,varmin,varmax);        
        newbee.cost=costfunctionslm(newbee.position,dm(i,:));
        if newbee.cost<=pop_A(i).cost
            pop_A(i)=newbee;
        else
            C(i)=C(i)+1;
        end
    end
    for i=1:npop
        if C(i)>=L1
            
      % pop_A(i).position=unifrnd(randint(varmin,varmax,varsize));
       index=(randint(1,1,[1,64]));
     pop_A(i).position=B1(index,:);
          pop_A(i).cost=costfunctionslm(pop_A(i).position,dm(i,:));     
       C(i)=0;
        end
    end
    for i=1:npop
        if pop_A(i).cost<=bestsol.cost
            bestsol=pop_A(i);
            
        end
    end
    bestcost(it)=bestsol.cost;
   % bestcost(MaxIt)=min([pop_A.cost]);
end
final_bestcost_ABC(i1)=min(bestcost);
end    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ABC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure,
hold on
[cy,cx] = cdfcalc(paprofdm);
%[cy,cx] = cdfcalc(paprtds)
semilogy(flipud(cx),cy(2:end),'r')

[cy,cx] = cdfcalc(paprslm);
semilogy(flipud(cx),cy(2:end),'b')
[cy,cx] = cdfcalc(paprofdm);
[cy,cx] = cdfcalc(final_bestcost_ABC)
semilogy(flipud(cx),cy(2:end),'g')
grid on
title('PAPR');
legend('ofdm','slm','slmabc');

% 
