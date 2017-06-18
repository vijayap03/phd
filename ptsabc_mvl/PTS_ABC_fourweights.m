clc; 
clear all; 
close all;  
V = 4;
N = 64;
P =8;
%fft_len = 4% FFt size  
QPSK_Set = [1 -1 j -j];  
Choose_Len =64; 
p=[1 -1 j -j]; % phase factor possible values
B1=[];
for b1=1:4
    for b2=1:4
        for b3=1:4
                                      
                B1 = [B1; [p(1)  p(b1)  p(b2) p(b3) ]];
                % all possible combinations
            end
        end
end
MAX_SYMBOLS = 100;
PAPR_PTSadqk = zeros(1,MAX_SYMBOLS);
framelen = 128;
for nSymbol=1:MAX_SYMBOLS; 
    K = 128; %SIZE OF OFDM Symbol   
    data = randint(1,K,4)+1;
    data_mod = QPSK_Set(data);
S(:,1:2:framelen) = data_mod(:,1:2:framelen);   
S(:,2:2:framelen) = -conj(data_mod(:,2:2:framelen)); 
S_t = S;
%original papr
ofdm_signal=ifft(S_t,[],2);
abs_ofdm = abs(ofdm_signal.^2); 
        Peak_Power_ofdm= max(abs_ofdm,[],2); 
        Mean_Power_ofdm = mean(abs_ofdm,2);
        PAPR_OFDM(:,nSymbol) = 10*log10(Peak_Power_ofdm./Mean_Power_ofdm);
m =1; 
Y = reshape(S_t,16,8);
A = zeros(N,P);
for i = 1:4:61;
   A(i,:)=[Y(m,1:2) zeros(1,6)];
   A(i+1,:)=[zeros(1,2) Y(m,3:4) zeros(1,4)];
   A(i+2,:)=[zeros(1,4) Y(m,5:6) zeros(1,2)];
   A(i+3,:)=[zeros(1,6) Y(m,7:8)];
    m =m+1;
end
A1 = reshape(A,4,128);
a1 = ifft(A1,[],2);
size(a1);
min_value = 10;  
for e =1:Choose_Len;  
temp_phase = (B1(e,:)).';
%size(temp_phase);
temp_max = max(abs(sum(a1.*repmat(temp_phase,1,K))));
if temp_max<min_value;  
min_value = temp_max;  
Best_e = e;
end
end
aa = sum(a1.*repmat((B1(Best_e,:)).',1,K));
Signal_Poweradqk = abs(aa.^2); 
        Peak_Poweradqk= max(Signal_Poweradqk,[],2); 
        Mean_Poweradqk = mean(Signal_Poweradqk,2);
        PAPR_PTSadqk(nSymbol) = 10*log10(Peak_Poweradqk./Mean_Poweradqk);
%end
%%%%%%%%%%%%%%%%%%%%%%%%%ABC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
nvar=4;
varsize=[1 nvar];
varmin=[-j,-j,-j,-j];
varmax=[j,j,j,j];
MaxIt=10;%%%%%%%%%generations=MaxIt%%%%%in this program%%%%%%
npop=8;
nonlooker=npop;
L=round(0.6*nvar*npop);
a_A=1;
empty_bee.position=[];
empty_bee.cost=[];
pop_A=repmat(empty_bee,npop,1);
bestsol.cost=PAPR_OFDM(nSymbol);
for i=1:npop;
    index=(randint(1,1,[1,64]));
     pop_A(i).position=B1(index,:);
    %pop_A(i).position=2*(randn(1,nvar,[varmin,varmax]))-1;
    pop_A(i).cost=costfunction(pop_A(i).position,a1);
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
        phi=a_A*unifrnd(-1,+1,varsize);
        %newbee.position=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position);
        new=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position);
        newbee.position=simplebounds(new,varmin,varmax);
             newbee.cost=costfunction(newbee.position,a1);
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
        phi=a_A*unifrnd(-1,+1,varsize);
        %newbee.position=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position)
        new=pop_A(i).position+phi.*(pop_A(i).position-pop_A(k).position);
        newbee.position=simplebounds(new,varmin,varmax);        
        newbee.cost=costfunction(newbee.position,a1);
        if newbee.cost<=pop_A(i).cost
            pop_A(i)=newbee;
        else
            C(i)=C(i)+1;
        end
    end
    for i=1:npop
        if C(i)>=L
            
      % pop_A(i).position=unifrnd(randint(varmin,varmax,varsize));
       index=(randint(1,1,[1,64]));
     pop_A(i).position=B1(index,:);
          pop_A(i).cost=costfunction(pop_A(i).position,a1);     
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
final_bestcost_ABC(nSymbol)=min(bestcost);
end    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ABC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cdf1,PAPR1]=ecdf(PAPR_PTSadqk);
 [cdf2, PAPR2] = ecdf(final_bestcost_ABC);
 [cdf3,PAPR3]=ecdf(PAPR_OFDM);
 semilogy(PAPR1,1-cdf1,'-b',PAPR2,1-cdf2,'-g',PAPR3,1-cdf3,'-r','linewidth',2);
% semilogy(PAPR2,1-cdf2,'-c','linewidth',2);
  %disp(max(PAPR2));
  %semilogy(PAPR2,1-cdf2,'-g','linewidth',2); 
  legend('pts','pts ABC','ofdm');
% % title('MIMO ADJACENT PTS WITH DIFFERENT MODULATIONS');  
 xlabel('PAPR0 [dB]');  
  xlim([1 10]);
  ylabel('CCDF2 (Pr[PAPR>PAPR0])');  
 grid on;
 figure;
 plot(bestcost,'Linewidth',2);
 
 
 
 
 
 
 
 
 %----------------------W=4-----------------------------
 

        
  
 
 
 
 
 
