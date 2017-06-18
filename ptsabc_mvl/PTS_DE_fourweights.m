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
        PAPR_OFDM(:,nSymbol) = 10*log10(Peak_Power_ofdm./Mean_Power_ofdm)
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
size(a1)
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
        PAPR_PTSadqk(nSymbol) = 10*log10(Peak_Poweradqk./Mean_Poweradqk)
%end
     
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11-----------------------%de%---------------------------
G=1;
NP=4;
F=1;
CR=0.8;
best_cost=PAPR_OFDM(nSymbol);    
              
 for j=1:NP
     index=(randint(1,1,[1,64]));
     pop(j,:)=B1(index,:);
     sum_pop = sum(a1.*repmat((pop(j,:)).',1,K));
     pop_Power = abs(sum_pop.^2); 
        pop_Peak_Power = max(pop_Power,[],2); 
        pop_Mean_Power = mean(pop_Power,2);
        papr_pop(j) = 10*log10(pop_Peak_Power./pop_Mean_Power)
        if papr_pop(j)< best_cost
            best_cost=papr_pop(j);
        end
  end
%best_cost=PAPR_OFDM;
V=[];
H=[];
for G=1:1
for i=1:NP
R=randperm(NP)
R(R==i)=[];
a=R(1);
b=R(2);
c=R(3);
  V(i,:)=B1(a,:)+F.*(B1(b,:)-B1(c,:))
 for j= 1:4
     if (V(i,j)>=exp(j*pi/4) )&& (V(i,j)<exp(3*j*pi/4))
         V(i,j)=j;
     else if V(i,j)>=exp(3*j*pi/4)&& V(i,j)<exp(5*j*pi/4)
         V(i,j)=-1;
         else if V(i,j)>=exp(5*j*pi/4)&& V(i,j)<exp(7*j*pi/4)
               V(i,j)=-j;
             else
                 V(i,j)=1;
             end
         end
     end
 
  end
 %disp('mutation value :');
 %disp(V(i,:));
 r=rand(1);
 %disp('random value is');
 %disp(r);
 if  r< CR
     U(i,:)=V(i,:);
else
     U(i,:)=pop(i,:);
%     disp('solution space value is:');
%     disp(U(i));
 end
 
 sum_U = sum(a1.*repmat((U(i,:)).',1,K));
U_Power = abs(sum_U.^2); 
        U_Peak_Power = max(U_Power,[],2); 
        U_Mean_Power = mean(U_Power,2);
        papr_U(i) = 10*log10(U_Peak_Power./U_Mean_Power)

        
 if papr_U(i)<papr_pop(i)
     pop(i,:)=U(i,:)
     if pop(i)< best_cost
         best_cost=papr_pop(i)
     end
 end
 H(i,:)=pop(i,:);
sum_H = sum(a1.*repmat((H(i,:).'),1,K));
H_Power= abs(sum_H.^2); 
        H_Peak_Power = max(H_Power,[],2); 
        H_Mean_Power = mean(H_Power,2);
        papr_H(i) = 10*log10(H_Peak_Power./H_Mean_Power)
end

disp('final papr is');
final_de(G)=min(papr_H) 
end
final_de_1(nSymbol)=min(final_de)
end
% %disp(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cuckoo%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
































[cdf1,PAPR1]=ecdf(PAPR_PTSadqk);
 [cdf2, PAPR2] = ecdf(final_de_1);
 [cdf3,PAPR3]=ecdf(PAPR_OFDM);
 semilogy(PAPR1,1-cdf1,'-b',PAPR2,1-cdf2,'-g',PAPR3,1-cdf3,'-r','linewidth',2);
% semilogy(PAPR2,1-cdf2,'-c','linewidth',2);
  %disp(max(PAPR2));
  legend('pts','pts_de','papr_ofdm');
% % title('MIMO ADJACENT PTS WITH DIFFERENT MODULATIONS');  
 xlabel('PAPR0 [dB]');  
  xlim([1 10]);
  ylabel('CCDF2 (Pr[PAPR>PAPR0])');  
 grid on;
 
 
 
 
 
 
 
 
 
 %----------------------W=4-----------------------------
 

        
  
 
 
 
 
 
