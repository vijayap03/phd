function papr=paprsso(S_t_sso)
N = 64;
P =8;
Y = reshape(S_t_sso,16,8);
A = zeros(N,P);
m=1;
for i1 = 1:4:61;
   A(i1,:)=[Y(m,1:2) zeros(1,6)];
   A(i1+1,:)=[zeros(1,2) Y(m,3:4) zeros(1,4)];
   A(i1+2,:)=[zeros(1,4) Y(m,5:6) zeros(1,2)];
   A(i1+3,:)=[zeros(1,6) Y(m,7:8)];
    m =m+1;
end    
 A1 = reshape(A,4,128);
a1 = ifft(A1,[],2);
a1_SP = a1;  
 BS_SP=BestSpider(a1_SP);
papr =Fitness(BS_SP,a1_SP);
end