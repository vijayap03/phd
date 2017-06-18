function z=costfunctionslm(u,a1)
PP1=(u.*a1);
N=16;
L=4;
pt1=abs(ifft([PP1(1,1:N/2) zeros(1,(L-1)*N) PP1(1,N/2+1:N)]));
        z=10*log10(max(abs(pt1).^2)/mean(abs(pt1).^2));
end