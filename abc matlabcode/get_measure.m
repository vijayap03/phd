function [papr,BER,SER] = get_measure(B,Pt1,Pt2,Pt3,Pt4,T)
% the number of possible phase factor combinations
k = 1;

final_signal = B(k,1)*Pt1+B(k,2)*Pt2+B(k,3)*Pt3+B(k,4)*Pt4;
[number,BER] = biterr(round(T*100),round(abs(real(final_signal))*100));
[number,SER] = symerr(round(T*100),round(abs(real(final_signal))*100));
meank = mean(abs(final_signal).^2);
peak = max(abs(final_signal).^2);
papr = 10*log10(peak/meank);

end