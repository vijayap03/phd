clc
clear all
close all
warning off
dbstop if error

%% Initialize Data

load ofdm % QPSK modulated, 64-element OFDM symbols.

%% All permutations of phase factor B

NN = 40;         % SNR
N = 512;         % number of subcarriers
L = 4;           % oversampling factor
i = 10;

%% Initialize the process

for k = 1:NN
    
    O1 = ofdm_symbol(i,1:32);
    O2 = zeros(1,(L-1)*N);
    O3 = ofdm_symbol(i,33:64);
    time_domain_signal=abs(ifft([O1 O2 O3]));
    Tds = time_domain_signal+k;
    
    %% Partition OFDM Symbol
    
    P1 = [ofdm_symbol(i,1:16) zeros(1,48)];
    P2 = [zeros(1,16) ofdm_symbol(i,17:32) zeros(1,32)];
    P3 = [zeros(1,32) ofdm_symbol(i,33:48) zeros(1,16)];
    P4 = [zeros(1,48) ofdm_symbol(i,49:64)];
    
    %% Transform Pi to Time Domain
    
    Pt1 = abs(ifft([P1(1:32) zeros(1,(L-1)*N) P1(33:64)]));
    Pt2 = abs(ifft([P2(1:32) zeros(1,(L-1)*N) P2(33:64)]));
    Pt3 = abs(ifft([P3(1:32) zeros(1,(L-1)*N) P3(33:64)]));
    Pt4 = abs(ifft([P4(1:32) zeros(1,(L-1)*N) P4(33:64)]));
    Pt = Pt1+Pt2+Pt3+Pt4;
    
    %% Measure
    
    [Bestbee] = ABC1(Pt1,Pt2,Pt3,Pt4);
    [Bestnest] = CS(Pt1,Pt2,Pt3,Pt4);
    [Best] = ABC_CS1(Pt1,Pt2,Pt3,Pt4);
    [papr(k),BER(k),SER(k)] = get_measure(Best,Pt1,Pt2,Pt3,Pt4,Pt);
    [papr2(k),BER2(k),SER2(k)] = get_measure(Bestbee,Pt1,Pt2,Pt3,Pt4,Pt);
    [papr1(k),BER1(k),SER1(k)] = get_measure(Bestnest,Pt1,Pt2,Pt3,Pt4,Pt);
    
end

figure,
hold on
[cy,cx] = cdfcalc(papr2);
semilogy(flipud(cx),cy(2:end),'b')
grid on
[cy,cx] = cdfcalc(papr1);
semilogy(flipud(cx),cy(2:end),'r')
[cy,cx] = cdfcalc(papr);
semilogy(flipud(cx),cy(2:end),'g')
title('PAPR');
legend('With ABC','With CS','With ABC CS');

BER1 = sort(BER1,'descend');
BER = sort(BER,'descend');
BER2 = sort(BER2,'descend');
SER1 = sort(SER1,'descend');
SER2 = sort(SER2,'descend');
SER = sort(SER,'descend');

range = 1:NN;
figure,
semilogy(range,BER,'b');
hold on
grid on
semilogy(range,BER1,'r');
semilogy(range,BER2,'g');
title('Bit Error Rate');
xlabel('SNR');
ylabel('BER');
legend('With ABC','With CS','With ABC CS');
hold off

range = 1:NN;
figure,
semilogy(range,SER,'b');
hold on
grid on
semilogy(range,SER1,'r');
semilogy(range,SER2,'g');
title('Symbol Error Rate');
xlabel('SNR');
ylabel('SER');
legend('With ABC','With CS','With ABC CS');
hold off

