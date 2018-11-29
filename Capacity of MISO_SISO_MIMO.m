% Reference Book David Tse and Pramod Viswanath(2005)
% Page 226, Figure 5.17

clc;
close all;
clear all;
% 4 Transmitters and 4 Receivers
Tx = 4;
Rx = 4;
% Symbol Periods
L = 4;     
% Signal to Noise Ratio
SNR = [-10:20];
snr = 10.^(SNR/10);

%% 1. SIMO(1 Transmitter and 4 Receivers)
Capacity_SIMO = zeros(1,length(snr));
% Loop to run till the total symbol periods
for i = 1:L
% Channel Matrix 
h_SIMO = (randn(Rx,1)+1j*randn(Rx,1))/sqrt(2);
% Loop to run for all SNR values
for K = 1:length(snr)
Capacity_SIMO(K) = Capacity_SIMO(K) + log2(1+ snr(K)*(h_SIMO'*h_SIMO));
end
end
% Channel Capacity
Capacity_SIMO = Capacity_SIMO/L;

%% 2. MISO(4 Transmitters and 1 Receiver) without knowledge of CSI
Capacity_MISO_without_CSI = zeros(1,length(snr));
% Loop to run till the total symbol periods
for ite = 1:L
% Channel Matrix     
h_MISO = (randn(1,Tx)+1j*randn(1,Tx))/sqrt(2);
% Loop to run for all SNR values
for K = 1:length(snr)
Capacity_MISO_without_CSI(K) = Capacity_MISO_without_CSI(K) + log2(1+ (snr(K)*(h_MISO*h_MISO'))/(Tx*Tx)); % For Nt=Tx i.e. 4
end
end
% Channel Capacity
Capacity_MISO_without_CSI = Capacity_MISO_without_CSI/L;

%% 3. MISO(4 Transmitters and 1 Receiver) with knowledge of CSI
Capacity_MISO_with_CSI = zeros(1,length(snr));
% Loop to run till the total symbol periods
for ite = 1:L
% Channel Matrix     
h_MISO = (randn(1,Tx)+1j*randn(1,Tx))/sqrt(2);
% Loop to run for all SNR values
for K = 1:length(snr)
Capacity_MISO_with_CSI(K) = Capacity_MISO_with_CSI(K) + log2(1+ snr(K)*(h_MISO*h_MISO')/Tx);
end
end
% Channel Capacity
Capacity_MISO_with_CSI = Capacity_MISO_with_CSI/L;

%% 4. MIMO(4 Transmitters and 4 Receivers) with Equal Power Distribution
Capacity_MIMO_norm4 = zeros(1,length(snr));
Capacity_MIMO_norm3 = zeros(1,length(snr));
Capacity_MIMO_norm2 = zeros(1,length(snr));
Capacity_MIMO_norm1 = zeros(1,length(snr));
% Loop to run till the total symbol periods
for ite = 1:L
% Channel Matrix     
h_MIMO4 = (randn(Rx,Tx)+1j*randn(Rx,Tx))/sqrt(2);
h_MIMO3 = (randn(3,3)+1j*randn(3,3))/sqrt(2);
h_MIMO2 = (randn(2,2)+1j*randn(2,2))/sqrt(2);
h_MIMO1 = (randn(1,1)+1j*randn(1,1))/sqrt(2);
for K = 1:length(snr)
Capacity_MIMO_norm4(K) = Capacity_MIMO_norm4(K) + log2(det(eye(Rx)+snr(K)*(h_MIMO4)*(h_MIMO4)'/Tx));
Capacity_MIMO_norm3(K) = Capacity_MIMO_norm3(K) + log2(det(eye(3)+snr(K)*(h_MIMO3)*(h_MIMO3)'/3));
Capacity_MIMO_norm2(K) = Capacity_MIMO_norm2(K) + log2(det(eye(2)+snr(K)*(h_MIMO2)*(h_MIMO2)'/2));
Capacity_MIMO_norm1(K) = Capacity_MIMO_norm1(K) + log2(det(eye(1)+snr(K)*(h_MIMO1)*(h_MIMO1)'));
end
end
% Channel Capacity
Capacity_MIMO_norm4 = Capacity_MIMO_norm4/L;
Capacity_MIMO_norm3 = Capacity_MIMO_norm3/L;
Capacity_MIMO_norm2 = Capacity_MIMO_norm2/L;
Capacity_MIMO_norm1 = Capacity_MIMO_norm1/L;

%% 5. MIMO(4 Transmitters and 4 Receivers) with Water Filling Algorithm
% Reference Book Andrea Goldsmith(2005), Page 327
Capacity_MIMO_WFA4 = zeros(1,length(snr));
Capacity_MIMO_WFA3 = zeros(1,length(snr));
Capacity_MIMO_WFA2 = zeros(1,length(snr));
Capacity_MIMO_WFA1 = zeros(1,length(snr));
% Loop to run till the total symbol periods
for ite = 1:L
% Channel Matrix     
h_MIMO4 = (randn(Rx,Tx)+1j*randn(Rx,Tx))/sqrt(2);
h_MIMO3 = (randn(3,3)+1j*randn(1,1))/sqrt(2);
h_MIMO2 = (randn(2,2)+1j*randn(1,1))/sqrt(2);
h_MIMO1 = (randn(1,1)+1j*randn(1,1))/sqrt(2);
% SVD of the channel matrix 'h_MIMO'
[u4,E4,v4]=svd(h_MIMO4);
[u3,E3,v3]=svd(h_MIMO3);
[u2,E2,v2]=svd(h_MIMO2);
[u1,E1,v1]=svd(h_MIMO1);
%Gamma calculation for optimal power allocation
e4=diag(E4);
e3=diag(E3);
e2=diag(E2);
e1=diag(E1);
y4=sort(e4.*e4); 
y3=sort(e3.*e3); 
y2=sort(e2.*e2); 
y1=sort(e1.*e1); 
% 4X4
if (4/(1+(1/y4(1)+1/y4(2)+1/y4(3)+1/y4(4))))<y4(1)
yo4=4/(1+(1/y4(1)+1/y4(2)+1/y4(3)+1/y4(4)));
elseif (3/(1+(1/y4(2)+1/y4(3)+1/y4(4))))<y4(2)
yo4=3/(1+(1/y4(2)+1/y4(3)+1/y4(4)));
elseif (2/(1+(1/y4(3)+1/y4(4))))<y4(3)
yo4=2/(1+(1/y4(3)+1/y4(4)));
else
yo4=4/(1+(1/y4(4)));
end
% 3X3
if (3/(1+(1/y3(1)+1/y3(2)+1/y3(3))))<y3(1)
yo3=3/(1+(1/y3(1)+1/y3(2)+1/y3(3)));
elseif (3/(1+(1/y3(2)+1/y3(3))))<y3(2)
yo3=3/(1+(1/y3(2)+1/y3(3)));
else 
yo3=1/(1+(1/y3(3)));
end
% 2X2
if (2/(1+(1/y2(1)+1/y2(2))))<y2(1)
yo2=2/(1+(1/y2(1)+1/y2(2)));
else 
yo2=1/(1+(1/y2(2)));
end
% 1X1
yo1=4/(1+(1/y1(1)));
% Loop to run for all SNR values
for K = 1:length(snr)
Capacity_MIMO_WFA4(K) = Capacity_MIMO_WFA4(K) + log2(1+(snr(K)*y4(3))/yo4)+log2(1+(snr(K)*y4(4))/yo4);
Capacity_MIMO_WFA3(K) = Capacity_MIMO_WFA3(K) + log2(1+(snr(K)*y3(2))/yo3)+log2(1+(snr(K)*y3(3))/yo3);
Capacity_MIMO_WFA2(K) = Capacity_MIMO_WFA2(K) + log2(1+(snr(K)*y2(1))/yo2)+log2(1+(snr(K)*y2(2))/yo2);
Capacity_MIMO_WFA1(K) = Capacity_MIMO_WFA1(K) + log2(1+(snr(K)*y1(1))/yo1);
end
end
% Channel Capacity
Capacity_MIMO_WFA4 = Capacity_MIMO_WFA4/L;
Capacity_MIMO_WFA3 = Capacity_MIMO_WFA3/L;
Capacity_MIMO_WFA2 = Capacity_MIMO_WFA2/L;
Capacity_MIMO_WFA1 = Capacity_MIMO_WFA1/L;

%% Plot of various Systems
figure
% Ploting SIMO capacity
plot(SNR,Capacity_SIMO,'b')
hold on
% Ploting MISO capacity with CSI Information
plot(SNR,Capacity_MISO_with_CSI,'k')
hold on
% Ploting MISO capacity without CSI Information
plot(SNR,Capacity_MISO_without_CSI,'r')
hold on 
% Ploting MIMO capacity Equal Power Allocation(4Tx - 4Rx)
plot(SNR,real(Capacity_MIMO_norm4),'g-')
hold on 
% Ploting MIMO capacity Equal Power Allocation(3Tx - 3Rx)
plot(SNR,real(Capacity_MIMO_norm3),'g--')
hold on
% Ploting MIMO capacity Equal Power Allocation(2Tx - 2Rx)
plot(SNR,real(Capacity_MIMO_norm2),'g-.')
hold on
% Ploting MIMO capacity Equal Power Allocation(1Tx - 1Rx)
plot(SNR,real(Capacity_MIMO_norm1),'g:')
hold on
% Ploting MIMO capacity with Water Filling Algorithm(4Tx - 4Rx)
plot(SNR,Capacity_MIMO_WFA4,'m-')
hold on
% Ploting MIMO capacity with Water Filling Algorithm(3Tx - 3Rx)
plot(SNR,Capacity_MIMO_WFA3,'m--')
hold on
% Ploting MIMO capacity with Water Filling Algorithm(2Tx - 2Rx)
plot(SNR,Capacity_MIMO_WFA2,'m-.')
hold on
% Ploting MIMO capacity with Water Filling Algorithm(1Tx - 1Rx)
plot(SNR,Capacity_MIMO_WFA1,'m:')
xlabel('SNR(dB)')
ylabel('Capacity')
title('Capacity/SNR(dB)')
legend('SIMO','MISO with CSI','MISO without CSI','MIMO with EQD(4X4)','MIMO with EQD(3X3)','MIMO with EQD(2X2)','MIMO with EQD(1X1)','MIMO with WFA(4X4)','MIMO with WFA(3X3)','MIMO with WFA(2X2)','MIMO with WFA(1X1)','Location','northwest')