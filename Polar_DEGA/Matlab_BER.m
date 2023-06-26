load Polar_AWAN_QPSK_A01_GIR001.txt
i=Polar_AWAN_QPSK_A01_GIR001(:,1);
j=Polar_AWAN_QPSK_A01_GIR001(:,2);

I = semilogy(i,j,'r-O');
set(I,'linewidth',2,'Markersize',8);
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_1024512_QPSK_AWAN_A001GIR001_DEGA.txt
% k=SPolar_1024512_QPSK_AWAN_A001GIR001_DEGA(:,1);
% l=SPolar_1024512_QPSK_AWAN_A001GIR001_DEGA(:,2);
% 
% K = semilogy(k,l,'r-*');
% set(K,'linewidth',2,'Markersize',8);
% hold on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_1024512_QPSK_AWAN_A001GIR001_knownimpulse.txt
% m=SPolar_1024512_QPSK_AWAN_A001GIR001_knownimpulse(:,1);
% n=SPolar_1024512_QPSK_AWAN_A001GIR001_knownimpulse(:,2);
% 
% M = semilogy(m,n,'b-^');
% set(M,'linewidth',2,'Markersize',8);
% hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_AWAN_BPSK_A01_GIR01_LT1.txt
% o=SPolar_AWAN_BPSK_A01_GIR01_LT1(:,1);
% p=SPolar_AWAN_BPSK_A01_GIR01_LT1(:,2);
% 
% O = semilogy(o,p,'b-*');
% set(O,'linewidth',2,'Markersize',8);
% hold on;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Polar_AWAN_QPSK_A01_GIR001.txt
% q=Polar_AWAN_QPSK_A01_GIR001(:,1);
% r=Polar_AWAN_QPSK_A01_GIR001(:,2);
% 
% Q = semilogy(q,r,'g-o');
% set(Q,'linewidth',2,'Markersize',8);
% hold on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Polar_AWAN_BPSK_A01_GIR001.txt
% s=Polar_AWAN_BPSK_A01_GIR001(:,1);
% t=Polar_AWAN_BPSK_A01_GIR001(:,2);
% 
% S = semilogy(s,t,'g-*');
% set(S,'linewidth',2,'Markersize',8);
% hold on;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_20481723_05test.txt
% q=SPolar_20481723_05test(:,1);
% r=SPolar_20481723_05test(:,2);
% 
% Q = semilogy(q,r,'m-o');
% set(Q,'linewidth',2,'Markersize',8);
% hold on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_20481723_BG05_CF.txt
% s=SPolar_20481723_BG05_CF(:,1);
% t=SPolar_20481723_BG05_CF(:,2);
% 
% S = semilogy(s,t,'m-*');
% set(S,'linewidth',2,'Markersize',8);
% hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_256128_test1.txt
% s=SPolar_256128_test1(:,1);
% t=SPolar_256128_test1(:,2);
% 
% S = semilogy(s,t,'k-o');
% set(S,'linewidth',2,'Markersize',8);
% hold on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_256128_test_2.5.txt
% s=SPolar_256128_test(:,1);
% t=SPolar_256128_test(:,2);
% 
% S = semilogy(s,t,'m-o');
% set(S,'linewidth',2,'Markersize',8);
% hold on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load SPolar_256128_test_4.txt
% s=SPolar_256128_test_4(:,1);
% t=SPolar_256128_test_4(:,2);
% 
% S = semilogy(s,t,'b-o');
% set(S,'linewidth',2,'Markersize',8);
% hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eb_No_dB =0:0.5:40;
% EbNo_lin = 10.^(Eb_No_dB/10);
% %theoryBer_Rayl = 0.5.*(1-sqrt(EbNo_lin./(EbNo_lin+1)));
% theoryBer_AWGN = 0.5*erfc(sqrt(2*EbNo_lin/2));

set(gca,'FontWeight','bold','fontsize',17);
%set(gca,'xtick',[0:0.1:4]);% 增加軸刻度

ax = gca;
ax.GridAlpha = 0.5;
ax.LineWidth = 0.5;
ax.MinorGridAlpha = 1.0;

%semilogy(k,l,'r-o');
%semilogy(Eb_No_dB,theoryBer_AWGN,'k--O');
grid on;
axis([0 20 10^-6 1]);
%axis([1 7 10^-1 1]);


xlabel('Eb/N0 (dB)');% x軸註解
ylabel('BER');% y軸註解
%xlabel('樣本數');% x軸註解
%ylabel('雜訊大小');% y軸註解

title('SystematicPolarCode(1024,512) QPSK BG');
%title('SystematicPolarCode QPSK AWAN A0.1 GIR0.01');
%title('Performance of 0.5 (256, 128) Polar Code BPSK AWGN');
%title('Bernoulli-Gaussian')

%legend('QPSK A0.1 GIR0.1','BPSK A0.1 GIR0.1','QPSK A0.01 GIR0.1','BPSK A0.01 GIR0.1','QPSK A0.1 GIR0.01','BPSK A0.1 GIR0.01'); % 圖形註解
legend('Bhattacharyya','DEGA','frozenbits for impulse'); % 圖形註解
%legend('P_b=0.03, IGR=100, Gaussian variance=0.5');

