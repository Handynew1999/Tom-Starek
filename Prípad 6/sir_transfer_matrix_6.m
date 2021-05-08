clear;
clc;
beta= 10^-1*[2;2;2];
gamma = 0.07; 
delta = 1/360; 
N = 10^6*[8.87;5.54;10.69];
I0=[1;1;0]; 
T = 1500; 
dt = 1; 
fprintf('Hodnota parametru R0 v AT je %.2f',beta(1)/gamma)
fprintf('Hodnota parametru R0 v SK je %.2f',beta(2)/gamma)
fprintf('Hodnota parametru R0 v CZ je %.2f',beta(3)/gamma)
[S1,I1,R1,S2,I2,R2,S3,I3,R3,R01,R02,R03,x1,x2,x3,y1,y2,y3,P1,P2,P3,P4,P5,P6,P7,P8,ti1,ti2,ti3,ti4,ti5,ti6,ti7,ti8,R021,R022] = sir_model_transfer_matrix_6(beta,gamma,N,I0,T,dt,delta);
tt = 0:dt:T-dt;
subplot(2,2,1);
plot(tt,S1/N(1),tt,I1/N(1),tt,R1/N(1),'LineWidth',2); grid on; hold on;
title('Rakúsko');
xlabel('Poèet dní'); ylabel('Poèet obyvate¾ov (%)');
legend('S1','I1','R1');
subplot(2,2,3);
plot(tt,R01,tt,R021,'LineWidth',2);  
title('Rt');
xlabel('Poèet dní'); ylabel('Hodnota Rt');
legend('Rt(log)','Rt(Numerický výpoèet)');
subplot(2,2,2);
plot(tt,S2/N(2),tt,I2/N(2),tt,R2/N(2),'LineWidth',2); grid on; hold on;
title('Slovensko');
xlabel('Poèet dní'); ylabel('Poèet obyvate¾ov (%)');
legend('S2','I2','R2');
subplot(2,2,4);
plot(tt,R02,tt,R022,'LineWidth',2); 
title('Rt');
xlabel('Poèet dní'); ylabel('Hodnota Rt');
legend('Rt(log)','Rt(Numerický výpoèet)');