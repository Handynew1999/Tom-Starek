function [S1,I1,R1,S2,I2,R2,S3,I3,R3,R01,R02,R03,x1,x2,x3,y1,y2,y3,P1,P2,P3,P4,P5,P6,P7,P8,ti1,ti2,ti3,ti4,ti5,ti6,ti7,ti8,R021,R022] = sir_model_transfer_matrix_6(beta,gamma,N,I0,T,dt,delta)
    P=0;
    S1 = zeros(1,T/dt);
    S1(1) = (N(1)-I0(1))*(1-0.5);
    I1 = zeros(1,T/dt);
    I1(1) = I0(1);
    R1 = zeros(1,T/dt);
    R1(1) = (N(1)-I0(1))*(0.5);
    R01= zeros(1,T/dt);
    R021= zeros(1,T/dt);
    S2 = zeros(1,T/dt);
    S2(1) = (N(2)-I0(2))*(1-0.2);
    I2 = zeros(1,T/dt);
    I2(1) = I0(2);
    R2= zeros(1,T/dt);
    R2(1)=(N(2)-I0(2))*(0.2);
    R02= zeros(1,T/dt);
    R022= zeros(1,T/dt);
    S3= zeros(1,T/dt);
    S3(1) = (N(3)-I0(3))*(1-P);
    I3 = zeros(1,T/dt);
    I3(1) = I0(3);
    R3 = zeros(1,T/dt);
    R3(1)=(N(3)-I0(3))*(P);
    R03= zeros(1,T/dt);
    x1= zeros(1,T/dt);
    x1(1)=I1(1)/N(1);
    y1= zeros(1,T/dt);
    y1(1)=S1(1)/N(1);
    x2= zeros(1,T/dt);
    x2(1)=I2(1)/N(2);
    y2= zeros(1,T/dt);
    y2(1)=S2(1)/N(2);
    x3= zeros(1,T/dt);
    x3(1)=I3(1)/N(3);
    y3= zeros(1,T/dt);
    y3(1)=S3(1)/N(3);
    tau=2/3;
    alfa=0;
    opatrenie11=0;
    opatrenie12=0;
    opatrenie13=0;
    opatrenie14=0;
    opatrenie21=0;
    opatrenie22=0;
    opatrenie23=0;
    opatrenie24=0;
    pomm1=0;
    pomm2=0;
    pomm3=0;
    P1=0;
    P2=0;
    P3=0;
    P4=0;
    P5=0;
    P6=0;
    P7=0;
    P8=0;
    ti1=0;
    ti2=0;
    ti3=0;
    ti4=0;
    ti5=0;
    ti6=0;
    ti7=0;
    ti8=0;
    for tt = 1:(T/dt)-1
        %PRVA KRAJINA%
        dS1 =-tau*(S1(tt)*I1(tt)*beta(1))/N(1) + delta*R1(tt)  - alfa*(1-tau)*((S1(tt)-y1(tt)*(271.43))*(6800*x2(tt)*beta(2)+(I1(tt)-x1(tt)*(271.43))*beta(1)))/(N(1)-(271.43)+(6800))-alfa*(1-tau)*y1(tt)*((271.43*((I2(tt)-x2(tt)*(6800))*beta(2)+271.43*x1(tt)*beta(1)))/(N(2)-(6800)+271.43));
        dI1 =   tau*(S1(tt)*I1(tt)*beta(1))/N(1) - gamma*I1(tt) + alfa*(1-tau)*((S1(tt)-y1(tt)*(271.43))*(6800*x2(tt)*beta(2)+(I1(tt)-x1(tt)*(271.43))*beta(1)))/(N(1)-(271.43)+(6800))+alfa*(1-tau)*y1(tt)*((271.43*((I2(tt)-x2(tt)*(6800))*beta(2)+271.43*x1(tt)*beta(1)))/(N(2)-(6800)+271.43));      
        dR1 = (gamma*I1(tt)-delta*R1(tt)) * dt;
        
        S1(tt+1) = S1(tt) + dS1;
        I1(tt+1) = I1(tt) + dI1;
        R1(tt+1) = R1(tt) + dR1;
        
        %VYPOCITANIE X a Y%
        x1(tt)=I1(tt)/N(1);
        y1(tt)=S1(tt)/N(1);
            
       %OPATRENIA PRE PRVU KRAJINU (DRUHA VLNA VACSIA AKO PRVA)%
        if I1(tt)/(S1(tt)+R1(tt))>=0.05 && opatrenie11==0
         opatrenie11=1;
         beta(1)=beta(1)*0.86;
         P1=I1(tt);
         ti1=tt;
        end
        


        
         %VYPOCITANIE R0%
         if tt>2 & tt<numel(R01)
          %R02(tt)=(I2(tt+1)-I2(tt-1))/(2*gamma*(I2(tt)))+1; %posunuté o
          %pol dòa cca
          R01(tt)=(I1(tt+1)-I1(tt))/(gamma*(I1(tt)))+1; %najpresnejšie
            end
         R021(tt)=(1/gamma)*((S1(tt)*beta(1))/N(1))*(tau+(1-tau)*alfa); 
        %DRUHA KRAJINA%
        dS2 =  -tau*(S2(tt)*I2(tt)*beta(2))/N(2) + delta*R2(tt) - alfa*(1-tau)*((S2(tt)-y2(tt)*(6800))*(271.43*x1(tt)*beta(1)+(I2(tt)-x2(tt)*(6800))*beta(2)))/(N(2)-(6800)+(271.43))-alfa*(1-tau)*y2(tt)*((6800*((I1(tt)-x1(tt)*(271.43))*beta(1)+6800*x2(tt)*beta(2)))/(N(1)-(271.43)+(6800)));
        dI2 =  tau*(S2(tt)*I2(tt)*beta(2))/N(2) - gamma*I2(tt) + alfa*(1-tau)*((S2(tt)-y2(tt)*(6800))*(271.43*x1(tt)*beta(1)+(I2(tt)-x2(tt)*(6800))*beta(2)))/(N(2)-(6800)+(271.43))+alfa*(1-tau)*y2(tt)*((6800*((I1(tt)-x1(tt)*(271.43))*beta(1)+6800*x2(tt)*beta(2)))/(N(1)-(271.43)+(6800)));
        dR2 = (gamma*I2(tt) - delta*R2(tt)) * dt;
        
        S2(tt+1) = S2(tt) + dS2;
        I2(tt+1) = I2(tt) + dI2;
        R2(tt+1) = R2(tt) + dR2;
        %VYPOCITANIE X a Y%
        x2(tt)=I2(tt)/N(2);
        y2(tt)=S2(tt)/N(2);
       
        %OPATRENIA PRE DRUHU KRAJINU%
        if I2(tt)/(S2(tt)+R2(tt))>=0.05 && opatrenie21==0
         opatrenie21=1;
         beta(2)=beta(2)*0.75;
         P5=I2(tt);
         ti5=tt;
        end
        

        
         %VYPOCITANIE R0%
        if tt>2 & tt<numel(R02)
          %R02(tt)=(I2(tt+1)-I2(tt-1))/(2*gamma*(I2(tt)))+1; %posunuté o
          %pol dòa cca
          R02(tt)=(I2(tt+1)-I2(tt))/(gamma*(I2(tt)))+1; %najpresnejšie
            end
         R022(tt)=(1/gamma)*((S2(tt)*beta(2))/N(2))*(tau+(1-tau)*alfa);  

    end
      
    end
