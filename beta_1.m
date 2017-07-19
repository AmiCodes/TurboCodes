%This function calculates BETA probabilities at each stage for all states,
%using GAMMA probabilities obtained previously. 

function [BETA]=beta_1(GAMMA,N)
    
    BETA=zeros(4,N);
    %Initialization assuming the final stage to be 00
    BETA(1,N)=1;BETA(2,N)=0;BETA(3,N)=0;BETA(4,N)=0;
    
    j=2*N-1;
    for i=N-1:-1:1
       BETA(1,i)=(GAMMA(1,j)*BETA(1,i+1))+(GAMMA(1,j+1)*BETA(2,i+1));
       BETA(2,i)=(GAMMA(2,j)*BETA(3,i+1))+(GAMMA(2,j+1)*BETA(4,i+1));
       BETA(3,i)=(GAMMA(3,j)*BETA(2,i+1))+(GAMMA(3,j+1)*BETA(1,i+1));
       BETA(4,i)=(GAMMA(4,j)*BETA(4,i+1))+(GAMMA(4,j+1)*BETA(3,i+1));
       j=j-2; 
       
       if (BETA(1,i)<10^(-20) && BETA(2,i)<10^(-20) &&...
               BETA(3,i)<10^(-20) && BETA(4,i)<10^(-20) )
           BETA(:,i)=10^(20)*BETA(:,i);         %Scaling beta if became very less      
       end
    end
    
    
end