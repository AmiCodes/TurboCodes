%This function calculates ALPHA probabilities at each stage for all states,
%using GAMMA probabilities obtained previously. 

function [ALPHA]=alpha_1(GAMMA,N)

    ALPHA=zeros(4,N);

    %Initialization of alpha assuming first state to be 00
    ALPHA(1,1)=1;ALPHA(2,1)=0;ALPHA(3,1)=0;ALPHA(4,1)=0;
    
    j=1;
    for i=2:N
        ALPHA(1,i)=((GAMMA(1,j)*ALPHA(1,i-1))+(GAMMA(3,j+1)*ALPHA(3,i-1)));
        ALPHA(2,i)=((GAMMA(3,j)*ALPHA(3,i-1))+(GAMMA(1,j+1)*ALPHA(1,i-1)));
        ALPHA(3,i)=((GAMMA(2,j)*ALPHA(2,i-1))+(GAMMA(4,j+1)*ALPHA(4,i-1)));
        ALPHA(4,i)=((GAMMA(4,j)*ALPHA(4,i-1))+(GAMMA(2,j+1)*ALPHA(2,i-1)));
        j=j+2;
        
        if (ALPHA(1,i)<10^(-20) && ALPHA(2,i)<10^(-20) &&...
                ALPHA(3,i)<10^(-20) && ALPHA(4,i)<10^(-20) )
            ALPHA(:,i)=10^(20)*ALPHA(:,i);   %Scaling Alpha if became very less  
        end    
    end

end