%function: Find angle between Normal and Y axis  
% since normal is not pointing outward direction I have made some 
% adjustments depending on the quadrant of normal.

function [th]=FD_find_angle(nd)

    N=length(nd);
    n1 = nd(1,1:N);    n2=nd(2,1:N); 
      
%increase the length of the normal to check its direction
    nn1=1000.*n1;     nn2=1000.*n2;
    
    th = zeros(1,N);
    for i=1:N        
        th(i)=atan2(abs(n2(i)),abs(n1(i)));
        
        %first Quadrant
        if nn1(i)>0 && nn2(i)> 0
            th(i)=pi/2+th(i); % adjustment 
        end
        %second Quadrant
        if nn1(i)<0 && nn2(i)> 0
             th(i)=-th(i)-pi/2; % adjustment 
        end
        %third Quadrant
        if nn1(i)<0 && nn2(i)< 0
            th(i)=th(i)-pi/2; % adjustment
        end
        %fourth Quadrant
        if nn1(i)>0 && nn2(i)< 0
            th(i)=pi/2-th(i); % adjustment 
        end
        
        if nn1(i)==0 || nn2(i)==0
            if nn1(i)==0 && nn2(i)>0
                th(i)=-pi;
            end
            if nn1(i)==0 && nn2(i)<0
                th(i)=0;
            end   
            
            if nn2(i)==0 && nn1(i)>0
                th(i)=pi/2;
            end
            if nn2(i)==0 && nn1(i)<0
                th(i)=-pi/2;
            end
        end     
    end
end