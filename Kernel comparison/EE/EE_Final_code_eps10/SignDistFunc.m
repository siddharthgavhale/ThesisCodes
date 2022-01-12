

%sign distance funciton for ellipse 

% function sdf=SignDistFunc(gpi,R)

clear;close;clc
 gpi=8; R=1;
N=2^gpi; L=10; dx = L/N;
x = -L/2:dx:L/2-dx; y = x; [X, Y]= meshgrid(x,y);

N1=2^8;

%load fron tracking solution
load(sprintf('%s_%d ','X_F',N1),'X_F');
load(sprintf('%s_%d ','Y_F',N1),'Y_F');

%extracing for qureter time or half time to extinction
F_X=X_F(R,:); F_Y=Y_F(R,:);

%erase time entry 
F_X(:,1)=[];F_Y(:,1)=[];

plot(F_X,F_Y), axis equal
kjsdh

sdf=zeros(N,N);

%first way
for i=1:N
    for j=1:N
         dis=999;
         %from each segment of front tracking solution find minimum distance
            for k=1:length(F_X)-1 
                di1 = sqrt((X(i,j)-F_X(k))^2+(Y(i,j)-F_Y(k))^2);
                di2 = sqrt((X(i,j)-F_X(k+1))^2+(Y(i,j)-F_Y(k+1))^2);

                t = -((F_X(k+1)-F_X(k))*(F_X(k)-X(i,j))+(F_Y(k+1)-F_Y(k))*(F_Y(k)-Y(i,j)))/((F_X(k+1)-F_X(k))^2+(F_Y(k+1)-F_Y(k))^2);
                if ( t>0 && t<1 )
                      di3 = abs((F_Y(k+1)-F_Y(k))*X(i,j)-(F_X(k+1)-F_X(k))*Y(i,j)+F_X(k+1)*F_Y(k)-F_Y(k+1)*F_X(k))/sqrt((F_X(k)-F_X(k+1))^2+(F_Y(k)-F_Y(k+1))^2);
                else
                       di3 = 1000;
                end

                    d = min([di1,di2,di3]);             
                    
                  if dis > d 
                     dis=d;
                  end
            end
    % check in side or outside
    [in,on] = inpolygon(X(i,j),Y(i,j),F_X,F_Y);
                    if on == 1 || in == 1 
                        if on == 1
                            sdf = 0;
                        else
                            sdf(i,j)= dis;
                        end
                    else
                        sdf(i,j)= -dis;
                    end
    end
end
% end
