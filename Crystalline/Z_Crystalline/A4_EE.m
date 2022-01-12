
%EE Error table: Convergence table:
%IC-- circle : crystalline, %Two phase : area preservation
%evolve first with 0-1 form, when particle is close to analytical solution
%the activate smoothening


function [gpi,dti] = EE(gpi,dti,epsi,plotting)

if plotting==1
vidObj = VideoWriter('Crystaline.avi'); 
open(vidObj)
end
       
% Initial setup for Numerical result    
N=2^gpi; L=10; dx = L/N;
x = -L/2:dx:L/2-dx; y = x; [X, Y]= meshgrid(x,y);

diff_time=1/(4*2^dti);%diffusion time

%Initial circle: in sign distance form at height 0.5
%it has area=4
U = 1.628379167095513 - (((X).^2 + (Y).^2).^(0.5));

%initial area, round off 
area_in=4/dx^2; area_in=round(area_in);

%Limit to terminate BMOstep loop
%if two consecutive shape is less than this limit then terminate
termi_limit = area_in*0.0001;

%Smoothening (subgrid improvement): when particle enough close to 
%analytical solution then use subgrid improvment
smooth_limit = termi_limit*10;

fprintf('space = %d\t time = %d \t \n', gpi, dti); 

%Kernel %Elsey Esedoglu Kernel
K=zeros(N,N); 
for i=1:N
    for j=1:N  
        g1=(exp(-(X(i,j)^2)/(4*diff_time)))*(exp((-Y(i,j).^2)/(4*diff_time*epsi.^2)));
        g2=(exp(-(Y(i,j)^2)/(4*diff_time)))*(exp((-X(i,j).^2)/(4*diff_time*epsi.^2)));
        K(i,j)=g1+g2;
    end
end
K=(1/(4*epsi*sqrt(pi))).*K;
K=(K./diff_time);
%Initial set up end________________________

%analytical solution in sign distance form 
A = 2 - abs(X+Y) - abs(X-Y);
%analytical solution in smooth 0-1 form
A=B_smooth(A,N);


%Numerical solution
movement=1000; count=0;

 while movement > termi_limit
    
count=count+1; %counting of BMOStep
 
U_prev=U;%previous time step matrix

%Convolution     
U = (dx^2).*ifftshift(ifft2(fft2(fftshift(K)).*(fft2(fftshift(U)))));

%with are preservation 
M=-U;
sortlist=sort(M(:)); 
delta=sortlist(area_in); %thresholding height

for i=1:N
    for j=1:N
        if M(i,j)<= delta
            U(i,j)=1;
        else
            U(i,j)=0;
        end
    end
end

%difference between previous and this timestep
movement= sum(abs(U_prev-U),'all');

%if 'movement' is less than limit then do smoothening
if smooth_limit > movement
U=B_smooth(U-0.5,N);
end


if plotting==1 
 contour(X,Y,A,[0.5 0.5],'b','linewidth',2), hold on
contour(X,Y,U,[0.5 0.5],'r','linewidth',2),
axis equal, grid on, grid minor
axis([-3 3 -3 3])
text(-2,-2,[' Time = ', num2str(count*diff_time)],'FontSize',20);
text(-2,2.4,[' gpi = ', num2str(gpi)],'FontSize',20);
text(-2,2,[' dti = ', num2str(dti)],'FontSize',20);
text(1,2.2,[' \epsilon = ', num2str(epsi)],'FontSize',20);
title('EE: Numerical:Red and Analytical:Blue')
 currFrame = getframe(gcf);  writeVideo(vidObj,currFrame); hold off
end


     
end%end of BMOsteps


%%save final solution (final matrix)
save(sprintf('%s_%d %d %d','EE_eps_gpi_dti',epsi*100,gpi,dti),'U');

if plotting==1 
  close(vidObj);
end

end



