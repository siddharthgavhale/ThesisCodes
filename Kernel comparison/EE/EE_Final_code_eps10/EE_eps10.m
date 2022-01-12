
%EE positive L^1 and L^inf error for fix epsilon=0.1
%IC--Wulff shape : Ellipse, %Two phase : Shrinking Wulff shap
%Error table: Convergence table and few figures:

tic;
clc,close all;
fileID = fopen('Error and convergence table.txt','w');
fprintf(fileID,'%13s  %11s %8s %12s %13s\n', 'final time', 'diff_time','dx','l^1','l^inf');
fprintf(fileID, '\n');

gpiSV=7; gpiEV=13; 
dtiSV=0; dtiEV=7;
eps=0.1;
eps_save=eps*100;

%To export for plotting
L_1=zeros(gpiEV-gpiSV+1,dtiEV-dtiSV+1); L_inf=zeros(gpiEV-gpiSV+1,dtiEV-dtiSV+1);
DX=zeros(1,gpiEV-gpiSV+1); DT=zeros(1,dtiEV-dtiSV+1);

for gpi=gpiSV:gpiEV %Grid Point ind
    
for dti=dtiSV:dtiEV %Diffusion Time ind
    
diff_time=1/(4*2^dti);
    
% Initial setup for Numerical result    
N=2^gpi; L=10; dx = L/N;
x = -L/2:dx:L/2-dx; y = x; [X,Y]= meshgrid(x,y);

ax=1;ay=2;
time_steps=2/diff_time; %final time Decide time steps.

fprintf('space = %d\t time = %d \t BMOsteps = %d \n', gpi, dti,time_steps); 

%Initialization of Wulff shape
%sign distance function
U = 1.5 - (((X./ax).^2 + (Y./ay).^2).^(0.5));

% %Import saved kernel
%  load(sprintf('%s_%d %d %d','eps_gpi_dti',eps_save,gpi,dti));
 
 %  %Kernel %Elsey Esedoglu Kernel
K=zeros(N,N); 
for i=1:N
    for j=1:N  
fun=@(t) ((ax^2 * ay^2)./((ax^2.*cos(t+pi/2).^2 + ay^2.*sin(t+pi/2).^2).^(3/2)))...
    .*(exp(-((X(i,j).*cos(t) + Y(i,j).*sin(t)).^2)/(4*diff_time)))...
    .*(exp(-((X(i,j).*sin(t) - Y(i,j).*cos(t)).^2)/(4*diff_time*eps.^2)));
       K(i,j)=integral(fun,0,2*pi); 
    end
end
 K=(1/(16*eps*sqrt(pi))).*K;
K=(K./diff_time);


%Based on the numerical area of kernel for eps=0.05
K_area = 4.292184;
 
%Error terms, l^inf and l^1 
er_inf=0; er_1=0; er=zeros(1,fix(time_steps));

%Numerical solution
for BMOstep=1:time_steps

%Convolution     
U = (dx^2).*ifftshift(ifft2(fft2(fftshift(K)).*(fft2(fftshift(U)))));

%Reinitialization and smoothening of U for next BMOstep
U= (U-(K_area/2))./(K_area); %values between -0.5 and 0.5, interface is at 0
UC = U; % a copy of U is needed

for i=2:N-1
    for j=2:N-1
        if sign(U(i,j)) == sign(U(i,j+1)) == sign(U(i+1,j)) == sign(U(i,j-1)) ... 
                == sign(U(i-1,j)) == sign(U(i-1,j-1)) == sign(U(i-1,j+1)) ...
                == sign(U(i+1,j-1)) == sign(U(i+1,j+1)) 
            U(i,j)= 0.5*sign(U(i,j));
        else % compute the area ratio of positive phase on 4 grid cells centered at (i,j)
                b=10^-7;   
w1=UC(i-1,j-1); w2=UC(i-1,j); w3=UC(i,j-1);w4=UC(i,j);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=areacell;

w1=UC(i-1,j); w2=UC(i-1,j+1); w3=UC(i,j);w4=UC(i,j+1);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=Areacell+areacell;

w1=UC(i,j-1); w2=UC(i,j); w3=UC(i+1,j-1);w4=UC(i+1,j);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=Areacell+areacell;


w1=UC(i,j); w2=UC(i,j+1); w3=UC(i+1,j);w4=UC(i+1,j+1);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=Areacell+areacell;            
            U(i,j) = Areacell/4.0 -0.5;
        end
    end
end
U=U+0.5;

%Adjusting boundary
U(1,:) = 0; U(N,:) = 0; 
U(:,1) = 0; U(:,N) = 0;

%Analytical solution(Front tracking soliution)
TT= BMOstep*diff_time;%time

% At qurter time(approx) to extinction, for approximate L^1 error
if TT == 1 
    %sign distance function using front tracking solution
    A1=SignDistFunc(gpi,1);
    
    %subgrid improvement
    A1=smooth(A1,N);

    %l^1 error 
    er_1=sum(abs(U-A1),'all').*dx^2;
end

% At Half time (approx) to extinction, to approximate L^inf error
if TT == 2 
    %sign distance function using front tracking solution
    A2=SignDistFunc(gpi,2);
    
    %subgrid improvement
    A2=smooth(A2,N);
    
    %l^inf error 
    er_inf=sum(abs(U-A2),'all').*dx^2;
end

     
end%end of BMOsteps

%save final solution (both front tracking and kernel)
save(sprintf('%s_%d %d %d','Ueps_gpi_dti',eps_save,gpi,dti),'U');
save(sprintf('%s_%d %d %d','Aeps_gpi_dti',eps_save,gpi,dti),'A2');

%Export data to text file
Export_data=[(BMOstep)*diff_time diff_time dx er_1 er_inf];
fprintf(fileID,'%12.6f %12.6f %12.8f %12.8f %12.8f\n',Export_data);

%to export for plotting
L_1(gpi-gpiSV+1,dti-dtiSV+1) = er_1;
L_inf(gpi-gpiSV+1,dti-dtiSV+1) = er_inf;

if gpi == gpiSV
    DT(dti-dtiSV+1)=diff_time;
end
    

end

%to plot DX value
DX(gpi-gpiSV+1)=dx;

fprintf(fileID, '\n');
end


%Savind data for future plots
save('L_1.mat'); save('L_inf.mat');
save('DX.mat'); save('DT.mat');
%__________________


%Convergence table for fix dx (last dx): 
conv_dx_er_1=zeros(1,size(L_1,2)-1);
conv_dx_er_inf=zeros(1,size(L_1,2)-1);
for i=2:size(L_1,2)
    %With. r. to l^1
    conv_dx_er_1(i-1)= log(L_1(gpiEV-gpiSV+1,i-1)/L_1(gpiEV-gpiSV+1,i))/log(2);
    %With r. to l^inf
    conv_dx_er_inf(i-1)= log(L_inf(gpiEV-gpiSV+1,i-1)/L_inf(gpiEV-gpiSV+1,i))/log(2);
end
save('conv_dx_er_inf.mat'); save('conv_dx_er_1.mat');

fprintf(fileID,'%10s','------------------------------------------------------------------------');
 fprintf(fileID,'\n \n \n \n %12s','value of fix dx is');
fprintf(fileID,'\n %12.6f',dx);
 fprintf(fileID,'\n \n \n %12s \n','Convergence rate of l^1 error for fix dx');
 fprintf(fileID,'%10s','-');
 fprintf(fileID,'\n %12.6f',conv_dx_er_1);

 fprintf(fileID,'\n \n \n \n %12s  \n','Convergence rate of l^inf error for fix dx ');
 fprintf(fileID,'%10s','-');
 fprintf(fileID,'\n %12.6f',conv_dx_er_inf);
 
 
 %Convergence table for fix dt (last dt): 
 conv_dt_er_1=zeros(1,size(L_1,1)-1);
 conv_dt_er_inf=zeros(1,size(L_1,1)-1);
for i=2:size(L_1,1)
    %With. r. to l^1
    conv_dt_er_1(i-1)= (log(L_1(i-1,dtiEV-dtiSV+1)/L_1(i,dtiEV-dtiSV+1)))/log(2);
    %With r. to l^inf
    conv_dt_er_inf(i-1)= log(L_inf(i-1,dtiEV-dtiSV+1)/L_inf(i,dtiEV-dtiSV+1))/log(2);
end
save('conv_dt_er_inf.mat'); save('conv_dt_er_1.mat');
fprintf(fileID,'\n \n \n \n  %10s','------------------------------------------------------------------------');
 fprintf(fileID,'\n \n \n \n %12s','value of fix dt is');
fprintf(fileID,'\n %12.6f',diff_time);

 fprintf(fileID,'\n \n \n %12s \n','Convergence rate of l^1 error for fix dt');
 fprintf(fileID,'%10s','-');
 fprintf(fileID,'\n %12.6f',conv_dt_er_1);

 fprintf(fileID,'\n \n \n \n %12s \n','Convergence rate of l^inf error for fix dx');
 fprintf(fileID,'%10s','-');
 fprintf(fileID,'\n %12.6f',conv_dt_er_inf);
 fprintf(fileID,'\n  %10s','------------------------------------------------------------------------');

 %______________________
 fclose(fileID);
 
 disp('Code done! Plotting start. Check Error and Convergence in txt file');
 



%Figure1 loglog plot of ( $l^1$error,diff time) for various dx 
figure()
for i=1:length(DX)
    dx=DX(i);
loglog(DT,L_1(i,:),'DisplayName', num2str(dx),'Linewidth',2),
title('loglog plot of l^1 vs diff time for various dx'), hold on
xlabel('diff. time'),ylabel('l^1'), lgd = legend;
lgd.FontSize = 14; 
grid on,  hold on
end

%Figure2 loglog plot of  ( $l^1$error,dx) for various diff time 
figure()
for i=1:length(DT)
    dt=DT(i);
loglog(DX,L_1(:,i),'DisplayName', num2str(dt),'Linewidth',2),
title('loglog plot of l^1 vs dx for various diff time'), 
ax = gca; ax.FontSize = 18;
hold on
xlabel('dx'),ylabel('l^1'), lgd = legend;
lgd.FontSize = 14; 
grid on,  hold on
end

%Figure 3
figure()
[DTT, DXX]= meshgrid(DT,DX);
surf(DXX,DTT,L_1), hold on,
xlabel('dx'), ylabel('diff time'),ax = gca; ax.FontSize = 18;
title('(L^1) error for different pair of dx and diff time')
view(2), colorbar

%same plots above mentioned for l^inf

%Figure 4 loglog plot of  ( $l^inf$error,diff time) for various dx 
figure()
for i=1:length(DX)
    dx=DX(i);
loglog(DT,L_inf(i,:),'DisplayName', num2str(dx),'Linewidth',2),
title('loglog plot of l^{inf} vs diff time for various dx'), 
ax = gca; ax.FontSize = 18; hold on
xlabel('diff. time'),ylabel('l^{inf}'), lgd = legend;
lgd.FontSize = 18; 
grid on,  hold on
end

%Figure 5 loglog plot of  ( $l^{inf}$error,dx) for various diff time 
figure()
for i=1:length(DT)
    dt=DT(i);
loglog(DX,L_inf(:,i),'DisplayName', num2str(dt),'Linewidth',2),
title('loglog plot of l^{inf} vs dx for various diff time'), 
ax = gca; ax.FontSize = 18; hold on
xlabel('dx'),ylabel('l^{inf}'), lgd = legend;
lgd.FontSize = 18; 
grid on,  hold on
end

%Figure 6
figure()
[DTT, DXX]= meshgrid(DT,DX);
surf(DXX,DTT,L_inf), hold on,
xlabel('dx'), ylabel('diff time'),
title('(L^{inf}) error for different pair of dx and diff time')
ax = gca; ax.FontSize = 18; hold on
view(2), colorbar, FontSize=18;

%Figure 7
figure()
%for fix dx l1 error
pDT=DT; pDT(1)=[];
plot(pDT,conv_dx_er_1,'b','DisplayName', num2str(dx),'Linewidth',2), 
lgd = legend; lgd.FontSize = 18; 
hold on
%for fix dt
pDX=DX; pDX(1)=[];
plot(pDX,conv_dt_er_1,'r','DisplayName', num2str(diff_time),'Linewidth',2), 
lgd = legend; lgd.FontSize = 18; 
%for fix dx l^inf error
plot(pDX,conv_dt_er_inf,'--mo','color','r','DisplayName', num2str(diff_time),'Linewidth',2),
plot(pDT,conv_dx_er_inf,'--mo','color','b','DisplayName', num2str(dx),'Linewidth',2),
hold off
title('Convergance rate Vs dx (red) for fix dt and Convergance rate Vs dt (blue) for fix dx. l^1 error: Solid line, l^inf error: dotted line'), 
ax = gca; ax.FontSize = 18; hold on; grid 
xlabel('log(dx) or log(dt)'),ylabel('Convergance rate'), hold off

%Figure 8 most important
figure()
th = 0:pi/50:2*pi;
px=zeros(1,length(DX));
py=zeros(1,length(DX));
rad=zeros(1,length(DX));
for i=1:length(DX)
    
[mint,idx] = min(L_1(i,:));
    
px(i)=DX(i); py(i)=DT(idx); rad(i)=mint;
    
%px(i)=log(DX(i)); py(i)=log(DT(idx)); rad(i)=log(minx);

x_c = rad(i) * cos(th) + px(i); y_c = rad(i) * sin(th) + py(i);

   plot(x_c,y_c,'DisplayName', num2str([px(i) py(i) rad(i)]),'Linewidth',1.5),
    lgd = legend; lgd.FontSize = 14; hold on
    
end
xlabel('dx'),ylabel('Diff_time'), axis equal 
plot(px,py,'--*'),
title('Optimal pairs (dx,dt). For given dx we found dt such that, (dx,dt) pair gives least error. radius of circle represent error at (dx,dt) center'), 
ax = gca; ax.FontSize = 18; hold on; grid 
hold off


load train
sound(y,Fs)
toc



