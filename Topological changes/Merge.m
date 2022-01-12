
%megrer

%particle on patterned (two) substrate code
%Characteristic type:
clc; clear; close all;

gpiSV=12; gpiEV=12; %range of grid points gpi=12
dtiSV=6; 
dtiEV=dtiSV;

L=10; %L*L domain

S1_st=0;%:-0.1:-2%"gammaSP-gammaSV" value for subsrtate 1. % Substrate-1_surface tensions
S2_st=-0.1;%:-0.1:-2;%"gammaSP-gammaSV" value for subsrtate 2. % Substrate-2_surface tensions

final_time = 50;

%anisotropy choice
% m fold= 1, 
%coming soon......crystalline=2, ellipse = 3, ellipse without smooth end =4,
% if anisotropy is m fold then provide more info (uncomment next to next line)
anisotropy = 1;
beta=0.3; phi=pi/2; m=2;%m-fold anisotropy,

movement_ini=9999; %to track how much particle change within one time step
syms th %symbolic theta
TT=0; %total time
iteration=0; %iteration count


for gpi=gpiSV:gpiEV %Grid Point ind
    
   %Initial setup of  domain for Numerical result  
   N=2^gpi;%grid points in one direction      
   dx = L/N; x = -L/2:dx:L/2-dx; 
   y = x; [X, Y]= meshgrid(x,y);
   

     LF=1/dx; dxf = 1/(dx*N);
     xf= -LF/2:dxf:LF/2-dxf; yf = xf; [XF, YF] = meshgrid(xf,yf); 

          
  %Initialization: particle(P), substrate(S) & vapor(V) 
  P=zeros(N,N); S1=zeros(N,N); S2=zeros(N,N);

    %substrate
    S2_dist=0.5;
     for i=1:N
            for j=1:N
                 if Y(i,j) <= 0 %substrate region           
                         S1(i,j)=1;
                     if -S2_dist <= X(i,j) && X(i,j) <= S2_dist
                         S1(i,j)=0;
                         S2(i,j)=1;
                     end
                end
            end
     end
    S=S1+S2;

    area=0;
        for i=1:N
            for j=1:N
                if X(i,j) <=-0.2   &&  Y(i,j) <= X(i,j) + 2 && Y(i,j) >= 0
                    P(i,j) = 1; area=area+1;
                end
                 if X(i,j) >=0.2   &&  Y(i,j) <= tan(-pi/4)*X(i,j) + 2 && Y(i,j) >= 0
                    P(i,j) = 1; area=area+1;
                end
            end
        end

    V=ones(N,N)-(P+S);
    P_initial=P;
     
    %Until which row, we have substrate? Since half domain fixed for S
    sub_level = N/2 + 1;
    
    for dti=dtiSV:dtiEV %Diffusion Time ind

        fprintf('space = %d\t time = %d\n', gpi, dti);     

        %Initial setup of time related para. for Numerical result  
        diff_time=1/(4*2^dti); %diffusion time
        Ini_dt = diff_time; % For exporting data (in modified case, final diffusion time will be different)
        time_steps=final_time/diff_time; %maximum time steps.
        movement=movement_ini;
        
     
    %Since only one kernel option given therefore not writtting separate
    %kernel
    Gau=(1/(4*pi*(diff_time))).*exp(-(X.^2+Y.^2)./(4*(diff_time))); 
     K=zeros(N,N);
        for i=1:N
            for j=1:N  
                  thr = atan2(YF(i,j),XF(i,j));
                  g = (1 + beta*cos(m*(thr+phi)))*sqrt((XF(i,j))^2 + (YF(i,j))^2);
                  K(i,j)=exp(-4*pi*pi*diff_time*g^2);   
     
            end
        end  
        
   
        
    %first convolution with substrate interface 
    ph1_S1 = (dx^2).*ifftshift(ifft2(fft2(fftshift(Gau))...
              .*(fft2(fftshift((S1_st).*S1)))))./sqrt(diff_time);
    ph1_S2 = (dx^2).*ifftshift(ifft2(fft2(fftshift(Gau))...
              .*(fft2(fftshift((S2_st).*S2)))))./sqrt(diff_time);
          
    ph1 = ph1_S1+ ph1_S2;
              

    %Main loop of thresholding
       while movement > 0 && final_time > TT 
           
           
           %convolution either with Kernel in Fourier form of regular form
                ph2 = ifftshift(ifft2((fftshift(K))...
                     .*(fft2(fftshift(V-P)))))./sqrt(diff_time);
                 ph2 = real(ph2);   
    
    
           M=ph1+ph2;


           %exclude all the values (from calculation of thresholding height)
           %which intersecting with substrate
           M(1:sub_level,:) = 99999999999; 
           sortlist=sort(M(:)); %sorting
           delta=sortlist(area); %thresholding height, it preserves area.
           
           %iteratin  details and total time
           iteration=iteration+1;
           TT= TT+diff_time;
    
               if iteration==49    
                   contour(X,Y,ph1+ph2,[delta,delta],':b','LineWidth',3); hold on
               end
               if iteration==55   
                   contour(X,Y,ph1+ph2,[delta,delta],'--m','LineWidth',3); hold on
               end
               
              if iteration==75   
                   contour(X,Y,ph1+ph2,[delta,delta],'-.k','LineWidth',3); hold on
              end
           
           %reinitialization
           P_new=zeros(N,N);
             for i=1:N
                for j=1:N
                    if M(i,j) <= delta
                        P_new(i,j)=1;
                    end
                end
             end
          V=ones(N,N)-P_new-S;
          
          %Change between, two diffusion time
          movement= sum(abs(P-P_new),'all');     
         
          P=P_new;  
          
 
       end %%end of "while movement > 0"
 
      
    
    end


end

%intial particle
contour(X,Y,P_initial,[0.5,0.5],'-k','LineWidth',3);  
%final shape
contour(X,Y,ph1+ph2,[delta,delta],'-k','LineWidth',3); hold on


set(gca,'FontSize',40)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);

hold on
xstart=-4;
ystart=-1;
width=-(xstart+S2_dist);
height=-ystart;

rectangle('Position',[xstart,ystart,width,height],'FaceColor','w','LineWidth',2)
for i = 1:40
    x = xstart+i*width/41;
    line([x x],[ystart ystart+height],'Color','k');
end
text(-1.8,-0.5,' S1 ','FontSize',80);

xstart=S2_dist;
ystart=-1;
width=(4-S2_dist);
height=-ystart;

rectangle('Position',[xstart,ystart,width,height],'FaceColor','w','LineWidth',2)
for i = 1:40
    x = xstart+i*width/41;
    line([x x],[ystart ystart+height],'Color','k');
end
text(1.2,-0.5,' S1 ','FontSize',80);

xstart=-S2_dist;
ystart=-1;
width=(2*S2_dist);
height=-ystart;

rectangle('Position',[xstart,ystart,width,height],'FaceColor','w',...
    'EdgeColor','b')
for i = 1:10
    y = ystart+i*height/11;
    line([xstart xstart+width],[y y],'Color','k');
end
text(-0.3,-0.5,' S2 ','FontSize',80);
 axis([-2.5 2.5 -1 2.5]) 

 %maximise the plot and save
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
saveas(gcf,'Fig8_2.png')



 







