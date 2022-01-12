

%particle on patterned (two) substrate code
%Characteristic type:
clc; clear; close all;

gpiSV=12; gpiEV=12; %range of grid points %gpiu=12

dtiSV=6; 
dtiEV=dtiSV; %range of diffusion time

L=10; %L*L domain

S1_st=0;%:-0.1:-2%"gammaSP-gammaSV" value for subsrtate 1. % Substrate-1_surface tensions
S2_st=2;%:-0.1:-2;%"gammaSP-gammaSV" value for subsrtate 2. % Substrate-2_surface tensions

final_time = 80;

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

          
  %Initialization: particle(P), substrate(S) & vapor(V) & calculate
  %initial area of particle (analytical,in grid scale)
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
    
     %particle: 
         width=0.35; len=2.5;
         for i=1:N
            for j=1:N
                if X(i,j)> -len && X(i,j) < len && Y(i,j) < width && Y(i,j) >= 0
                    P(i,j) = 1; area=area+1;
                end
            end
         end    

        V=ones(N,N)-(P+S);
        P_initial=P; 
        V_initial=V; 
  
      
    %Until which row, we have substrate? Since half domain fixed for S
    sub_level = N/2 + 1;
    
         for dti=dtiSV:dtiEV %Diffusion Time ind

            fprintf('space = %d\t time = %d\n', gpi, dti);     

            %Initial setup of time related para. for Numerical result  
            diff_time=1/(4*2^dti); %diffusion time
            Ini_dt = diff_time; % For exporting data (in modified case, final diffusion time will be different)
            time_steps=final_time/diff_time; %maximum time steps.
            movement=movement_ini;

        %Initilize P and V
        P = P_initial;
        V = V_initial;

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

        first_entry=1; delta=0; temp8=1; ph2=zeros(N,N);
        %Main loop of thresholding
           while movement > 0 && final_time > TT && first_entry==1

               prev_M=ph1+ph2; prev_delta=delta;

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
               iteration=iteration+1; TT= TT+diff_time;

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

             if first_entry==1
                 if  sum(P(N/2:N,N/2)) == 0
                     first_entry=0;
                 end
             end
           end %%end of "while movement > 0"

    %first plot before topological change
    contour(X,Y,prev_M,[prev_delta,prev_delta],':b','LineWidth',3); hold on
    b1=iteration-1;

    %area of left particle and area of right particle

    P_L=zeros(N,N);
    P_R=zeros(N,N);

    P_L(1:N,1:N/2-1)=P(1:N,1:N/2-1);
    P_R(1:N,N/2+1:N)=P(1:N,N/2+1:N);


    area_L=sum(P_L,'all');
    area_R=sum(P_R,'all');

    P=P_L+P_R;


           %two particle evolution
            while movement > 0 && final_time > TT


               %convolution either with Kernel in Fourier form of regular form
                    ph2_L = ifftshift(ifft2((fftshift(K))...
                         .*(fft2(fftshift(V-P_L)))))./sqrt(diff_time);
                     ph2_L = real(ph2_L); 

                     ph2_R = ifftshift(ifft2((fftshift(K))...
                         .*(fft2(fftshift(V-P_R)))))./sqrt(diff_time);
                     ph2_R = real(ph2_R);


               M_R=ph1+ph2_R;
               M_L=ph1+ph2_L;

               %exclude all the values (from calculation of thresholding height)
               %which intersecting with substrate
               M_L(1:sub_level-1,:) = 99; 
               M_R(1:sub_level-1,:) = 99; 


               %right particle
               M_R(1:N,1:N/2-1)= 99;
               M_L(1:N,N/2+1:N)= 99;


               sortlist_L=sort(M_L(:)); %sorting
               sortlist_R=sort(M_R(:)); %sorting

               delta_L=sortlist_L(area_L); %thresholding height, it preserves area.
               delta_R=sortlist_R(area_R); %thresholding height, it preserves area.

               %iteratin  details and total time
               iteration=iteration+1; TT= TT+diff_time;

                       if iteration == b1+3
                       contour(X,Y,M_L,[delta_L,delta_L],'--m','LineWidth',3); 
                       contour(X,Y,M_R,[delta_R,delta_R],'--m','LineWidth',3);
                       end
                       if iteration == 70
                       contour(X,Y,M_L,[delta_L,delta_L],'-.k','LineWidth',3);
                       contour(X,Y,M_R,[delta_R,delta_R],'-.k','LineWidth',3);
                       end
                       axis equal, 
                       axis([-2.5 2.5 0 2.5]) 


               %reinitialization
               P_new_L=zeros(N,N);
               P_new_R=zeros(N,N);
                 for i=1:N
                    for j=1:N/2
                        if M_L(i,j) <= delta_L
                            P_new_L(i,j)=1;
                        end
                    end
                    for j=N/2:N
                        if M_R(i,j) <= delta_R
                            P_new_R(i,j)=1;
                        end
                    end
                 end
              V=ones(N,N)-P_new_L-P_new_R-S;

              %Change between, two diffusion time
              movement_L= sum(abs(P_L-P_new_L),'all');   
              movement_R= sum(abs(P_R-P_new_R),'all');
              if movement_L <= movement_R
                  movement=movement_R;
              else
                  movement = movement_L;
              end


              P_L=P_new_L;  P_R=P_new_R; 

           end %%end of "while movement > 0"


           %reset to zero
           iteration=0; TT= 0;

        end

end

    %intial particle
    contour(X,Y,P_initial,[0.5,0.5],'-k','LineWidth',3);  
    %final shape
    contour(X,Y,M_L,[delta_L,delta_L],'-k','LineWidth',3); hold on
    contour(X,Y,M_R,[delta_R,delta_R],'-k','LineWidth',3);


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
     axis([-2.75 2.75 -1 2.5]),
     grid off

%maximise the plot and save
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
saveas(gcf,'Fig8_1.png')
hold off

 







