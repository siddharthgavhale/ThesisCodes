
%EJZ PD %IC-- S shape

tic; clc, clear, close all;

plotting=0; %want plotting? yes=1, No=0
Area_preserve=0;% Preserve area: Yes=1; No=0
Front_track_plt=0; %want front tracking solution in plot
final_time=1.251;

gpiSV=13; gpiEV=13; 
dtiSV=9; dtiEV=9;

%anisotropy details
bet=.05; mf=4;
K_area=2.04333944;

if plotting==1
vidObj = VideoWriter('S shape EJZ PD.avi');  open(vidObj); end


for gpi=gpiSV:gpiEV %Grid Point ind
    
    %Initial setup for Numerical result    
     N=2^gpi; L=10; dx = L/N;
     x = -L/2:dx:L/2-dx; y = x; [X, Y]= meshgrid(x,y);
        
    %Initialization of Wulff shape
    U=zeros(N,N); area_in=0;
     for i=1:N
         for j=1:N
              if  X(i,j)>= 0 %right side of domain
                  if Y(i,j) <= 0 %  down-right side of domain (4th quadrant)
                       if  (X(i,j))^2 + (Y(i,j)+2)^2 <= 4 && ...
                                (X(i,j))^2 + (Y(i,j)+2)^2 >= 1
                            U(i,j) =1; area_in=area_in+1;
                       end
                  else  %up-right side quarter circle
                       if  (X(i,j))^2 + (Y(i,j)-1)^2 <= 4 && ...
                               (X(i,j))^2 + (Y(i,j)-1)^2 >= 1 && Y(i,j) >= 1
                            U(i,j) =1; area_in=area_in+1;
                       end
                  end        
              else
                  if Y(i,j) >= -1 %  up-left side of domain (2nd quadrant)
                       if  (X(i,j))^2 + (Y(i,j)-1)^2 <= 4 && ...
                                (X(i,j))^2 + (Y(i,j)-1)^2 >= 1
                            U(i,j) =1; area_in=area_in+1;
                       end
                  else  %down-left side : quarter circle
                       if  (X(i,j))^2 + (Y(i,j)+2)^2 <= 4 && ...
                               (X(i,j))^2 + (Y(i,j)+2)^2 >= 1 && Y(i,j) <= -2
                            U(i,j) =1; area_in=area_in+1;
                       end
                  end
              end

         end
     end
     U_ini=U; % same initial shape for THIS gpi
    
    for dti=dtiSV:dtiEV %Diffusion Time ind
       
        diff_time=1/(4*2^dti);%diffusion time
        time_steps=final_time/diff_time; %final time Decide time steps
        fprintf('space = %d\t time = %d \t \n', gpi, dti); 
        TT=0;%total  time count
        if Front_track_plt == 1
            %front tracking solution
            load('X_F_102 .mat')
            load('Y_F_102 .mat')
            ftp=1; %fix value
        end
        
        %Kernel EJZ  PD
        m0 = 0.34029423827512594869;
        m2 = 0.37372385180378853153;
        K=zeros(N,N); 
        for i= 1:N
            for j=1:N
                th = atan2(Y(i,j),X(i,j)) - pi/2; % Polar form: angle 
                rt = sqrt((Y(i,j)^2 + X(i,j)^2)/diff_time); % Polar form

                ani= 1+ bet*cos(mf*th);% anisotropy
                stiff= 1+ bet*cos(mf*th) - bet*mf^2*cos(mf*th);% stiff=g+g"
                mu=ani;%mobility = anisotropy

                beta = sqrt(m2/(2*m0*mu*stiff)); 
                epw = rt*beta;

                if 0 < epw && epw < 2
                    eta = exp(-1/((epw^2)*(epw-2)^2));

                    alpha = beta/(4*m0*mu);

                    K(i,j)=alpha*eta;
                end

            end
        end
        K=(K./diff_time);
        
        U = U_ini; 
        movement=999;
        
        %Numerical solution
        for BMOstep=1:time_steps
 
        U_prev=U;
        
        %Convolution     
        U = (dx^2).*ifftshift(ifft2(fft2(fftshift(K)).*(fft2(fftshift(U)))));
        
        if Area_preserve ==1
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
            area=area_in;
        else
            for i=1:N
              for j=1:N 
                 if U(i,j) >= K_area/2
                    U(i,j) =1;
                 else
                    U(i,j)=0;
                 end
              end
            end           
               area = sum(U,'all'); %check area
        end

        movement= sum(abs(U_prev-U),'all');
        TT=BMOstep*diff_time; %total time
        
        if Front_track_plt==1 &&  plotting==1  
            X_F(ftp,1)

            while X_F(ftp,1) < TT
                 ftp=ftp+1;
            end
            XX_F=X_F(ftp,:); %copy perticular row 
            XX_F(1)=[];%exclude time entry
            YY_F=Y_F(ftp,:);
            YY_F(1)=[];%exclude time entry
        end
            
             
        if plotting==1
            if Front_track_plt==1 
                plot(XX_F,YY_F,'b','linewidth',2); hold on; end
            contour(X,Y,U,[0.5 0.5],'r','linewidth',2),  hold on
            axis equal, axis([-4 4 -4 4]), %grid on, grid minor;
            text(-3.5,-3.5,[' Time = ', num2str(BMOstep*diff_time)],'FontSize',20);
            text(-3.5,2.5,[' gpi = ', num2str(gpi)],'FontSize',20);
            text(-3.5,2,[' dti = ', num2str(dti)],'FontSize',20);
            title('Blue: Front tracking, Red: Numerical');
            currFrame = getframe(gcf);  writeVideo(vidObj,currFrame); hold off
         end
        
        %save file for plotting
        if TT==0.25
         save(sprintf('%s_%d %d %d','EJZPD_U_gpi_dti_TT',gpi,dti,TT*10000),'U'); 
        end
        if TT==0.75
         save(sprintf('%s_%d %d %d','EJZPD_U_gpi_dti_TT',gpi,dti,TT*10000),'U'); 
        end
        if TT==1.25
         save(sprintf('%s_%d %d %d','EJZPD_U_gpi_dti_TT',gpi,dti,TT*10000),'U'); 
         break
        end


        end%end of BMOsteps
        

    end

end


if plotting==1 
  close(vidObj);
end

load train
sound(y,Fs)
toc



