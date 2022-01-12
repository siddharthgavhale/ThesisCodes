
%EJZ FD %IC-- S shape

plotting=0; %want plotting? yes=1, No=0
Area_preserve=0;% Preserve area: Yes=1; No=0
Front_track_plt=0; %want front tracking solution in plot
final_time=1.251;

gpiSV=13; gpiEV=13; 
dtiSV=9; dtiEV=9;

%anisotropy details
bet=.05; mf=4;

if plotting==1
vidObj = VideoWriter('S shape EJZ FD.avi');  open(vidObj); end

for gpi=gpiSV:gpiEV %Grid Point ind
    
    %Initial setup for Numerical result    
     N=2^gpi; L=10; dx = L/N;
     x = -L/2:dx:L/2-dx; y = x; [X, Y]= meshgrid(x,y);
    
    %Fourier domain to define a kernel 
     LF=1/dx; dxf = 1/(dx*N);
     xf= -LF/2:dxf:LF/2-dxf; yf = xf; [XF, YF] = meshgrid(xf,yf);
    
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
     
     first_save_done=0;
     second_save_done=0;
     third_save_done=0;
     
    
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
        
        %kernel construction
        %constants,  % Please check zeta_constants.m
        S0 =0.224959360360906; S2 =0.115241082929453;      
        K=zeros(N,N); eps=0.001;
        cn= (1+eps)*(8*S0*S2); % we choose mobility equal to anisotropy

         for i= 1:N
            for j=1:N
                 thr = atan2(YF(i,j),XF(i,j));
                 gk=(1 + bet*cos(mf*thr))*sqrt((XF(i,j))^2 + (YF(i,j))^2);

                 mu=gk;%mobility

                 alpha= (gk*sqrt(cn) + (sqrt(cn*gk*gk - 8*S0*S2*gk*mu)))*pi/(S2*sqrt(cn));
                 beta = (gk*sqrt(cn) - (sqrt(cn*gk*gk - 8*S0*S2*gk*mu)))*pi/(S2*sqrt(cn));

            xi = sqrt(diff_time)*alpha;
            if  xi <= -2 || 2  <= xi 
                k1 = exp(-xi^2);
            else
                if  -1 <= xi && xi<= 1
                     k1=1;
                else
                     k1= exp(-(-((5*xi^6)/27)+ ((14*xi^4)/9) -((23*xi^2)/9) + 32/27));    
                end
            end

            eta = sqrt(diff_time)*beta;
              %_________________________________________________________
            if  eta <= -2 ||  2 <= eta 
                k2 = exp(-eta^2);
            else
                if  -1 <= eta && eta<= 1
                     k2=1;
                else
                     k2= exp(-(-((5*eta^6)/27)+ ((14*eta^4)/9) -((23*eta^2)/9) + 32/27));
                end
            end

                 K(i,j)=(k1+k2)/2;
            end
         end

        K_area=1;
        U = U_ini; 
        movement=999;
        
        %Numerical solution
        for BMOstep=1:time_steps*cn
 
        U_prev=U;
        
        %Convolution     
        U = ifftshift(ifft2((fftshift(K)).*(fft2(fftshift(U)))));
        U = real(U); % to exclude imaginary part(Due to matlab limitation)
        
        
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
            X_F(ftp,1);

            while X_F(ftp,1) < TT/cn
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
            text(-3.5,-3.5,[' Actual time = ', num2str(BMOstep*diff_time)],'FontSize',20);
            text(-3.5,-2.5,[' Scaled time = ', num2str(BMOstep*diff_time/cn)],'FontSize',20);
            text(-3.5,2.5,[' gpi = ', num2str(gpi)],'FontSize',20);
            text(-3.5,2,[' dti = ', num2str(dti)],'FontSize',20);
            title('Blue: Front tracking, Red: Numerical');
            currFrame = getframe(gcf);  writeVideo(vidObj,currFrame); hold off
         end
        

        %save file for plotting
        if 0.25 >= TT/cn 
            save(sprintf('%s_%d %d %d','EJZFD_U_gpi_dti_TT',gpi,dti,2500),'U');
        end
        
        if 0.75 >= TT/cn 
            save(sprintf('%s_%d %d %d','EJZFD_U_gpi_dti_TT',gpi,dti,7500),'U');
        end
        
        if 1.25 > TT/cn 
            save(sprintf('%s_%d %d %d','EJZFD_U_gpi_dti_TT',gpi,dti,12500),'U');
        end
         

        end%end of BMOsteps
        

    end

    
end


if plotting==1 
  close(vidObj);
end

load train
sound(y,Fs)




