
%EE kernel
%particle on substrate code
%Characteristic type: without modification
clear;

disp('EE')

%File to export data
fileID = fopen('Z4 EE.txt','w');
fprintf(fileID,'\n \n%5s  %8s %9s %9s %11s %15s %16s\n',...
   'final time', 'Ini dt','End dt','dx','Error','Cont. A(N)','Cont. A(A)');
fprintf(fileID, '\n');

plotting=0;%want plotting yes=1, No=0

gpiSV=9; gpiEV=12; %range of grid points %gpiSV=9; gpiEV=12;
dtiSV=0; dtiEV=7; %range of diffusion time %dtiSV=0; dtiEV=7;

gammaSP=1.5; gammaSV=1; % surface tensions
L=10; %L*L domain
final_time = 50;

%anisotropy choice
% m fold= 1, 
%coming soon......crystalline=2, ellipse = 3, ellipse without smooth end =4,
% if anisotropy is m fold then provide more info (uncomment next to next line)
anisotropy = 1;
beta=0.05; phi=.0; m=4;%m-fold anisotropy,

%kernel: Gaussian=1, EJZ PD=2, EE=3, Bonnetier=4, EJZ FD=5
kernel_choice=3;
epsilon = 0.01; %for EE kernel

%for which gpi value need contact angle plot
ca_plt=10;

%initial shape of particle
%semi-circle=1, square=2, rectangel=3, half asteroid=4, half ellipse=5,
%WUlff shape=6
initial_shape=2;

%Contact angle parameters (CAP)
%How many points (representing particle bouondary) should be exclude
%(CAP)from calculation (near substrate)
exclude_pts=3;

%(CAP)Consider all the points which have verticle distance less that 'THIS' (from substrate)
dist_above = 0.15; % this parameter must be change according to final shape and anisotropy

movement_ini=9999; %to track how much particle change within one time step
syms th %symbolic theta
Simulation_details=zeros((gpiEV-gpiSV+1)*(dtiEV-dtiSV+1),9); %to save all the detalis for error analysis
temp1=1;% to count total computation pairs gpi and dti
TT=0; %total time
iteration=0; %iteration count
DXX=zeros(1,gpiEV-gpiSV+1);
DTT=zeros(1,dtiEV-dtiSV+1);

%ays1 & ays2 are solution of anisotropic young equaiton (contact angle)
%and parametric form of Wulff shape (xw,yw)
% Area of stable shape (Winterbottom form) Wint_area
if anisotropy ~= 1
    beta=[]; phi=[]; m=[];
end
[ays1,ays2,xw,yw,Wint_area] = B_solution_young_eq (anisotropy,beta,phi,m,th,gammaSP,gammaSV);

L_1=zeros(gpiEV-gpiSV+1,dtiEV-dtiSV+1);%To export for plotting
if plotting==1
   vidObj = VideoWriter('Evolution of particle.avi'); open(vidObj); end

for gpi=gpiSV:gpiEV %Grid Point ind
    
   %Initial setup of  domain for Numerical result  
   N=2^gpi;%grid points in one direction      
   dx = L/N; x = -L/2:dx:L/2-dx; 
   y = x; [X, Y]= meshgrid(x,y);
   DXX(gpi-gpiSV+1)=dx;
   
   %to save contact angle for plottng
   if gpi==ca_plt
       CA_EE=zeros(1,dtiEV-dtiSV+1);
   end
   
   if kernel_choice == 4 || kernel_choice == 5 %for bonnetier and EJZ FD
    %Fourier domain to define a kernel 
     LF=1/dx; dxf = 1/(dx*N);
     xf= -LF/2:dxf:LF/2-dxf; yf = xf; [XF, YF] = meshgrid(xf,yf); 
   else
       XF=[];YF=[];% assign empty set to avoid error whhile calling kernel 
   end
          
  %Initialization: particle(P), substrate(S) & vapor(V) & calculate
  %initial area of particle (analytical,in grid scale)
  [P,S,V,area] = C_initialization(initial_shape,N,dx,X,Y);
   P_initial=P; %Will  be used for  different dti & To plot Initial shape 
   V_initial=V; %Will  be used for  different dti 
         
   %Analytical solution for given gpi: 
   % xw & yw parametric form (with inflation of area).
   % Inflated WUlff shape in 0-1 matrix form 'Wm'  
   Ratio_area= area/(Wint_area/dx^2); %To inflate the analytical solution
   xw_n=xw.*sqrt(Ratio_area); yw_n=yw.*sqrt(Ratio_area);%analytical solution for given set up

    Wm=zeros(N,N);%This matrix represent analyttical solution in 0-1 form
    for i=1:N
        for j=1:N
             [in,on] = inpolygon(X(i,j),Y(i,j),xw_n,yw_n);
                   if on == 1 || in == 1 
                       Wm(i,j)=1; end
        end
    end
    
    %Until which row, we have substrate? Since half domain fixed for S
    sub_level = N/2 + 1;
    
    for dti=dtiSV:dtiEV %Diffusion Time ind

        fprintf('space = %d\t time = %d\n', gpi, dti);     

        %Initial setup of time related para. for Numerical result  
        diff_time=1/(4*2^dti); %diffusion time
        Ini_dt = diff_time; % For exporting data (in modified case, final diffusion time will be different)
        time_steps=final_time/diff_time; %maximum time steps.
        movement=movement_ini;
        DTT(dti-dtiSV+1)=diff_time;
        
    %Initilize P and V
    P = P_initial;
    V = V_initial;
     
    %Since only one kernel option given therefore not writtting separate
    %kernel
    Gau=(1/(4*pi*(diff_time))).*exp(-(X.^2+Y.^2)./(4*(diff_time))); 
    K=F_kernel(anisotropy,kernel_choice,beta,phi,m,diff_time,N,X,Y,XF,YF,epsilon); 
      
    %first convolution with substrate interface 
    ph1 = (dx^2).*ifftshift(ifft2(fft2(fftshift(Gau))...
              .*(fft2(fftshift((gammaSP-gammaSV).*S)))))./sqrt(diff_time);
              

    %Main loop of thresholding
       while movement > 0 && final_time > TT
          
           %convolution either with Kernel in Fourier form of regular form
           if kernel_choice == 4 || kernel_choice == 5
                ph2 = ifftshift(ifft2((fftshift(K))...
                     .*(fft2(fftshift(V-P)))))./sqrt(diff_time);
                 ph2 = real(ph2);   
           else
                ph2 = (dx^2).*ifftshift(ifft2(fft2(fftshift(K))...
                     .*(fft2(fftshift(V-P)))))./sqrt(diff_time);
           end
    
          M=ph1+ph2;

           %exclude all the values (from calculation of thresholding height)
           %which intersecting with substrate
           M(1:sub_level,:) = 99999999999; 
           sortlist=sort(M(:)); %sorting
           delta=sortlist(area); %thresholding height, it preserves area.
           
           %iteratin  details and total time
           iteration=iteration+1; TT= TT+diff_time;
    
           %%%plotting during evolution
           if plotting==1 
                contour(X,Y,P_initial,[0.5,0.5],'g','LineWidth',2);  hold on             
                contour(X,Y,ph1+ph2,[delta,delta],'b','LineWidth',2); hold on
                plot(xw_n,yw_n,'Color','black','LineWidth',2); 
                line([-3.5,3.5],[0,0],'Color','red','LineWidth',4)
                     text(1,2,[' Time = ',num2str(TT)],'FontSize',20);
                     text(-2.7,2.2,[' dti = ',num2str(dti)],'FontSize',20);
                     text(-2.75,2.5,[' gpi = ',num2str(gpi)],'FontSize',20);
                     text(1,2.4,[' \deltat = ',num2str(diff_time)],'FontSize',20);
                axis equal, grid on, xlabel('X'), 
                axis([-3.5 3.5 0 3.5]) 
                currFrame = getframe(gcf); writeVideo(vidObj,currFrame); hold off;
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
  
       
    %Finding contact angle (CA) left and right, CA_L & CA_R
    [CA_L, CA_R] = D_contact_angle...
                    (dist_above,exclude_pts,M,delta,dx,X,Y,ays1,ays2,sub_level);
    
     %save left contact angle for plotting: along with other kernels
    if gpi==ca_plt
        CA_EE(dti-dtiSV+1)=CA_L;
    end
                
                
    %Error analysis    
         % Creat matrix for Error
         Error=(P-Wm); 
         
         %exclude value at substrate and one more row above it
         Error([sub_level sub_level+1],:)=zeros(2,N);
         Error=sum(abs(Error),'all').*dx^2;
         Error=Error/(area*dx^2);
         
    %Export data to text file
        Export_data=[TT Ini_dt diff_time dx Error [CA_R,CA_L] (ays1) (ays2)];
        fprintf(fileID,'%7.3f %11.5f %10.5f %10.5f %10.4f %9.2f & %5.2f %8.2f & %5.2f\n',Export_data);
         
    %To export for plotting
        L_1(gpi-gpiSV+1,dti-dtiSV+1) = Error;
        Simulation_details(temp1,:) =  Export_data;
        temp1=temp1+1;
        
       %reset to zero
       iteration=0; TT= 0;
    
    end

    fprintf(fileID, '\n');
    
    %save left contact angle for plotting: along with other kernels
    if gpi==ca_plt
      save('CA_EE.mat','CA_EE');
    end

end

if plotting==1
close(vidObj);
end

fclose(fileID);

figure()
linS = {'--',':','-.','-','--*',':o','-.v','--^'};
%saved data contain all the infromrmation, extracting last 4 values of
%error by starting specific kk
kk=1;
for i=kk:length(DXX)
    dx=DXX(i);
 
name=['dx = ' num2str(round(dx,5))];
%line type is different and need different format in matlab hence if-else
%condition
    if i <= 1
loglog(DTT,L_1(i,:),'DisplayName', name,'linestyle',linS{i},'Linewidth',7.5),
    else
loglog(DTT,L_1(i,:),linS{i},'DisplayName', name,'Linewidth',7.5),
    end
    
hold on
set(gca,'FontSize',60)
xlabel('\deltat','FontSize',55),ylabel('error','FontSize',55),
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
lgd = legend('Location','southeast');
lgd.FontSize = 60; 

end
legend('AutoUpdate','off')

title('\fontsize{40}  Elsey-Esedoglu (\epsilon = 0.01)');
xlim([10^-3 10^0])

%maximise the plot and save
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
saveas(gcf,'Fig6_4.png')

hold off
 







