
%Front tracking Sov for S shape


plotting =0;  %Yes=1, No=0

%anisotropy choice
% m fold= 1, crystalline=2, ellipse = 3,non-convex wulff=4, stadium = 5
%isotropic =6
% if anisotropy is m fold then provide more info (uncomment next to next line)
anisotropy = 1;
bet=.05; mf=4;%m-fold anisotropy,

N=100; %Number of points to represent curve (approximate) 
final_time=1.2501;

%terminate the main loop, 
%when the curve is smaller than square whose length is "termination"
termination = 2*10^-1;

%want to fix lower limit for diffusion time?
tau_ll = 1; %yes=1 , No =0;
%if yes, then what is limit?
tau_limit = 10^-4;

syms c t n phi(t) dphi(t)
itaration=0;

%Shape function (phi)
% eps = 0.5; % phi(t) = 1 - eps + eps*sqrt(1 - eps + eps*(t^2));
%phi(t)=(norm(t))^(2/3); %for the length discrepancy
phi(t)= 1; %(norm(t))^(1/3); %for the area discrepancy
dphi(t) = 0; %diff(phi(t));

%Values for relaxation function (Omega)
k1=1; k2=0.01;
lambda=0.01;%keep lambda close to 0 for large diffusion time at beggining 

Total_diffusion_time = 0;

g_s = (4*pi+2)/N; %total lentgh of S shape divide by number of points
z = zeros(2,N+1); %Initial setup (length of z will change later)

fileID = fopen('Front tracking details.txt','w');
fprintf(fileID,'_____________________________________________\n\n');
if plotting==1
 vidObj = VideoWriter('Front tracking.avi');open(vidObj); end

if anisotropy ~= 1
    bet=[]; mf=[]; end


%Initial setup : S shape

%uper big circle
m = 1; ax= 2; ay=2;
for i = 0:(pi/3)*g_s: 3*pi/2
    z(1,m)= ax*cos(i);
    z(2,m)=1+ ay*sin(i);
    m = m + 1;
end

%Lower   small  circle
 ax= 1; ay=1;
for i = pi/2:-(2*pi/3)*g_s: -pi
    z(1,m)= ax*cos(i);
    z(2,m)= - 2+ ay*sin(i);
    m = m + 1;
end

%lower horizonttal line
for i = 0: 2*g_s : 1
    z(1,m)= -1-i;
    z(2,m)= -2;
    m = m + 1;
end

%Lower big  circle
 ax= 2; ay=2;
for i = -pi :(pi/3)*g_s: pi/2
    z(1,m)= ax*cos(i);
    z(2,m)= - 2+ ay*sin(i);
    m = m + 1;
end

%uper small circle
 ax= 1; ay=1;
for i = 3*pi/2:-(2*pi/3)*g_s: 0
    z(1,m)= ax*cos(i);
    z(2,m)=1+ ay*sin(i);
    m = m + 1;
end

%lower horizonttal line
for i = 0: 2*g_s : 1
    z(1,m)= 1+i;
    z(2,m)= 1;
    m = m + 1;
end

z_ini=z; % for plotting
N=m-1; %change size of N according to distribution of the points

fprintf('No. of points = %d \n',N)

X_F=zeros(1,N+2); % will represent front, X coordinate 
Y_F=zeros(1,N+2); % will represent front, Y coordinate 

fprintf(fileID,'%12s \n', 'No. of points');
fprintf(fileID,'%12.0f  \n',N);
fprintf(fileID,'\n \n %12s  %11s %11s\n', 'Iteration', 'diff_time','TOTAL time');


while Total_diffusion_time < final_time
   
%Step 1..
%The discretized arc(p), Length of arc(L) and tangential vector(tv)   
 [r, rd, L, tv] = FA_discretizedArc(z,N); 
 
%tangential angle "nu"
 [nu, nud] = FB_TangentialAngle(tv,N);

% average normal vector(avt) and normal vector field
 nd = FC_avgTVandNormalvf(tv,N); 
  
%Curvature and avarage curvature
  k = zeros(1,N); kd = zeros(1,N);
  k(1) = (nu(2)-nu(N+3))./(2*r(1)); 
  k(2:N) = (nu(3:N+1)-nu(1:N-1))./(2*r(2:N)); 
  kd(1:N-1) = (k(2:N) + k(1:N-1))/2; kd(N) = (k(1) + k(N))/2;  
   
%theta % angle between normal and Y axis 
  th=FD_find_angle(nd); 

%flow equation   %beta= mu*(g"+g)k
[sigma,stiff] = FE_flow_eq(anisotropy,th,bet,mf,N);

%mobility
mu=sigma;
   

% mu=sigma;  
W= mu.*stiff; 
beta = W.*k;
 
%Step 2..finding alpha
avekbeta = dot(k.*beta,r(1:N))/L; %<kbeta>:the average curvature*beta for omega value
omega=k1+ k2.*avekbeta;%relaxation function (Omega)
 
%The tangential component of the velocity vector
[alpha]=FF_findAlpha(r,rd,beta,k,kd,phi,dphi,L,N,omega);

%step 3..new time curve %Semi-implicit numerical scheme
a = alpha./(2*rd); b = W./rd;
aplus(1:N) = b(1:N)/r(2:N+1) + a(1:N);
aminus = b./r(1:N) - a;
a0 = -(aminus + aplus);
abs_alp = abs(alpha); 

%Diffusion time
%Diffusion time
if tau_ll == 1
   tau = tau_limit;
else
   tau = min(r)/((4 + 4*lambda)*(max(W)/min(r) + max(abs_alp)/2));
end

%The matrix for semi-impliit numerical scheme   
X = zeros(N,N); 
  for i = 2:N-1
    X(i,i-1) = -aminus(i)*tau;
    X(i,i) = 1 - a0(i)*tau;
    X(i,i+1) = -aplus(i)*tau;
  end
%Boundary points
    X(1,1) = 1 - a0(1)*tau; X(1,2) = -aplus(1)*tau; X(1,N) = -aminus(1)*tau;
    X(N,1) = -aplus(N)*tau; X(N,N-1) = -aminus(N)*tau; X(N,N) = 1 - a0(N)*tau; 
    
% V=(g+g")k+0, Linear simultaneous equation
 Ztmp = zeros(2,N); 
 for i = 1:2
    Ztmp(i,1:N) = z(i,1:N);% + nd(i,1:N).*tau.*F(1:N); 
 end   

%New curve at next timestep
  z(1,1:N) = X\Ztmp(1,1:N)'; z(2,1:N) = X\Ztmp(2,1:N)';
  z(:,N+1) = z(:,1); %Periodic condition
  
  
Total_diffusion_time = Total_diffusion_time +tau;
itaration=itaration+1;

%plotting

if plotting ==1 
     for i = 1:N
       quiver(z(1,i),z(2,i),nd(1,i),nd(2,i),'linewidth',1); hold on; end

     plot(z_ini(1,:),z_ini(2,:),'b','LineWidth',1); hold on; %inital shape
     plot(z(1,:),z(2,:),'red','LineWidth',2.5); hold on ;%current front
 
     axis equal, xlabel('X axis'),axis([-4 4 -4 4]), grid 
     text(0,3.5,[' Time = ', num2str(Total_diffusion_time)],'FontSize',20);
     text(-3.5,3.5,[' \deltat = ', num2str(tau)],'FontSize',20);
     title('Front tracking')
     currFrame = getframe(gcf);  writeVideo(vidObj,currFrame); hold off
end
    

%save data
Export_data=[itaration tau Total_diffusion_time];
fprintf(fileID,'%12.0f %12.6f %12.8f\n',Export_data);

X_F(itaration,:)=[Total_diffusion_time z(1,:)];
Y_F(itaration,:)=[Total_diffusion_time z(2,:)];

%breaking condition
if abs(min(z(1,:))-max(z(1,:))) < termination  &&...
        abs(min(z(2,:))-max(z(2,:))) < termination
    break
end

end

%save whole front tracking result
save(sprintf('%s_%d ','X_F',N),'X_F');
save(sprintf('%s_%d ','Y_F',N),'Y_F');




%save specific result for plotting purpose
ftp=1;
%saving for time 0.25
while X_F(ftp,1) < 0.25
      ftp=ftp+1;
end
XX_F1=X_F(ftp,:); %copy perticular row 
YY_F1=Y_F(ftp,:); 


%saving for time 0.75
while X_F(ftp,1) < 0.75
      ftp=ftp+1;
end
XX_F2=X_F(ftp,:); %copy perticular row 
YY_F2=Y_F(ftp,:); 



%saving for time 1.25
while X_F(ftp,1) < 1.25
      ftp=ftp+1;
end
XX_F3=X_F(ftp,:); %copy perticular row 
YY_F3=Y_F(ftp,:); 

XX_F=[XX_F1;XX_F2;XX_F3];
YY_F=[YY_F1;YY_F2;YY_F3];

X_F=XX_F;
Y_F=YY_F;

save('X_F3.mat','X_F'); 
save('Y_F3.mat','Y_F'); 



if plotting ==1
 close(vidObj); end
fclose(fileID);







