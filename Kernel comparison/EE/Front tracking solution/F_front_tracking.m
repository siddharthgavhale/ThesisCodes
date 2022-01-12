
%Front tracking Sov with repeat option and analytical solution for elliptic
%anisotropy

  function A = F_front_tracking(ax,ay,final_time,front_pts,gpi,eps,XX,YY)

N=front_pts;
N1=2^gpi;


plotting=0;

if plotting==1
    vidObj = VideoWriter('Front tracking.avi');open(vidObj); end

%Initial setup for analytical solution
z = zeros(2,N+1); 

%Initial Curve
m = 1;
for i = -pi/2 + 2*pi/N:2*pi/N:3*pi/2 + 2*pi/N
    z(1,m)=ax*cos(i);
    z(2,m)=ay*sin(i);
    m = m + 1;
end

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

while Total_diffusion_time <= final_time

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
 sigma(1:N)=sqrt((ay.*cos(th(1:N))).^2 + (ax.*sin(th(1:N))).^2);
 stiff(1:N)= (ax^2*ay^2)./(sigma(1:N)).^3;% stiff=g+g"

% EE_kernel mobility
 mu = zeros(1,N); 
 mu_p=zeros(1,N);
  for i=1:N 
  fun = @(t) ((ax.^2.*ay.^2)./...
              (sqrt((ay.*cos(t-pi./2)).^2 + (ax.*sin(t-pi./2)).^2)).^3)...
             ./(((1-eps.^2) .* (cos(t-th(i))).^2) + eps^2).^(0.5);
   mu_p(i)= integral(fun,0,2*pi);
   mu(i)=1/(0.25*mu_p(i));
  end
   

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
tau = min(r)/((4 + 4*lambda)*(max(W)/min(r) + max(abs_alp)/2));

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
  
  
Total_diffusion_time=Total_diffusion_time +tau;


itaration=itaration+1;

      % %plotting
      if plotting==1
     plot(z(1,:),z(2,:),'red','LineWidth',1.5);
     axis equal, xlabel('X axis'),axis([-4.5 4.5 -4.5 4.5])
      text(1.5,2,[' Time = ', num2str(Total_diffusion_time)],'FontSize',20);
     text(-3.5,-3,[' Diff time = ', num2str(tau)],'FontSize',20);
     title('EJZ FD. Red:Analytical, Blue:Numerical')
     currFrame = getframe(gcf);  writeVideo(vidObj,currFrame); hold off
      end
     

end

F_X=z(1,:);
F_Y=z(2,:);


%sign distance function of closed curve
sdf=zeros(N,N);

%first way
for i=1:N
    for j=1:N
         dis=999;
         %from each segment of front tracking solution find minimum distance
            for k=1:length(F_X)-1 
                di1 = sqrt((XX(i,j)-F_X(k))^2+(YY(i,j)-F_Y(k))^2);
                di2 = sqrt((XX(i,j)-F_X(k+1))^2+(YY(i,j)-F_Y(k+1))^2);

                t = -((F_X(k+1)-F_X(k))*(F_X(k)-XX(i,j))+(F_Y(k+1)-F_Y(k))*(F_Y(k)-YY(i,j)))/((F_X(k+1)-F_X(k))^2+(F_Y(k+1)-F_Y(k))^2);
                if ( t>0 && t<1 )
                      di3 = abs((F_Y(k+1)-F_Y(k))*XX(i,j)-(F_X(k+1)-F_X(k))*YY(i,j)+F_X(k+1)*F_Y(k)-F_Y(k+1)*F_X(k))/sqrt((F_X(k)-F_X(k+1))^2+(F_Y(k)-F_Y(k+1))^2);
                else
                       di3 = 1000;
                end

                    d = min([di1,di2,di3]);             
                    
                  if dis > d 
                     dis=d;
                  end
            end
    % check in side or outside
    [in,on] = inpolygon(XX(i,j),YY(i,j),F_X,F_Y);
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
A=sdf;

%subgrid improvement
AC = A; % a copy of A is needed
for i=2:N1-1
    for j=2:N1-1
        if sign(A(i,j)) == sign(A(i,j+1)) == sign(A(i+1,j)) == sign(A(i,j-1)) ... 
                == sign(A(i-1,j)) == sign(A(i-1,j-1)) == sign(A(i-1,j+1)) ...
                == sign(A(i+1,j-1)) == sign(A(i+1,j+1)) 
            A(i,j)= 0.5*sign(A(i,j));
        else % compAte the area ratio of positive phase on 4 grid cells centered at (i,j)
            
    b=10^-7;   
w1=AC(i-1,j-1); w2=AC(i-1,j); w3=AC(i,j-1);w4=AC(i,j);
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

w1=AC(i-1,j); w2=AC(i-1,j+1); w3=AC(i,j);w4=AC(i,j+1);
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

w1=AC(i,j-1); w2=AC(i,j); w3=AC(i+1,j-1);w4=AC(i+1,j);
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


w1=AC(i,j); w2=AC(i,j+1); w3=AC(i+1,j);w4=AC(i+1,j+1);
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
            A(i,j) = Areacell/4.0 -0.5;
        end
    end
end
A=A+0.5;

%AdjAsting boAndary
A(1,:) = 0; A(N1,:) = 0; 
A(:,1) = 0; A(:,N1) = 0;


if plotting==1
 close(vidObj); end

 end





