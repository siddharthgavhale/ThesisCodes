
%kernel construction 

function K = F_kernel(anisotropy,kernel_choice,beta,phi,m,...
                             diff_time,N,X,Y,XF,YF,eps)
                          
 if anisotropy == 1 % m=fold
     if kernel_choice == 1 %gaussian
        K = (1/(4*pi*(diff_time))).*exp(-(X.^2+Y.^2)./(4*(diff_time)));
     end
     
     if kernel_choice == 2 %EJZ PD
        %Kernel EJZ  PD

        m0 = 0.34029423827512594869;
        m2 = 0.37372385180378853153;
        K=zeros(N,N); 
        for i= 1:N
            for j=1:N
                th = atan2(Y(i,j),X(i,j)) - pi/2; % Polar form: angle 
                rt = sqrt((Y(i,j)^2 + X(i,j)^2)/diff_time); % Polar form

                 ani= 1+ beta*cos(m*(th+phi));% anisotropy
                 stiff= 1 + (beta*cos(m*(th+phi))) - (beta*m^2*cos(m*(th+phi)));% stiff=g+g"
                         
                         
                 mu=ani;%mobility = anisotropy

                Beta = sqrt(m2/(2*m0*mu*stiff));
                epw = rt*Beta;

                if 0 < epw && epw < 2
                    eta = exp(-1/((epw^2)*(epw-2)^2));
                    alpha = Beta/(4*m0*mu);
                    K(i,j)=alpha*eta;
                end

                    
            end
        end
        K=(K./diff_time);
     end 
     
     
     if kernel_choice == 3 %EE kernel
         K=zeros(N,N); 
        for i=1:N
            for j=1:N  
        fun=@(t) (1+ beta*cos(m*(t+phi-(pi/2))) - beta*m^2*cos(m*(t+phi-(pi/2))))...
            .*(exp(-((X(i,j).*cos(t) + Y(i,j).*sin(t)).^2)/(4*diff_time)))...
            .*(exp(-((X(i,j).*sin(t) - Y(i,j).*cos(t)).^2)/(4*diff_time*eps.^2)));
               K(i,j)=integral(fun,0,2*pi); 
            end
        end
         K=(1/(16*eps*sqrt(pi))).*K;
         K=(K./diff_time);
     end
     
     if kernel_choice == 4 %Bonnetier kernel
        K=zeros(N,N);
        for i=1:N
            for j=1:N  
                  thr = atan2(YF(i,j),XF(i,j));
                  g = (1 + beta*cos(m*(thr+phi)))*sqrt((XF(i,j))^2 + (YF(i,j))^2);
                  K(i,j)=exp(-4*pi*pi*diff_time*g^2);   
            end
        end   
     end
     
     if kernel_choice == 5 %EJZ FD
        %constants,  % Please check zeta_constants.m
        S0 =0.224959360360906; S2 =0.115241082929453;      
        K=zeros(N,N); eps=0.01;
        cn= (1+eps)*(8*S0*S2); % we choose mobility equal to anisotropy

         for i= 1:N
            for j=1:N
                 thr = atan2(YF(i,j),XF(i,j));
                 gk=(1 + beta*cos(m*(thr+phi)))*sqrt((XF(i,j))^2 + (YF(i,j))^2);
    %             gk=(1 + beta*cos(m*(thr+phi)));
                 mu=gk;%mobility
        
                 alpha= (gk*sqrt(cn) + (sqrt(cn*gk*gk - 8*S0*S2*gk*mu)))*pi/(S2*sqrt(cn));
                 Beta = (gk*sqrt(cn) - (sqrt(cn*gk*gk - 8*S0*S2*gk*mu)))*pi/(S2*sqrt(cn));

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

            eta = sqrt(diff_time)*Beta;
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
         
     end
         
     
 end

    
end