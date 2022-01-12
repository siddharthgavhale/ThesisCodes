
%flow equation   %beta= mu.*(g"+g)k
 
 function [sigma,stiff] = FE_flow_eq(anisotropy,th,bet,mf,N)
     
 if anisotropy == 1  %m fold 
    sigma(1:N) =  bet.*cos(mf.*th(1:N)) + 1;
    stiff(1:N) =  bet.*cos(mf.*th(1:N)) + 1 - bet.*mf.^2..*cos(mf.*th(1:N));% stiff=g+g" 
 end

 if anisotropy == 2 %crystalline
     epsi=0.01;
    sigma(1:N) =  (epsi.^2 + cos(th(1:N)).^2).^(1/2) + (epsi.^2 + sin(th(1:N)).^2).^(1/2);
    stiff(1:N) = cos(th(1:N)).^2/(epsi.^2 + sin(th(1:N)).^2).^(1/2)...
                   - cos(th(1:N)).^2/(epsi.^2 + cos(th(1:N)).^2).^(1/2)...
                   + sin(th(1:N)).^2/(epsi.^2 + cos(th(1:N)).^2).^(1/2)...
                   - sin(th(1:N)).^2/(epsi.^2 + sin(th(1:N)).^2).^(1/2)...
                   + (epsi.^2 + cos(th(1:N)).^2).^(1/2) + ...
                   (epsi.^2 + sin(th(1:N)).^2).^(1/2)...
                    - (cos(th(1:N)).^2.*sin(th(1:N)).^2)/(epsi.^2 +...
                    cos(th(1:N)).^2).^(3/2)...
                    - (cos(th(1:N)).^2.*sin(th(1:N)).^2)/...
                    (epsi.^2 + sin(th(1:N)).^2).^(3/2);% stiff=g+g" 
 end
 
 if anisotropy == 3 %ellipse
     sigma(1:N) = (4.*cos(th(1:N)).^2 + sin(th(1:N)).^2).^(1./2);
     stiff(1:N) = (3.*sin(th(1:N)).^2)./(4.*cos(th(1:N)).^2 +...
                    sin(th(1:N)).^2).^(1./2) - (3.*cos(th(1:N)).^2)...
                    ./(4.*cos(th(1:N)).^2 + sin(th(1:N)).^2).^(1./2) +...
                    (4.*cos(th(1:N)).^2 + sin(th(1:N)).^2).^(1./2) - ...
                    (9.*cos(th(1:N)).^2.*sin(th(1:N)).^2)./(4.*cos(th(1:N)).^2 ...
                    + sin(th(1:N)).^2).^(3./2);
 end
 
 if anisotropy == 4 % non convex wulff 
    sigma(1:N) = sin(2.*th(1:N)) + 1./2;
    stiff(1:N) = 1./2 - 3.*sin(2.*th(1:N));
 end
 
 if anisotropy == 5 % strip shape
    sigma(1:N) = 1./4 + (5./4).*(abs(sin(th(1:N)+pi./4)));
    stiff(1:N) = (5.*abs(sin(th(1:N) + pi./4)))./4 + ...
                  (5.*dirac(sin(th(1:N) + pi./4)).*cos(th(1:N) + pi./4).^2)./2 ...
                  - (5.*sign(sin(th(1:N) + pi./4)).*sin(th(1:N) + pi./4))./4 + 1./4;
 end
 
  if anisotropy == 6  %Isotropic 
      bet=0; mf=1;
    sigma(1:N) =  bet.*cos(mf.*th(1:N)) + 1;
    stiff(1:N) = bet.*cos(mf.*th(1:N)) + 1 - bet.*mf.^2..*cos(mf.*th(1:N));% stiff=g+g" 
 end
 
     
      
 end