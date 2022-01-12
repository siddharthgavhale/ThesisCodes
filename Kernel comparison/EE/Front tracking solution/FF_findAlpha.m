
function  [alpha]=FF_findAlpha(p,rd,beta,k,kd,phi,dphi,L,N,omega)

pN = zeros(1,N); pN = p(1:N);

%second derivative of beta
ddbeta(2:N-1)=((beta(3:N)-beta(2:N-1))./rd(2:N-1)-(beta(2:N-1)-beta(1:N-2))./rd(1:N-2))./p(2:N-1);
ddbeta(1)=((beta(2)-beta(1))./rd(1)-(beta(1)-beta(N))./rd(N))./p(1);
ddbeta(N)=((beta(1)-beta(N))./rd(N)-(beta(N)-beta(N-1))./rd(N-1))./p(N);

f = zeros(1,N);
f(1:N) = ((ddbeta(1:N) + (k(1:N).^2).*beta(1:N)).*dphi(k(1:N))) - k(1:N).*beta(1:N).*phi(k(1:N));

avef = sum(f*pN')/L;

phik = zeros(1,N); phik(1:N) = phi(k(1:N));
avephi = sum(phik*pN')/L;

psi = zeros(1,N);
psi(1:N) = pN(1:N).*phi(k(1:N))*avef/avephi - f(1:N).*pN(1:N) + (L*avephi/N - phi(k(1:N)).*pN(1:N))*omega; 
psi(1)=0;

Psic = cumsum(psi); %partial sum of  psi

     %Adjusted tangential velocity
     alpha = zeros(1,N); 
     alpha(1) = - (Psic*rd')/(L * phi(kd(1)));
     temp = phi(kd(1))*alpha(1);
     alpha(2:N) = (temp + Psic(2:N))./phi(kd(2:N));
        
end