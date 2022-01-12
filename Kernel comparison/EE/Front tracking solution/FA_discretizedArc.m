
%length of arc and points to represent arc
function [r, rd, L, tv]=FA_discretizedArc(z,N)

%Distance between each points
 r = zeros(1,N+1);
 r(2:N) = ((z(1,2:N)-z(1,1:N-1)).^2 + (z(2,2:N)-z(2,1:N-1)).^2).^(1/2);
 r(1) = norm(z(:,1)-z(:,N)); r(N+1) = r(1);  
 
%the average arc-lgth, which will be used for adjusted tangential velocity
 rd = zeros(1,N);
 rd(1:N-1) = (r(1:N-1) + r(2:N))/2; rd(N) = (r(N) + r(1))/2;
 
%total length of the curve
 L = sum(r(1:N)); 

% tangent vector
 tv = zeros(2,N+1); 
 for i = 1:2
    tv(i,2:N) = (z(i,2:N)-z(i,1:N-1))./r(2:N); 
 end
    tv(:,1) = (z(:,1)-z(:,N))./r(1);
    tv(:,N+1) = tv(:,1);
end
