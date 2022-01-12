%average tangential velocity and normal vector field

function nd=FC_avgTVandNormalvf(tv,N)

%averaged tangential vector 
 avt = zeros(2,N); 
 for i = 1:2
   avt(i,1:N-1) = (tv(i,1:N-1) + tv(i,2:N))./2;
 end
   avt(:,N) = (tv(:,1) + tv(:,N))./2;
  
%normal vector   
    nd = zeros(2,N); 
    nd(:,1:N) = [0 -1;1 0]*avt(:,1:N); %rotated by the apropriate matrix
    %nd(:,1) = [0 -1;1 0]*avt(:,N);
end