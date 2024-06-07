%%%%% This function file computes the receiver vector to get seismograms at
%%%%% a receiver location (x0r,y0r)


function f = gerRecvector(HFM,var_para)
% Defining the l_i using a Gaussian function
G =@(x,y,x0,y0,sigma) 1/(2*pi*var_para.sigm^2)*exp(-( (x-x0).^2 + (y-y0).^2)./(2*sigma^2) );

sigma = 10;% width of the Gaussian
var_para.sigx = 120;
var_para.sigy = 50;
var_para.x0r;
var_para.y0r;
fx = G(HFM.P(:,1),HFM.P(:,2),var_para.x0r,var_para.y0r,sigma);%
fy = G(HFM.P(:,1),HFM.P(:,2),var_para.x0r,var_para.y0r,sigma);%
fx = fx./trapz(fx); % Normalizing
fy = fy./trapz(fy); % Normalizing
np=size(HFM.P,1); nn=2*np;
% f at centers of mass
f1=(fx(HFM.elem(:,1))+fx(HFM.elem(:,2))+fx(HFM.elem(:,3)))/3; 
f2=(fy(HFM.elem(:,1))+fy(HFM.elem(:,2))+fy(HFM.elem(:,3)))/3; 
% triangles area
 area=getarea(HFM.P,HFM.elem);
% % assembly
f1=f1.*area/3; f2=f2.*area/3;
f=zeros(nn,1);

f(1:2:nn)=full(sparse(HFM.elem(:,1),1,f1,np,1)+sparse(HFM.elem(:,2),1,f1,np,1)...
          +sparse(HFM.elem(:,3),1,f1,np,1));
f(2:2:nn)=full(sparse(HFM.elem(:,1),1,f2,np,1)+sparse(HFM.elem(:,2),1,f2,np,1)...
          +sparse(HFM.elem(:,3),1,f2,np,1)); 

f = f(HFM.free_dofs,1);

end     

   

function [area,grad1,grad2,grad3]=getarea(p,t)

x21=p(t(:,2),1)-p(t(:,1),1); y21=p(t(:,2),2)-p(t(:,1),2); 
x32=p(t(:,3),1)-p(t(:,2),1); y32=p(t(:,3),2)-p(t(:,2),2);
x31=p(t(:,3),1)-p(t(:,1),1); y31=p(t(:,3),2)-p(t(:,1),2);
% triangles area
area=abs((x21.*y31-y21.*x31)/2);  
% gradients of basis functions
grad1=.5*[-y32./area x32./area]; grad2=.5*[y31./area -x31./area]; grad3=.5*[-y21./area x21./area]; 
end