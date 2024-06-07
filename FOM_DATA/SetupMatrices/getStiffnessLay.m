% Copyright (c) 2016, Jonas Koko
% All rights reserved.
% For details: Koko J., Fast MATLAB assembly of FEM matrices in 2D and 3D using cell-array approach, International Journal of Modeling, Simulation, and Scientific Computing Vol. 7, No. 3 (2016)

function K = getStiffnessLay(HFM,var_para)
    p = HFM.P;

    K_lam = cell(5,1);
    K_mu = cell(5,1);
    Ki = cell(5,1);
for layer = 1:5
    id = find(HFM.elem(:,4)==layer);
    elem = HFM.elem(id,1:3);
    
    np=size(p,1); nn=2*np; nt=size(elem,1);

    % Gradients of basis functions
    [area,grad1,grad2,grad3] = getGradients(p,elem);

    % Adding alternating horizontal and vertical components
    B=cell(3,6); [B{:,:}]=deal(sparse(nt,1)); 

    B{1,1}=grad1(:,1); B{1,3}=grad2(:,1); B{1,5}=grad3(:,1);
    B{2,2}=grad1(:,2); B{2,4}=grad2(:,2); B{2,6}=grad3(:,2);
    B{3,1}=grad1(:,2); B{3,2}=grad1(:,1); B{3,3}=grad2(:,2); 
    B{3,4}=grad2(:,1); B{3,5}=grad3(:,2); B{3,6}=grad3(:,1);

    % Plane strain condition Stifness matrix C
    C1 =[1 1       0; 
       1      1  0; 
       0        0         0];
    K1 = getstfmat(B,C1,p,elem,area);

    C2 =[2     0       0; 
         0      2      0; 
         0        0      1];
    K2 = getstfmat(B,C2,p,elem,area);

    K_lam{layer} = var_para.lami(layer)*K1;
    K_mu{layer} = var_para.mui(layer)*K2;   
    Ki{layer} = K_lam{layer} + K_mu{layer};

   
end
K = Ki{1} + Ki{2} + Ki{3} + Ki{4} + Ki{5};
end


function K = getstfmat(B,C,p,elem,area)
np=size(p,1); nn=2*np; nt=size(elem,1);


E=cell(3,6);                           
for i=1:3
    for j=1:6                    
         E{i,j}=C(i,1)*B{1,j}+C(i,2)*B{2,j}+C(i,3)*B{3,j};
    end
end

it1 = [1 1 2 2 3 3]; it2 = [1 2 1 2 1 2];
K = sparse(nn,nn);                      

% Entries below the diagonal
for i=1:6
    ik=2*(it1(i)-1)+it2(i); it=2*(elem(:,it1(i))-1)+it2(i);
    for j=1:i-1
        jl = 2*(it1(j)-1)+it2(j); jt=2*(elem(:,it1(j))-1)+it2(j);
        Rij = B{1,ik}.*E{1,jl}+B{2,ik}.*E{2,jl}+B{3,ik}.*E{3,jl};
        K = K+sparse(it,jt,area.*Rij,nn,nn);
    end
end

% Adding upper triangular part (Exploiting symmetry)
K = K + K.';      

% Adding the diagonal part
for i=1:6                              
    ik=2*(it1(i)-1)+it2(i); it=2*(elem(:,it1(i))-1)+it2(i);
    Rij=B{1,ik}.*E{1,ik}+B{2,ik}.*E{2,ik}+B{3,ik}.*E{3,ik};
    K=K+sparse(it,it,area.*Rij,nn,nn);
end

end

function [area,grad1,grad2,grad3]=getGradients(p,elem)

x21=p(elem(:,2),1)-p(elem(:,1),1); y21=p(elem(:,2),2)-p(elem(:,1),2); 
x32=p(elem(:,3),1)-p(elem(:,2),1); y32=p(elem(:,3),2)-p(elem(:,2),2);
x31=p(elem(:,3),1)-p(elem(:,1),1); y31=p(elem(:,3),2)-p(elem(:,1),2);
% Triangulation area
area=abs(x21.*y31-y21.*x31)/2;  
% Gradients of the linear (P1) basis 
grad1=.5*[-y32./area x32./area]; 
grad2=.5*[y31./area -x31./area]; 
grad3=.5*[-y21./area x21./area]; 
end