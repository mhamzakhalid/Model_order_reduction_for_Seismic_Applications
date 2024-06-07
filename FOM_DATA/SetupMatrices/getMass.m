% Copyright (c) 2016, Jonas Koko
% All rights reserved.
% For details: Koko J., Fast MATLAB assembly of FEM matrices in 2D and 3D using cell-array approach, International Journal of Modeling, Simulation, and Scientific Computing Vol. 7, No. 3 (2016)


function M = getMass(HFM,rho)
p = HFM.P;
np = size(HFM.P,1);
M = sparse(2*np,2*np); 
for layer = 1:5
id = find(HFM.elem(:,4)==layer);
elem = HFM.elem(id,1:3);

nt=size(elem,1); 
area = getarea(p,elem); 

globalid = 2*elem(:,[1 1 2 2 3 3])-kron(ones(nt,1),[1,0,1,0,1,0]); % indexing for hor. and vert. comp
Yid_3D = reshape(repmat(globalid,1,6)',6,6,nt); % Reshaping for global sparse matrix
Xid_3D = permute(Yid_3D,[2 1 3]); % Reshaping for transpose part
Mloc =rho(layer)*[2 0 1 0 1 0;%
           0 2 0 1 0 1;
           1 0 2 0 1 0;
           0 1 0 2 0 1;
           1 0 1 0 2 0;
           0 1 0 1 0 2]/12; % Local mass matrix
Mloc3D = repmat(Mloc,[1 1 nt]);  % Replication into an array     6 X 6 X No. elem
[nx,ny,nz] = size(Mloc3D); % Computing the size
area3D = reshape(area',1,1,nz); % Reshaping in array  1 X 1 X No. elem
area3D = area3D(ones(nx,1),ones(ny,1),:); % Reshaping for each elem 6 X 6 X No. elem
Mglob = Mloc3D .* area3D; % Muliplication with the area

M = M + sparse(Xid_3D(:),Yid_3D(:),Mglob(:),2*np,2*np); % Global sparse mass matrix

end

end


function [area]=getarea(p,elem)

x21=p(elem(:,2),1)-p(elem(:,1),1); y21=p(elem(:,2),2)-p(elem(:,1),2); 
x31=p(elem(:,3),1)-p(elem(:,1),1); y31=p(elem(:,3),2)-p(elem(:,1),2);
% Area of the triangle
area=abs((x21.*y31-y21.*x31)/2);  

end

