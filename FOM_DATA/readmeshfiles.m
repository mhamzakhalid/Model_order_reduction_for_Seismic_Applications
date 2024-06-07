function HFM = readmeshfiles()
%oldfolder = cd;
%cd ../

%% File loading
filename_table = 'mesh_levelPaper23_table';
filename_nodes= 'mesh_levelPaper23_nodes';
filename_triangles = 'mesh_levelPaper23_triangles';
filename_left = 'mesh_levelPaper23_left';
filename_right = 'mesh_levelPaper23_right';
filename_top = 'mesh_levelPaper23_top';
filename_bottom = 'mesh_levelPaper23_bottom';
filename_interface0 = 'mesh_levelPaper23_interface0';
filename_interface1 = 'mesh_levelPaper23_interface1';
filename_interface2 = 'mesh_levelPaper23_interface2';
filename_interface3 = 'mesh_levelPaper23_interface3';

%% Reading material properties
% format (5 X 7)
% - id layer
% - rho/density (kg/m3)
% - vs (km/s)
% - vp (km/s)
% - mu/G/shear modulus
% - kappa/bulk modulus
% - lambda(lame parameter)


fid = fopen(filename_table,'rt');
materProp = textscan(fid, '%f%f%f%f%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);

HFM.materProp = materProp;
% format (nt X 4)
% - First three: triangulation tags (indentifier)
% - Last: layer tag (indentifier)

fid = fopen(filename_triangles,'rt');
triang = textscan(fid, '%f%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);
triang = cell2mat(triang);
HFM.nt = triang(1,1);
HFM.elem = triang(2:end,:);

% format (np X 2)
% - First: triangulation tags (indentifier)
% - Last: nodes

fid = fopen(filename_nodes,'rt');
nodes = textscan(fid, '%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);
nodes = cell2mat(nodes);
HFM.nd = nodes(1,1);
HFM.P = nodes(2:end,:);

% format (bdl X 3) [Tri_id, edge1, edge2]

fid = fopen(filename_left,'rt');
bdleft = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);
bdleft = cell2mat(bdleft);
    %bd_ids_left = unique(HFM.elem(bdleft(:,1),1));
bd_ids_left = getboundaryVertIDs(HFM.elem,bdleft);

% format (bdr X 3)

fid = fopen(filename_right,'rt');
bdright = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);
bdright = cell2mat(bdright);
        %bd_ids_right = unique(HFM.elem(bdright(:,1)));
bd_ids_right = getboundaryVertIDs(HFM.elem,bdright);

% format (bdt X 3)

fid = fopen(filename_top,'rt');
bdtop = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);
bdtop = cell2mat(bdtop);
       %bd_ids_top = unique(HFM.elem(bdtop(:,1)));
bd_ids_top = getboundaryVertIDs(HFM.elem,bdtop);

% format (bdb X 3)

fid = fopen(filename_bottom,'rt');
bdbot = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);
bdbot = cell2mat(bdbot);
      %bd_ids_bot = unique(HFM.elem(bdbot(:,1)));
bd_ids_bot = getboundaryVertIDs(HFM.elem,bdbot);



fid = fopen(filename_interface0,'rt');
interface0 = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);

HFM.interface0 = cell2mat(interface0);

fid = fopen(filename_interface1,'rt');
interface1 = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);

HFM.interface1 = cell2mat(interface1);

fid = fopen(filename_interface2,'rt');
interface2 = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);

HFM.interface2 = cell2mat(interface2);

fid = fopen(filename_interface3,'rt');
interface3 = textscan(fid, '%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ');
fclose(fid);

HFM.interface3 = cell2mat(interface3);


HFM.bd = cell(4,1);
HFM.bd{1} = bdleft;
HFM.bd{2} = bdright;
HFM.bd{3} = bdtop;
HFM.bd{4} = bdbot;

HFM.bdbot = bd_ids_bot;
HFM.bdtop = bd_ids_top;
HFM.bdleft = bd_ids_left;
HFM.bdright = bd_ids_right;

rec_x1 = 3296;
rec_y1 = 10000;
[rec_dist1, rec_id1] = min((HFM.P(:,1) - rec_x1).^2 + (HFM.P(:,2) - rec_y1).^2);
rec_x2 = 6615;
rec_y2 = 10000;
[rec_dist2, rec_id2] = min((HFM.P(:,1) - rec_x2).^2 + (HFM.P(:,2) - rec_y2).^2);
rec_x3 = 8911;
rec_y3 = 10000;
[rec_dist3, rec_id3] = min((HFM.P(:,1) - rec_x3).^2 + (HFM.P(:,2) - rec_y3).^2);
rec_x4 = 1296;
rec_y4 = 10000;
[rec_dist4, rec_id4] = min((HFM.P(:,1) - rec_x4).^2 + (HFM.P(:,2) - rec_y4).^2);
var.receiver = [rec_id1,rec_id2,rec_id3,rec_id4];

[~,recElemIDs,~ ] = intersect(HFM.elem(:,1),var.receiver');
HFM.recElemIDs =  recElemIDs;

%cd(oldfolder)

end


function bd_ids = getboundaryVertIDs(elem,boundary)
edge1 = [elem(:,1),elem(:,2)];
edge2 = [elem(:,2),elem(:,3)];
edge3 = [elem(:,3),elem(:,1)];
edge1BD = edge1(boundary(:,1),:);
edge2BD = edge2(boundary(:,1),:);
edge3BD = edge3(boundary(:,1),:);
v1 = intersect(edge1BD(:,1:2),boundary(:,2:3),'rows');
v2 = intersect(edge2BD(:,1:2),boundary(:,2:3),'rows');
v3 = intersect(edge3BD(:,1:2),boundary(:,2:3),'rows');

bd_ids = unique([v1;v2;v3]);
if isempty(bd_ids)
    edge1 = fliplr(edge1);
    edge2 = fliplr(edge2);
    edge3 = fliplr(edge3);
    edge1BD = edge1(boundary(:,1),:);
    edge2BD = edge2(boundary(:,1),:);
    edge3BD = edge3(boundary(:,1),:);
    v1 = intersect(edge1BD(:,1:2),boundary(:,2:3),'rows');
    v2 = intersect(edge2BD(:,1:2),boundary(:,2:3),'rows');
    v3 = intersect(edge3BD(:,1:2),boundary(:,2:3),'rows');
    bd_ids = unique([v1;v2;v3]);

end
end