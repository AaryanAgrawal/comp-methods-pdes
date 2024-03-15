% Finction to plot field defined on a 2D triangular mesh
% Matt McGarry, 2/11/2022
% Input: x = nnx1 vector of the x position of nodes (One column of the .nod file
%        x = nnx1 vector of the x position of nodes
%        in = nex3 nodal incidence list (get from the .elm file)
%        v = nnx1 vector of nodal field values;
%        Vec(optional) = Vector field: can be nnx2, defined at the nodes, or nex2, defined at the element centroids.  

% Example:
% junk=load('hw44.nod');x=junk(:,2);y=junk(:,3);
% junk=load('hw44.ele');in=junk(:,2:4);
% V=x.^2-0.1*y.^3; % Make up field to plot;
% Vec_nodal=rand(size(x,1),2);Vec_centroid=rand(size(in,1),2); % Make up vector field
% [h]=Plot2dTriMesh(x,y,in,V);
% [h]=Plot2dTriMesh(x,y,in,V,Vec_nodal);
% [h]=Plot2dTriMesh(x,y,in,V,Vec_centroid);

% Notes for 105 students: You can also use scatteredinterpolant, but that
% does not have the option of using a defined triangularization, it
% automatically uses a Delauney traingularization so the elements may be
% different from the ones you used for the FEM problem. 
% Trisurf plots only the elements you tell it to. 
function Plot2dTriMesh(x,y,in,v)

nlevels=10; % Number of color levels to plot

trisurf(in,x,y,zeros(size(x)),v,'facecolor','interp','linestyle',':'); % Plot flat surface colored by v
hold on 
% ANother option:
% trisurf(in,x,y,v,v,'facecolor','interp','linestyle',':'); % Plot surface height and color as v
colormap(parula)
view(0,90); % look straight at the plot
xlabel('x');
ylabel('y');
colorbar;

axis equal tight

end