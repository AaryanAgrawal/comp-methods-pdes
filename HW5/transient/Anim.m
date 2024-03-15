clear; close all

u = readmatrix("transient.dat")';
u = u(2:end, :)

nodes = readmatrix("npeltr4.dat");
xNodes = nodes(:, 2);
yNodes = nodes(:, 3);
elemsJunk = readmatrix("epeltr4.dat");
elems = elemsJunk(:, 2:6);

Vec = zeros(size(elems,1),2);

%Triangulate existing Mesh
[elemsTri, elemsTriMat] = meshProb5(elems);

%%Find tumor elements
tumorElem = elemsTri(elemsTriMat==8, :);
tumorNode = unique(tumorElem);
tumorxNodes = xNodes(tumorNode);
tumoryNodes = yNodes(tumorNode);
% tumoruI = uI(tumorNode);
% tumorVec = zeros(size(tumorElem, 1), 2);

figure;

for t = 1:25:2499
    trisurf(elemsTri,xNodes,yNodes,zeros(size(xNodes)),u(:, t),'facecolor','interp','linestyle',':'); % Plot flat surface colored by v
    hold on 
    colormap(parula)
    view(0,90); % look straight at the plot
    xlabel('x');
    ylabel('y');
    axis equal tight
    k = convhull(tumorxNodes, tumoryNodes);
    plot(tumorxNodes(k), tumoryNodes(k), 'k', 'LineWidth',1);
    xlim([-0.25 0.25])
    ylim([-0.15 0.15])
    %clim([ min(u, 2) max(u, 2) ]);
    c = colorbar;
    c.Label.String = 'Temp (C)';
    title("T - Baseline Case")
    legend('', 'Tumor', 'FontSize', 6)
    hold off
    pause(0.001);
    clf;
end