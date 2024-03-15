function [elemsTri, elemsTriMat] = meshProb5(elems)

elemsTri = elems(:, 1:3);
elemsTriMat = elems(:, 5);
elemsQuadIndex = find(elems(:, 3) ~= elems(:, 4));

for i = elemsQuadIndex
    elemsRow = [elems(i, 1), elems(i, 3), elems(i, 4)];
    elemsTri = [elemsTri; elemsRow];
    elemsTriMat = [elemsTriMat; elemsTriMat(i)];
end


end