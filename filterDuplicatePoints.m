function [pointsA, pointsB] = filterDuplicatePoints(pointsA, pointsB)
    temp = pointsA.Location;
    [A,iA,icA] = unique(temp(:,1));
    temp = pointsB.Location;
    [B,iB,icB] = unique(temp(:,1));
    ind = intersect(iA,iB);
    pointsA = pointsA(ind, :);
    pointsB = pointsB(ind, :);