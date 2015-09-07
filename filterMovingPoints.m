function [pointsA, pointsB] = filterMovingPoints(pointsA, pointsB)
    
    ind = zeros(size(pointsA,1),1);
    for i=1:length(ind)
        try
            coord = pointsA(i).Location;
            if (staticPixels(round(coord(2)),round(coord(1))))
                ind(i) = 1;
            end
        end
    end
    if (sum(ind) > 30) % if enough keypoints will remain after this filtering, then apply it
        pointsA = pointsA(ind == 1);
        pointsB = pointsB(ind == 1);    
    end