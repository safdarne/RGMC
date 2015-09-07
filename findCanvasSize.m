function [minX, minY, maxX, maxY] = findCanvasSize(tforms, imageSize)
    clear Hcumulative
    minX = 10e5; minY = 10e5; maxX = -10e5; maxY = -10e5;
    Hcumulative2 = eye(3);
    for j=1:size(tforms,2)
        if (~isempty(tforms{j}))
            Hcumulative2 = tforms{j}  *  Hcumulative2; 
            [x, y] = transformPointsForward(projective2d(Hcumulative2),[1, imageSize(2), imageSize(2), 1]',[1, 1, imageSize(1), imageSize(1)]');
            minX = min([minX; x]);
            minY = min([minY; y]);
            maxX = max([maxX; x]);
            maxY = max([maxY; y]);
        end                
    end     

    minX = round(max(-3000, minX));
    minY = round(max(minY, -imageSize(1) * 0.5));
    maxX = round(min(3000, maxX));
    maxY = round(min(maxY, (1 + 0.5) * imageSize(1)));