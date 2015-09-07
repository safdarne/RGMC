function [tform, selectedFrom1, selectedFrom2, bestObjFunc,  status] ...
    = estimateTransform(matched_points1, matched_points2, ...
    imgA, imgB, motionHist, myFlag, prevTrans, maxNumTrials)


points1 = matched_points1.Location;
points2 = matched_points2.Location;
status  = 0;

[isFound, tmatrix, inliers, selectedFrom1, selectedFrom2, bestObjFunc] = msac(maxNumTrials, ...
         points1, points2, imgA, imgB, motionHist, myFlag, prevTrans);

tform = projective2d(tmatrix);

%==========================================================================
function [isFound, tform, inliers, selectedFrom1, selectedFrom2, bestObjFunc] = msac(maxNumTrials, ...
         points1, points2, imgA, imgB, motionHist, myFlag, prevTrans)
%%
numPts = size(points1, 1);
idxTrial = 1;
numTrials = int32(maxNumTrials);
bestObjFunc = 10e8;
bestTempDis = 10e8;

bestTForm = [];
selectedFrom1 = [];
selectedFrom2 = [];

maxPossible = 1;
for i = 0:3
    maxPossible = maxPossible * (numPts-i);
end
maxPossible = maxPossible / factorial(4);

addpath ('.\vgg_multiview\')
[s_p,ang_p,t_p,~] = convertTformToSRT(prevTrans.T);

se = strel('disk',2,4);
imgAedge = edge(imgA);
imgBedge = edge(imgB);

minX = min(points1(:,2));
maxX = max(points1(:,2));
minY = min(points1(:,1));
maxY = max(points1(:,1));
expX = maxX - minX;
expY = maxY - minY;
bestTForm = zeros(3,3);
while idxTrial <= min(numTrials, maxPossible)
    %% Make sure we do not select nearby points. This increases accuracy of homography estimation
    indices0 = randperm(numPts, 4);
    distantPoints = ones(size(points1,1),1);
    
    for i=1:4
        tempInd = find(distantPoints==1);
        if (length(tempInd)>0)
            ind = randperm(length(tempInd), 1);
            indices(i) = tempInd(ind);
            distantPoints(find(abs(points1(:,2)-points1(tempInd(ind),2)) < expX * 0.2 & abs(points1(:,1)-points1(tempInd(ind),1)) < expY * 0.2)) = 0;
        else
            indices(i) = indices0(i);
        end            
    end

    %%
    H = vgg_H_from_x_lin(points1(indices,:)',points2(indices,:)'); %1ms   
    tform = H';
    T = projective2d(tform);
    notCollinear = 1;        
    entr = pointsEntropy(points1(indices,:), size(imgA,1), size(imgA,2));%1ms     
    if (entr > 1.5) % Make sure the points are not on a line
        if (notCollinear)
            idxPt = 1;        
            b = 0.9;
            ND = numPts;
            threshInlier = b*(log(size(imgA,1)^2*size(imgA,1)^2)-log(16*b^4)); 

            inliers = 0;outliers = 0;
            %%            
            accDis = 0;
            disArr = points2 - transformPointsForward(T, points1(:, :));            
            disArr = abs(disArr(:,1)) + abs(disArr(:,2)); 
            inliers = length(find(disArr < threshInlier));
            outliers = numPts - inliers;
            disArr(find(disArr >= threshInlier)) = threshInlier;
            accDis = sum(disArr);
            
            if (myFlag) % if in initial stage of checking each cluster, the keypoint are more consistent, so less tolerance
                inlierBound = 0.7;
            else
                inlierBound = 0.3;
            end            
            if (inliers / (inliers+outliers) > inlierBound)
                if (myFlag) % if in initial stage of checking each cluster, the keypoint are more consistent, so less tolerance
                    errorBound = bestTempDis * 1.2;
                else
                    errorBound = bestTempDis * 3;
                end
                if (accDis < errorBound)   
                    bestTempDis = accDis;
                    [s,ang,t,~] = convertTformToSRT(tform);
                    
                    %% Do imwarp and ECC on edges directrly                       
                    [edgeX, edgeY] = find(imgBedge == 1);
                   
                    myedge = zeros(size(imgB));
                    edgeIndices = zeros(2,size(edgeX,1));               

                    [edgeIndices(2,:),edgeIndices(1,:)]=(transformPointsForward(T,edgeY,edgeX));                
                    edgeIndices = round(edgeIndices);

                    edgeIndices = edgeIndices(:,find((edgeIndices(1,:) > 0)&(edgeIndices(2,:) > 0)));                
                    edgeIndices(1,find(edgeIndices(1,:) > size(imgB,1))) = size(imgB,1);
                    edgeIndices(2,find(edgeIndices(2,:) > size(imgB,2))) = size(imgB,2);                
                    myedge(sub2ind(size(imgB), edgeIndices(1,:), edgeIndices(2,:))) = 1;

                    c = 0.001;
                    f = imgAedge .* (1 - motionHist);                   
                    g = myedge .* (1 - motionHist);      
                    ecc = 2 * sum(f(:) .* g(:)) / (sum(f(:)) + sum(g(:)) + c);   
                    
                    %% .52 is the original
                    objFunc = -log(gaussmf(ecc,[.08 .52])/(sqrt(2*pi)*.04)) ...
                        +accDis/numPts*20 ... 
                        +min(100,-log(gaussmf(abs(s-s_p),[2e-3 0])/(sqrt(2*pi)*2e-3)))...
                        +min(100,-log(gaussmf(abs(ang-ang_p)*180/pi,[2e-1 0])/(sqrt(2*pi)*2e-1)))...
                        +min(100,-log(1/(2*3.5)*exp(-abs(t(1)-t_p(1))/3.5)))...
                        +min(100,-log(1/(2*2.5)*exp(-abs(t(2)-t_p(2))/2.5)));

                    if objFunc < bestObjFunc
                        bestObjFunc = objFunc;
                        bestTForm = tform;
                        selectedFrom1 = indices;
                        selectedFrom2 = indices;                        
                        [-log(gaussmf(ecc,[0.1 0.45])),accDis, ecc, bestObjFunc];
                    end
                end
            end
        end
    end
    idxTrial = idxTrial + 1;
end

tform = bestTForm;
isFound = 1;
inliers = 1;

