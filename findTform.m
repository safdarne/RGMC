function [tform, diffImg, finalObj] = findTform(imgA, imgB, motionHistory, prevTrans, displayFlag, T_C, T_M)   
    tau = 0.3;   
    C = 50;
    K = 8;
    tau_s = 50;
    pointsA = detectSURFFeatures(imgA,'MetricThreshold',tau_s);
    pointsB = detectSURFFeatures(imgB,'MetricThreshold',tau_s);

    [featuresA, pointsA] = extractFeatures(imgA, pointsA);
    [featuresB, pointsB] = extractFeatures(imgB, pointsB);
    try
        indexPairs = matchFeatures(featuresA, featuresB,'MatchThreshold',50);
    catch
        disp('error matching points')
    end
    pointsA = pointsA(indexPairs(:, 1), :);
    pointsB = pointsB(indexPairs(:, 2), :);

    % Make sure there are not alsmot duplicate points, detected at
    % different scales
    [pointsA, pointsB] = filterDuplicatePoints(pointsA, pointsB);
    
    motionHistory = motionHistory > tau;
    minSize = 4;
    se = strel('disk', minSize, minSize);
    motionHistory = imdilate(motionHistory, se);    
    staticPixels = 1 - motionHistory;  
    if (displayFlag)
        figure(222);imagesc(staticPixels);title('Static pixels')
    end
    
    % Make sure points are not located on moving part of the image
    [pointsA, pointsB] = filterMovingPoints(pointsA, pointsB);
    
    % Cluster and display motion vectors
    vec = pointsB.Location - pointsA.Location;
    k = min(K, size(vec, 1));
    idx = kmeans(vec,k,'EmptyAction','drop');
    showCluster(imgA, imgB, pointsA, pointsB, idx, k, displayFlag);    
    %% Analyse each cluster
    finalObj = 1000000;
    ind = [1 : k];
    c = 0.001;            
    bestObjFunc = 100000 * ones(1,k);
    selectedIndices = [];
    fisrtPlot = 1;
    if (k>1)
        for j = 1 : k
            i = ind(j);
            indices = (idx == i);        
            if (sum(indices)>3)
                    [tform, selectedFrom1, selectedFrom2, bestObjFunc(j)] = estimateTransform(...
                        pointsB(indices), pointsA(indices),imgA, imgB,motionHistory, 1, prevTrans, T_C);             

                    if (bestObjFunc(j) < finalObj)
                        finalTform = tform;
                        finalObj = bestObjFunc(i);
                        imgBp = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB)));  
                        diffImg = abs(imgA - imgBp);
                        bestCluster = j;                                   
                    end          
            else 
                % Not enough matches, just check the translation
                if (sum(indices) > 1)
                    transl = pointsA(indices).Location-pointsB(indices).Location;
                    newtform =  projective2d([1 0 0; 0 1 0; mean(transl(:,1)) mean(transl(:,2)) 1]);
                    imgBp = imwarp(imgB, newtform, 'OutputView', imref2d(size(imgB)));

                    pointsBm = [];
                    pointsAm = [];
                    bestObjFunc(i) = 9999;                         
                            
                    imgAedge = edge(imgA);
                    imgBedge = edge(imgB);                                                        

                    [edgeX, edgeY] = find(imgBedge == 1);
                    myedge = zeros(size(imgB,1), size(imgB,2));
                    edgeIndices = zeros(2,size(edgeX,1));               
                    [edgeIndices(2,:),edgeIndices(1,:)]=(transformPointsForward(newtform,edgeY,edgeX));                
                    edgeIndices = round(edgeIndices);
                    edgeIndices = edgeIndices(:,find((edgeIndices(1,:) > 0)&(edgeIndices(2,:) > 0)));                
                    edgeIndices(1,find(edgeIndices(1,:) > size(imgB,1))) = size(imgB,1);
                    edgeIndices(2,find(edgeIndices(2,:) > size(imgB,2))) = size(imgB,2);                
                    myedge(sub2ind(size(imgB), edgeIndices(1,:), edgeIndices(2,:))) = 1;

                    c = 0.001;                       
                    f = imgAedge.*(1-staticPixels);                   
                    g = myedge.*(1-staticPixels);      
                    ecc = 2*sum(f(:).*g(:))/(sum(f(:))+sum(g(:))+c);                               
                    % ECC > 0.5 is the best, force it to get the best score                    
                    bestObjFunc(i) = -log(gaussmf(min(ecc,.5),[.08 .52])/(sqrt(2*pi)*.04))+100;                     
                    if (bestObjFunc(i) < finalObj)
                        finalTform = newtform;
                        finalECC = bestObjFunc(i);
                        diffImg = abs(imgA - imgBp);
                        bestCluster = j;  
                    end                    
                end
            end
        end
    end

    indices = (idx == bestCluster);
    %% Merge background clusters
    if k > 1
        [~, i] = min(bestObjFunc);
        normObj = exp(-bestObjFunc);
        normObj = normObj/sum(normObj);
        [~,ttt]=sort(normObj,'descend');
        i = ttt;
    else
        i = 1;
    end

    clear eccFinal
    finalObj = bestObjFunc(i(1));

    bestTform = finalTform;
    %
    continueMerging = true;
    if (length(i) == 1) %Fine tune the tform
        indices = (idx == i);
        [tform, ~, ~, ~, ~, ObjFunc] = etimateTransform(...
            pointsB(indices), pointsA(indices), imgA, imgB, motionHistory, 0, prevTrans, T_M);

                if (ObjFunc < finalObj)
                    finalTform = tform;
                    finalObj = ObjFunc;
                    imgBp = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB))); 
                    diffImg = abs(imgA - imgBp);
                end

    else    %Regularize the clusters and merge them
        for jj = 2 : length(i)
            if (continueMerging)
                indices = zeros(size(idx));
                for kk = 1 : jj % in a greedy fashion, add the best matched clusters
                    clusterindices = (idx == i(kk));
                    % if there are many sample, use at most 50 of those
                    if (sum(clusterindices) > C)
                        tempindices = find(clusterindices == 1);
                        clusterindices = zeros(length(clusterindices),1);
                        temprand = randperm(length(tempindices));
                        tempindices = tempindices(temprand(1:50));
                        clusterindices(tempindices) = 1;
                        clusterindices = logical(clusterindices);
                    end                
                    indices = indices | clusterindices;    
                end

                if (sum(indices) > 3)
                    bestObj = 10e6;
                    ObjFunc = [];

                    [bestTform, selectedFrom1, selectedFrom2, bestObj] = estimateTransform(...
                        pointsB(indices), pointsA(indices), imgA, imgB, staticPixels, 0, prevTrans, T_M);

                    if (bestObj < finalObj)
                        finalTform = bestTform;
                        finalObj = bestObj;
                        imgBp = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB))); 
                        diffImg = abs(imgA - imgBp);
                    else
                        continueMerging = false;
                        numberClustersMerged = kk - 1;
                    end
                    try
                        tform = finalTform;
                    catch
                        finalTform = projective2d(eye(3));
                        tform = finalTform;
                    end
                end
            end
        end
    end
    try
        tform= finalTform;
    catch
        disp('Error finding the transformation')
    end
    try
        disp([num2str(numberClustersMerged), ' clusters merged.'])
    end


