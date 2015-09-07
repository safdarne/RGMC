function showCluster(imgA, imgB, pointsA, pointsB,idx, k, displayFlag)
    theme = 'ymcrgkbwymcrgkbwymcrgkbw';    
    if (displayFlag)
        themeIndex = 1;
        for i = 1 : k
            indices = (idx == i);    
            figure(44); showMatchedFeatures(imgA, imgB, pointsA(indices,:), pointsB(indices,:),themeIndex==1,'PlotOptions', {'ro', 'g+', theme(themeIndex)});    
            title('Clusters of the matched keypoints')
            themeIndex = themeIndex + 1;
        end
    end