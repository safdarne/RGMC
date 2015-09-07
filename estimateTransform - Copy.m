function [tform, inlier_points1, inlier_points2, selectedFrom1, selectedFrom2, bestObjFunc,  status] ...
    = my___estimateGeometricTransform(matched_points1, matched_points2, ...
    imgA, imgB, motionHist, myFlag, prevTrans, transform_type, varargin)

% List of status code
statusCode = struct(...
    'NoError',           int32(0),...
    'NotEnoughPts',      int32(1),...
    'NotEnoughInliers',  int32(2));

% Parse and check inputs
[points1,points2,sampleSize,maxNumTrials,confidence,threshold,status] ...
    = parseInputs(statusCode, matched_points1, matched_points2, ...
        transform_type, varargin{:});
sampleSize = 4;
classToUse = getClassToUse(points1, points2);

% Compute the geometric transformation
if status == statusCode.NoError
    [isFound, tmatrix, inliers, selectedFrom1, selectedFrom2, bestObjFunc] = msac(sampleSize, maxNumTrials, ...
        confidence, threshold, points1, points2, imgA, imgB, motionHist, myFlag, prevTrans);
    
        try
            tformtemp = projective2d(tmatrix);
        catch
            'error'
        end
    
    if ~isFound
        status = statusCode.NotEnoughInliers;
    end
else
    tmatrix = zeros([3,3], classToUse);
end

% Extract inlier points
if status == statusCode.NoError
    inlier_points1 = matched_points1(inliers, :);
    inlier_points2 = matched_points2(inliers, :);
else
    inlier_points1 = matched_points1([]);
    inlier_points2 = matched_points2([]);
    tmatrix = zeros([3,3], classToUse);
end

% Report runtime error if the status output is not requested
reportError = (nargout ~= 4);

if isTestingMode()
    % Return tform as a matrix for internal testing purposes
    tform = tmatrix;
else
    sampleSize = 4;%%%%
    if sampleSize < 4  % similarity or affine
        tform = affine2d(tmatrix);
    else               % projective
        
        tform = projective2d(tmatrix);

    end
end

%==========================================================================
% Parse and check inputs
%==========================================================================
function [points1, points2, sampleSize, maxNumTrials, confidence, ...
    maxDistance, status] ...
    = parseInputs(statusCode, matched_points1, matched_points2, ...
    transform_type, varargin)

isSimulationMode = isempty(coder.target);
if isSimulationMode
    % Instantiate an input parser
    parser = inputParser;
    parser.FunctionName = 'estimateGeometricTransform';
    parser.CaseSensitive = true;
    
    % Specify the optional parameters
    parser.addParamValue('MaxNumTrials', 1000);
    parser.addParamValue('Confidence',   99);
    parser.addParamValue('MaxDistance',  1.5);
    
    % Parse and check optional parameters
    parser.parse(varargin{:});
    r = parser.Results;
    
    maxNumTrials = r.MaxNumTrials;
    confidence   = r.Confidence;
    maxDistance  = r.MaxDistance;
    
else
    % Instantiate an input parser
    parms = struct( ...
        'MaxNumTrials',       uint32(0), ...
        'Confidence',         uint32(0), ...
        'MaxDistance',        uint32(0));
    
    popt = struct( ...
        'CaseSensitivity', true, ...
        'StructExpand',    true, ...
        'PartialMatching', false);
    
    % Specify the optional parameters
    optarg       = eml_parse_parameter_inputs(parms, popt,...
        varargin{:});
    maxNumTrials = eml_get_parameter_value(optarg.MaxNumTrials,...
        1000, varargin{:});
    confidence   = eml_get_parameter_value(optarg.Confidence,...
        99, varargin{:});
    maxDistance  = eml_get_parameter_value(optarg.MaxDistance,...
        1.5, varargin{:}); 
end

% Check required parameters
sampleSize = checkTransformType(transform_type);
points1 = checkAndConvertPoints(matched_points1);
points2 = checkAndConvertPoints(matched_points2);
status  = checkPointsSize(statusCode, sampleSize, points1, points2);

% Check optional parameters
checkMaxNumTrials(maxNumTrials);
checkConfidence(confidence);
checkMaxDistance(maxDistance);

%==========================================================================
function status = checkPointsSize(statusCode, sampleSize, points1, points2)

coder.internal.errorIf( size(points1,1) ~= size(points2,1), ...
    'vision:estimateGeometricTransform:numPtsMismatch');

coder.internal.errorIf( ~isequal(class(points1), class(points2)), ...
    'vision:estimateGeometricTransform:classPtsMismatch');

if size(points1,1) < sampleSize
    status = statusCode.NotEnoughPts;
else
    status = statusCode.NoError;
end

%==========================================================================
function points = checkAndConvertPoints(matched_points)
if isnumeric(matched_points)
    checkPointsAttributes(matched_points);
    points = matched_points;
elseif isa(matched_points, 'vision.internal.FeaturePoints')
    points = matched_points.Location;
else % MSERRegions
    points = matched_points.Centroid;
end

%==========================================================================
function checkPointsAttributes(value)
validateattributes(value, {'numeric'}, ...
    {'2d', 'nonsparse', 'real', 'size', [NaN, 2]},...
    'estimateGeometricTransform', 'MATCHED_POINTS');

%==========================================================================
function r = checkMaxNumTrials(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'integer', 'positive', 'finite'},...
    'estimateGeometricTransform', 'MaxNumTrials');
r = 1;

%========================================================================== 
function r = checkConfidence(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'positive', 'finite', '<', 100},...
    'estimateGeometricTransform', 'Confidence');
r = 1;

%==========================================================================
function r = checkMaxDistance(value)
validateattributes(value, {'numeric'}, ...
    {'scalar', 'nonsparse', 'real', 'positive', 'finite'},...
    'estimateGeometricTransform', 'MaxDistance');
r = 1;

%==========================================================================
function sampleSize = checkTransformType(value)
list = {'similarity', 'affine', 'projective'};
validatestring(value, list, 'estimateGeometricTransform', ...
    'TransformType');

switch(value(1))
    case 's'
        sampleSize = 2;
    case 'a'
        sampleSize = 3;
    otherwise
        sampleSize = 4;
end

%==========================================================================
function c = getClassToUse(points1, points2)
if isa(points1, 'double') || isa(points2, 'double')
    c = 'double';
else
    c = 'single';
end

%==========================================================================
function flag = isTestingMode
isSimulationMode = isempty(coder.target);
coder.extrinsic('vision.internal.testEstimateGeometricTransform');
if isSimulationMode
    flag = vision.internal.testEstimateGeometricTransform;
else
    flag = eml_const(vision.internal.testEstimateGeometricTransform);
end

%==========================================================================
% Algorithm for computing the fundamental matrix.
%==========================================================================
function T = computeSimilarity(points1, points2, classToUse)
numPts = size(points1, 1);
constraints = zeros(2*numPts, 5, classToUse);
constraints(1:2:2*numPts, :) = [-points1(:, 2), points1(:, 1), ...
    zeros(numPts, 1), -ones(numPts,1), points2(:,2)];
constraints(2:2:2*numPts, :) = [points1, ones(numPts,1), ...
    zeros(numPts, 1), -points2(:,1)];
[~, ~, V] = svd(constraints, 0);
h = V(:, end);
T = coder.nullcopy(zeros(3, classToUse));
T(:, 1:2) = [h(1:3), [-h(2); h(1); h(4)]] / h(5);
T(:, 3)   = [0; 0; 1];

%==========================================================================
function T = computeAffine(points1, points2, classToUse)
numPts = size(points1, 1);
constraints = zeros(2*numPts, 7, classToUse);
constraints(1:2:2*numPts, :) = [zeros(numPts, 3), -points1, ...
    -ones(numPts,1), points2(:,2)];
constraints(2:2:2*numPts, :) = [points1, ones(numPts,1), ...
    zeros(numPts, 3), -points2(:,1)];
[~, ~, V] = svd(constraints, 0);
h = V(:, end);
T = coder.nullcopy(zeros(3, classToUse));
T(:, 1:2) = reshape(h(1:6), [3,2]) / h(7);
T(:, 3)   = [0; 0; 1];

%==========================================================================
function T = computeProjective(points1, points2, classToUse)
numPts = size(points1, 1);
p1x = points1(:, 1);
p1y = points1(:, 2);
p2x = points2(:, 1);
p2y = points2(:, 2);
constraints = zeros(2*numPts, 9, classToUse);
constraints(1:2:2*numPts, :) = [zeros(numPts,3), -points1, ...
    -ones(numPts,1), p1x.*p2y, p1y.*p2y, p2y];
constraints(2:2:2*numPts, :) = [points1, ones(numPts,1), ...
    zeros(numPts,3), -p1x.*p2x, -p1y.*p2x, -p2x];
[~, ~, V] = svd(constraints, 0);
h = V(:, end);
T = reshape(h, [3,3]) / h(9);

%==========================================================================
function num = computeLoopNumber(sampleSize, confidence, pointNum, inlierNum)
num = int32(ceil(log10(1 - 0.01 * confidence)...
    / log10(1 - (inlierNum/pointNum)^sampleSize)));

%==========================================================================
function [ptsNorm, normMatrix] = normalizePts(pts, classToUse)
ptsNorm = cast(pts, classToUse);
cent = mean(ptsNorm, 1);
ptsNorm(:, 1) = ptsNorm(:, 1) - cent(1);
ptsNorm(:, 2) = ptsNorm(:, 2) - cent(2);

weight = std(ptsNorm(:),[],1);
if weight > 0
    weight = sqrt(2) / weight;
else
    weight = ones(1, classToUse);  % Just pick a value
end

ptsNorm(:, 1) = ptsNorm(:, 1) * weight;
ptsNorm(:, 2) = ptsNorm(:, 2) * weight;

normMatrix = [...
    1/weight,     0,            0;...
    0,            1/weight,     0;...
    cent(1),      cent(2),      1];

%==========================================================================
function tform = computeTForm(sampleSize, points1, points2, indices, classToUse)
[samples1, normMatrix1] = normalizePts(points1(indices, :), classToUse);
[samples2, normMatrix2] = normalizePts(points2(indices, :), classToUse);

switch(sampleSize)
    case 2
        tform = computeSimilarity(samples1, samples2, classToUse);
    case 3
        tform = computeAffine(samples1, samples2, classToUse);
    otherwise % 4
        tform = computeProjective(samples1, samples2, classToUse);
        tform = tform / tform(end);
end
tform = normMatrix1 \ (tform * normMatrix2);

%==========================================================================
function dis = evaluateTForm(sampleSize, threshold, tform, point1, point2, classToUse)
pt1 = cast(point1, classToUse);
pt2 = cast(point2, classToUse);
pt = pt1 * tform(1:2, 1:2) + tform(3, 1:2);
if sampleSize == 4
    denom = pt1 * tform(1:2, 3) + tform(3, 3);
    if abs(denom) > eps(classToUse)
        pt = pt ./ denom;
    else % Mark this point invalid by setting it to a location far away from pt2
        pt(:) = pt2 + threshold;
    end
end
dis = norm(pt - pt2);
dis = min(dis, threshold);

%==========================================================================
function [isFound, tform, inliers, selectedFrom1, selectedFrom2, bestObjFunc] = msac(sampleSize, maxNumTrials, ...
        confidence, threshold, points1, points2, imgA, imgB, motionHist, myFlag, prevTrans)
%%
numPts = size(points1, 1);
idxTrial = 1;
numTrials = int32(maxNumTrials);
bestObjFunc = 100000000;%maxDis;
bestTempDis = 1000000000;

bestTForm = [];
selectedFrom1 = [];
selectedFrom2 = [];

maxPossible = 1;
for i=0:sampleSize-1
    maxPossible = maxPossible * (numPts-i);
end
maxPossible = maxPossible/factorial(sampleSize);

% Create a random stream. It uses a fixed seed for the testing mode and a
% random seed for other mode.
if isTestingMode()
    rng('default');
end


addpath ('.\vgg_multiview\')
[s_p,ang_p,t_p,~] = convertTformToSRT(prevTrans.T);

sampleSize = 4;
se = strel('disk',2,4);
imgAedge = edge(imgA);
imgBedge = edge(imgB);
% imgAedge = imdilate(edge(imgA),se);
% imgBedge = imdilate(edge(imgB),se);

minX = min(points1(:,2));
maxX = max(points1(:,2));
minY = min(points1(:,1));
maxY = max(points1(:,1));
expX = maxX-minX;
expY = maxY-minY;
bestTForm = zeros(3,3);
while idxTrial <= min(numTrials, maxPossible)
    %% Make sure we do not select nearby points. This increases accuracy of homography estimation
    indices0 = randperm(numPts, sampleSize);
    distantPoints = ones(size(points1,1),1);
    
    for i=1:sampleSize
        tempInd = find(distantPoints==1);
        if (length(tempInd)>0)
            ind = randperm(length(tempInd), 1);
            indices(i) = tempInd(ind);
            distantPoints(find(abs(points1(:,2)-points1(tempInd(ind),2)) < expX*.2 & abs(points1(:,1)-points1(tempInd(ind),1)) < expY*.2)) = 0;
        else
            indices(i) = indices0(i);
        end            
    end

%     figure(26);clf;plot(points1(:,1),points1(:,2),'.');hold on;plot(points1(distantPoints==1,1),points1(distantPoints==1,2),'co')
%     plot(points1(indices,1),points1(indices,2),'rs')
    %%
    H = vgg_H_from_x_lin(points1(indices,:)',points2(indices,:)'); %1ms   
    tform = H';
    T = projective2d(tform);
    notCollinear = 1;        
    entr = pointsEntropy(points1(indices,:), size(imgA,1), size(imgA,2));%1ms     
    if (entr > 1.5) % make sure the points are not on a line
        if (notCollinear)
            accDis = 0;
            idxPt = 1;        
            b = 0.9;
            ND = numPts;
            threshInlier = b*(log(size(imgA,1)^2*size(imgA,1)^2)-log(16*b^4)); 

            inliers = 0;outliers = 0;
            %%
            
            accDis = 0;
            disArr=points2-transformPointsForward(T,points1(:, :));            
            disArr = abs(disArr(:,1))+abs(disArr(:,2)); %////?????
            inliers = length(find(disArr<threshInlier));
            outliers = numPts - inliers;
            disArr(find(disArr>=threshInlier)) = threshInlier;
            accDis = sum(disArr);
            
            if (myFlag) % if in initial stage of checking each cluster, the keypoint are more consistent, so less tolerance
                inlierBound = 0.7;
            else
                inlierBound = 0.3;
            end            
            if (inliers/(inliers+outliers)>inlierBound)
                if (myFlag) % if in initial stage of checking each cluster, the keypoint are more consistent, so less tolerance
                    errorBound = bestTempDis*1.2;
                else
                    errorBound = bestTempDis*3;
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
                    f = imgAedge.*(1-motionHist);                   
                    g = myedge.*(1-motionHist);      
                    ecc = 2*sum(f(:).*g(:))/(sum(f(:))+sum(g(:))+c);   
                    
                    %% .52 is the original
                    objFunc = -log(gaussmf(ecc,[.08 .52])/(sqrt(2*pi)*.04)) ...
                        +accDis/numPts*20 ... 
                        +min(100,-log(gaussmf(abs(s-s_p),[2e-3 0])/(sqrt(2*pi)*2e-3)))...
                        +min(100,-log(gaussmf(abs(ang-ang_p)*180/pi,[2e-1 0])/(sqrt(2*pi)*2e-1)))...
                        +min(100,-log(1/(2*3.5)*exp(-abs(t(1)-t_p(1))/3.5)))...
                        +min(100,-log(1/(2*2.5)*exp(-abs(t(2)-t_p(2))/2.5)));
                        %-log(gaussmf(entr,[1 -1.2])); %1ms

%                     objFunc = -log(gaussmf(ecc,[.04 .52])/(sqrt(2*pi)*.04))+accDis/numPts*20 ... 
%                         +min(100,-log(gaussmf(abs(s-s_p),[2e-4 0])/(sqrt(2*pi)*2e-3)))...
%                         +min(100,-log(gaussmf(abs(ang-ang_p)*180/pi,[2e-3 0])/(sqrt(2*pi)*2e-1)))...
%                         +min(100,-log(1/(2*3.5)*exp(-abs(t(1)-t_p(1))/3.5)))...
%                         +min(100,-log(1/(2*.5)*exp(-abs(t(2)-t_p(2))/.5)));
%                         ;%-log(gaussmf(entr,[1 -1.2])); %1ms

                    if objFunc < bestObjFunc
                        bestObjFunc = objFunc;
                        bestTForm = tform;
                        inlierNum = numPts - bestObjFunc / threshold;
                        selectedFrom1 = indices;
                        selectedFrom2 = indices;                        
                        [-log(gaussmf(ecc,[.1 .45])),accDis, ecc, bestObjFunc];
                    end
                end
            end
        end
    end
    idxTrial = idxTrial + 1;
end

tform = bestTForm;
isFound = 1;
distances = 0;
inliers = 1;

