objArr = [];
tforms = {};
allFilesError = {};
files = dir(datapath);
totalFramesProcd = 0;

disp (fileName)
saveName = fileName(1:strfind(fileName,'.')-1);
filename = strcat(datapath,fileName);
vidObj = VideoReader(filename);    

% Read the first frame to set the output canvas size
% This canvas is just for temporary display, for the final video
% written to the final the canvas will be calculated automatically
imgB = read(vidObj,1);
imgB = imresize(imgB, resizeFactor);
tempref = imref2d(([round((posSpanY-negSpanY) * size(imgB,1)),round((posSpanX-negSpanX) * size(imgB,2))]));
tempref.XWorldLimits = [negSpanX * size(imgB,2) posSpanX * size(imgB,2)];
tempref.YWorldLimits = [negSpanY * size(imgB,1) posSpanY * size(imgB,1)];                

% Initial homography sets the first frame in the middle of the canvas
Hcumulative = eye(3);
Hcumulative(3,1) = size(imgB,2)  *  -0.5;
Hcumulative(3,2) = size(imgB,1)  *  -0.5;

overlaidIm = imwarp(zeros(size(imgB)), projective2d(Hcumulative), 'OutputView', tempref);
minStartX = size(overlaidIm,2);
maxStartX = 1;
minStartY = size(overlaidIm,2);
maxStartY = 1;
se = strel('disk',5);                       

ii = startFrame;
imgB = double(read(vidObj,ii)) / 255;
imgB = imresize(imgB, resizeFactor);
M = zeros(size(imgB,1),size(imgB,2));    
mytform.T = eye(3);   
prevTrans.T = eye(3);   
totalFramesProcd = 0;