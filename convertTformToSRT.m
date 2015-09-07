function [s,ang,t,R] = convertTformToSRT(H)
%Convert a 3-by-3 affine transform to a scale-rotation-translation
%transform.
%  [H,S,ANG,T,R] = cvexTformToSRT(H) returns the scale, rotation, and
%  translation parameters, and the reconstituted transform H.

% Extract rotation and translation submatrices
R = H(1:2,1:2);
t = H(3, 1:2);
% Compute theta from mean of stable arctangents
ang = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
% Compute scale from mean of two stable mean calculations
s = mean(R([1 4])/cos(ang));
% Reconstitute transform
R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
