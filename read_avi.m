function [mov] = read_avi(str) 
xyloObj = VideoReader([str,'.avi']);
% 
%nFrames = xyloObj.NumberOfFrames;
nFrames=20;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;


% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'));

% Read one frame at a time.
for k = 1 : nFrames
    mov(k).cdata = read(xyloObj, k);
end

% % Size a figure based on the video's width and height.
% hf = figure;
% set(hf, 'position', [150 150 vidWidth vidHeight])
% 
% % Play back the movie once at the video's frame rate.
% movie(hf, mov, 1, xyloObj.FrameRate);
% close all;
end