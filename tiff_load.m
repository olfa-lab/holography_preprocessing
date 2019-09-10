function [ Figs,num_green_images,name,path ] = tiff_load(name,path,redchannel)
%tiff_load Load in a stack of images using Matlab tiff class.
%   This function gets a stack of alternating tiff figures of the red and green
%   channels and output a matrix containing all the figures.
%   Jon Gill 2019

if nargin < 2
    [name,path]=uigetfile('*.tif','Please choose the Tiff Stack file','E:\StimTrials');
end

if nargin < 3
    redchannel = 0;
end

fname=fullfile(path,name);
info = imfinfo(fname);
num_images = numel(info);
num_green_images=num_images/(redchannel+1);
Figs=zeros(info(1).Width,info(1).Height,num_green_images);
n=1;
TifLink = Tiff(fname, 'r');
for k = 1:redchannel+1:num_images
    TifLink.setDirectory(n);
    Figs(:,:,n) = TifLink.read();
    n=n+1;
end
TifLink.close();
end
