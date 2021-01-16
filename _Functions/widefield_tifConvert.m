function [ new_filename ] = widefield_tifConvert(pn, fn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
disp('load first file')
%[filename,pathname]=uigetfile('.tif');
%underscore_index=find(filename=='_');
%new_filename=[filename(1:underscore_index(end)-1) '.tif'];
%number_idx=[find(filename=='.')-5:find(filename=='.')-1];
if nargin == 0
[fn,pn] = uigetfile('.tif');
end

firstfile = fn;
% find indices and initialize
underscore_index = find(firstfile=='_');
new_filename = [firstfile(1:underscore_index(1)-1) '.tif'];
number_idx = [find(firstfile=='.')-5:find(firstfile=='.')-1];
%number_idx = [find(firstfile=='.')-5:find(firstfile=='.')-1];

numFrames = 0;
stop_flag = 0;
filename = firstfile;
while stop_flag==0
    if exist(filename,'file')==2
        numFrames = numFrames+1;
        filename(number_idx)=num2str(str2num(filename(number_idx))+1,'%05d');
        %filename(number_idx)=num2str(str2num(filename(number_idx))+1,'%05d');
    elseif exist(filename,'file')==0
        stop_flag = 1;
    end
end

filename = firstfile;
header = imfinfo(filename);
xPixels = header(1).Width;
yPixels = header(1).Height;
image_matrix = zeros([yPixels xPixels numFrames],'uint16');


disp('Reading single-page Tif files...')
for i = 1:numFrames
    image_matrix(:,:,i) = imread(filename);
    filename(number_idx)=num2str(str2num(filename(number_idx))+1,'%05d');
    %filename(number_idx)=num2str(str2num(filename(number_idx))+1,'%05d');

    if rem(i,10)==0
        subroutine_progressbar(i/numFrames);
    end
end
subroutine_progressbar(1); 
close all

% write to tif file (tif files >4GB supported)
disp('Writing to multi-page Tif file...')
cd(pn)
options.big = true;
options.overwrite = true;
subroutine_saveastiff(image_matrix, new_filename, options);  
   

% delete single files
disp('Deleting single-page Tif files...')
try imread(new_filename,numFrames);
    filename = firstfile;
    for i = 1:numFrames
        delete(filename)
        filename(number_idx)=num2str(str2num(filename(number_idx))+1,'%05d');
        %filename(number_idx)=num2str(str2num(filename(number_idx))+1,'%05d');

        if rem(i,10)==0
            subroutine_progressbar(i/numFrames);
        end
    end
    subroutine_progressbar(1);
    close all
    disp('Conversion complete.')
catch
    disp('Error: new file appears not to have saved properly, single-page Tif files preserved')
end


