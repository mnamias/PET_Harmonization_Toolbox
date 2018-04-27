function [image, dicom_headers] = read_directory()
% function [image, dicom_headers] = read_directory()
% Copyright 2018 Mauro Namías - mnamias@gmail.com
% reads a directory containing dicom files of the same series, sorts them
% by slice position and returns a volume (image) and the corresponding
% dicom headers in a structure.

%% Select and scan input directory
%% Select and scan input directory
if ismac
    % Code to run on Mac plaform
    dire = [uigetdir '/'];
elseif isunix
    % Code to run on Linux plaform
    dire = [uigetdir '/'];
elseif ispc
    % Code to run on Windows platform
    dire = [uigetdir '\'];
else
    disp('Platform not supported')
end

D_PET = dir(dire);

% skips non-archive entries
i = 1;
while(D_PET(i).isdir == 1)
    i = i+1;
end

name = (D_PET(i).name);


%% Reads all the files that belong to the same series (same series UID)

D_PET = dir(dire);
good_flags = zeros(1,length(i+1:length(D_PET)));

%% Initialize dicom_headers
dicom_headers_temp = dicominfo([dire D_PET(i).name]);

for index = i:length(D_PET)
    dicom_headers(index).dicom_headers = dicom_headers_temp;
end

vUID = dicom_headers(i).dicom_headers.SeriesInstanceUID;  % Valid series UID

%% Read dicom headers

tic
parfor index = i:length(D_PET)
dicom_headers(index).dicom_headers = dicominfo([dire D_PET(index).name]);
    disp(['reading file: ' int2str(index) ' / ' int2str(length(D_PET)-2) ]); 
end
toc
%%

%% Checks for valid UIDs
for index = i:length(D_PET)   

[sa,sb] =  size(fieldnames(dicom_headers(index).dicom_headers));
[saa,sbb] = size(fieldnames(dicom_headers(i).dicom_headers));
        
if (strcmp(dicom_headers(index).dicom_headers.SeriesInstanceUID,vUID) )
    good_flags(index) = 1;
    disp(['Including file: ' int2str(index) ' / ' int2str(length(D_PET)-2) ]);
end
    
    
end

dicom_headers_ok = dicom_headers(good_flags==1);
clear dicom_headers;
dicom_headers = dicom_headers_ok;
clear dicom_headers_ok;

%% Reads the image files, applies rescale slope and intercept
image = single(zeros(dicom_headers(1).dicom_headers.Rows,dicom_headers(1).dicom_headers.Columns,length(dicom_headers)));
clear di;
for i = 1:length(dicom_headers)
   i
    dicom_headers(i).dicom_headers.SmallestImagePixelValue = 0;
                    dicom_headers(i).dicom_headers.LargestImagePixelValue = 32767;
   dicom_headers(i).dicom_headers = orderfields(dicom_headers(i).dicom_headers,dicom_headers(1).dicom_headers); %%

    dummy = dicom_headers(i).dicom_headers;
di(i) = dummy;

if(isfield(di(i),'RescaleSlope'))
else
di(i).RescaleSlope = 1;
end


if(isfield(di(i),'RescaleIntercept'))
else
di(i).RescaleIntercept = 0;
end



end

clear dicom_headers;
dicom_headers = di;
clear di;


%% read valid images
tic
for index = 1:length(dicom_headers)
    image_d = dicomread(dicom_headers(index).Filename);
    included_files(index).name = dicom_headers(index).Filename;
    disp(['Reading image: ' int2str(index)]);
    image(:,:,index) = single(single(image_d)*dicom_headers(index).RescaleSlope + dicom_headers(index).RescaleIntercept);    
end
toc

%%

[a,b,c] = size(image);
if (c==1)
delete(included_files(1).name);
return
end


%% sorts the images in consecutive slice order

dicom_headers_back = dicom_headers;

ii = 0;
for index_t = 1:length(included_files)
        ii = ii+1;
   dicom_headers(ii).pixels = image(:,:,ii); 
   dicom_headers(ii).file_name = included_files(ii).name;
end


[dicom_headers_sorted index] = nestedSortStruct(dicom_headers, 'InstanceNumber', 1);

ii = 0;
for index_t = 1:length(included_files)
     ii = ii+1;
   image(:,:,ii) = dicom_headers_sorted(ii).pixels;
   included_files(ii).name = dicom_headers(ii).file_name;
end


dicom_headers_sorted = rmfield(dicom_headers_sorted,'pixels'); %remove pixel fields
dicom_headers_sorted = rmfield(dicom_headers_sorted,'file_name'); %remove file name field

clear dicom_headers;

dicom_headers = dicom_headers_sorted;
clear dicom_headers_sorted;