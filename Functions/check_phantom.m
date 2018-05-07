function results = check_phantom(image, dicom_headers, showfig)
%% Name of code: check_phantom.m
% Version: 1.0 (April 27, 2018)
% Level of code: function. This function is from the PET harmonization toolbox GUI
% Copyright (c) 2018, Mauro Namías - mnamias@gmail.com
% All rights reserved.
% This code is intended to accompany the paper:
% Namías et. al
% A novel approach to quantitative harmonization in PET
% PMB (2018)
%
% function results = check_phantom(image, dicom_headers, showfig)
% Performs an alignment QC of the cylindrical phantom images 
%
% Inputs: 
%   im3D: the phantom slices (i.e.: 192x192x95)
%   dicom_headers: structure with the dicom headers for each slice
%   showfig: if 1, plot figures. 
%
% Outputs:
% results: MATLAB structure with the offsets (mm) and tilts (degrees) of the phantom
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% Variable initialization
slice_thk = dicom_headers(1).SliceThickness;
slice_spacing = abs( dicom_headers(1).ImagePositionPatient(3) - dicom_headers(2).ImagePositionPatient(3) );

pixel_sizes = [dicom_headers(1).PixelSpacing' min(slice_thk, slice_spacing) ];
results = [];

%% find z_profile
image = double(image);
z_profile = sum(image,1);
z_profile = sum(z_profile,2);
z_profile = squeeze(z_profile);


%% Determine if the image should be flipped in the axial direction or not
z_profile = squeeze(sum(sum(image,1),2));

bkg = z_profile - z_profile;
bkg( z_profile < max(z_profile(:)) * 0.1 )= 1;

bkg_indexes = find(bkg);
bkg_indexes2 = bkg_indexes(2:end)-bkg_indexes(1:end-1);
slice = find(bkg_indexes2-1);

f1 = sum( bkg( 1:slice - 1 ));
f2 = sum( bkg( slice + 1:end ));



    % not
    z_profile = squeeze(sum(sum(image,1),2));
    bkg = z_profile - z_profile;
    bkg( z_profile < max(z_profile(:)) * 0.1 )= 1;
    bkg_indexes = find(bkg);
    
    bkg_indexes2 = bkg_indexes(2:end)-bkg_indexes(1:end-1);
    slice = find(bkg_indexes2-1); % phantom surrounded by air on both sides (wide axial FOV scanners)
          
    if(not(isempty(slice)))

        f1 = sum( bkg( 1:slice - 1 ));
        f2 = sum( bkg( slice + 1:end ));
        
        if(f2 > f1)
            image = flipdim(image, 3);
            clc
            disp('Flipping image')
            z_profile = squeeze(sum(sum(image,1),2));          
        end
    end
    
        if( z_profile(1) > z_profile(end) & isempty(slice))
            image = flipdim(image, 3);
            clc
            disp('Flipping image')
            z_profile = squeeze(sum(sum(image,1),2));

      
        end

%%


%% Verify phantom alignment and tilt
z_profile = medfilt1(double(z_profile),5);
z_analysis_range = find( z_profile > max(z_profile(:)) * 0.9);
z_analysis_range = z_analysis_range(2:end-1);

mask = fspecial('disk',size(image,1)/2);
mask = mask(1:size(image,1),1:size(image,2));
mask = mask./max(mask(:));
mask(1,:) = 0;
mask(:,1) = 0;
mask(end,:) = 0;
mask(:,end) = 0;

se = strel('disk',9);
mask = imerode(mask,se);

for i = min(z_analysis_range(:)):max(z_analysis_range(:))
    slice = image(:,:,i);
    h = fspecial('gaussian',9,5);
    filtered_slice = imfilter(slice, h) .* mask;

    s = regionprops( filtered_slice >  max(filtered_slice(:))*0.2  , filtered_slice, {'Centroid','WeightedCentroid'});
    centroids(i - min(z_analysis_range(:)) + 1 , 1) = s(1).WeightedCentroid(1);
    centroids(i - min(z_analysis_range(:)) + 1 , 2) = s(1).WeightedCentroid(2);    
    
end

FOV_center = [size(image,1)/2 size(image,1)/2] ;

centroids = bsxfun(@minus, centroids, FOV_center);

x_centroids = centroids(:,1);
y_centroids = centroids(:,2);

x_centroids = smooth(x_centroids,5) * pixel_sizes(1);
y_centroids = smooth(y_centroids,5) * pixel_sizes(2);

results.x_offset = mean(x_centroids);
results.y_offset = mean(y_centroids);

z_positions = (min(z_analysis_range(:)):max(z_analysis_range(:))) * pixel_sizes(3);

p = polyfit(z_positions', x_centroids, 1);



x_tilt = atan(p(1))*180/pi;

p = polyfit(z_positions', y_centroids, 1);
y_tilt = atan(p(1))*180/pi;

if(showfig)
figure
subplot(2,1,1)
plot(z_analysis_range * pixel_sizes(3), x_centroids)
xlabel('Axial Position [mm]')
ylabel(' X-offset [mm]' )
title(['Phantom centering and alignment. X-tilt = ' num2str(x_tilt) ' degrees'])

subplot(2,1,2)
plot(z_analysis_range*pixel_sizes(3), y_centroids)
xlabel('Axial Position [mm]')
ylabel(' Y-offset [mm]' )
title(['Phantom centering and alignment. Y-tilt = ' num2str(y_tilt) ' degrees'])
end

results.x_tilt = x_tilt;
results.y_tilt = y_tilt;

%% Verify top edge of the phantom positioning 
acq_times_full = {dicom_headers.AcquisitionTime};
for i = 1:length(acq_times_full)
    temp = cell2mat(acq_times_full(i));
    acq_matrix(i) = mat2cell(temp(1:6),1,6);
end

if(strcmp( dicom_headers(1).Manufacturer, 'SIEMENS'))
    acq_times = smooth(datenum((acq_matrix),'hhmmss'),17);
    acq_times_diff = diff(acq_times);
    acq_times_diff(acq_times_diff < 1e-3) = 0;
    [pks,locs] = findpeaks(acq_times_diff,'MINPEAKDISTANCE', 10)  
    number_beds = length(locs)+1;
    c = size(image,3);

    acq_times = datenum((acq_matrix),'hhmmss');
    acq_times_diff = diff(acq_times);
    acq_times_diff(acq_times_diff < 1e-3) = 0;
    acq_times_diff(acq_times_diff > 1e-3) = 1;
    nslices = find( acq_times_diff);
    n_slices = nslices(round( length(nslices)/(number_beds-1)));

else
[acq_times, ia, ic] = unique(acq_matrix);
%acq_times_full = cell2mat(acq_times_full);
number_beds = length(acq_times);
    n_slices = sum(sum(ic==number_beds));
    
end


if(number_beds > 1)
    c = size(image,3);
    overlap = round((n_slices * number_beds - c) / (number_beds - 1));

    overlap_start = n_slices-overlap;
    overlap_end = n_slices;
    ideal_position = (overlap_end + overlap_start)/2;

    diff_z_profile = diff(z_profile);
    [max_value max_index] = max(diff_z_profile);

    z_offset = (max_index-ideal_position) * pixel_sizes(3);
    results.z_offset = z_offset;
else
    results.z_offset = [];
end



