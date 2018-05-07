function [faxis, MTF, lsf, esf] = LSFcalc_Z(im3D, dicom_headers, options);
%% Name of code: LSFcalc_Z.m
% Version: 1.0 (April 20, 2018)
% Level of code: function. This function is from the PET harmonization toolbox GUI
% Copyright (c) 2018, Mauro Namías - mnamias@gmail.com
% All rights reserved.
% This code is intended to accompany the paper:
% Namías et. al
% A novel approach to quantitative harmonization in PET
% PMB (2018)
%
% function [faxis, MTF, lsf, esf] = LSFcalc_Z(im3D, dicom_headers, options);
%
% Estimates the MTF, LSF and ESF in the Z direction (axial direction) from a PET cylindrical phantom scan .
%
% Inputs: 
%   im3D: the phantom slices (i.e.: 192x192x95)
%   dicom_headers: structure with the dicom headers for each slice
%   options: MATLAB structure with processing options
%
% Options: processing options
%   options.delta_x : Number of sagittal slices to average
%   options.over_factor_Z : Oversampling factor for LSF_Z analysis
%   options.thr_Z : Threshold to make the ESF_Z symmetrical (0 to 0.5). If set to zero, the ESF remains unchanged
%   options.use_antialias_Z : Use an antialiasing filter in FFT domain (0 or 1)
%   options.use_window_Z : use a multiplicative window before computing FFT (1 or 0)
%   options.pad_Z: padding of the LSF in mm (default = 40 mm)
%
% Outputs:
% faxis: frequency axis of the estimated MTF (in mm-1)
% MTF: estimated axial Modulation Transfer Function (unitless)
% lsf: estimated Line Spread Function (unitless, area normalized to 1)
% esf: estimated Edge Spread Function (unitless)
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

clc;

fprintf('****************************************************************************\n');
fprintf('PET axial MTF & LSF Calculation\n\n')
fprintf('This code is intended to accompany the paper:\n'); 
fprintf('Namías et. al  "A\n')  
fprintf('"A novel approach to quantitative harmonization in PET". \n')
fprintf('M Namías et al 2018 Phys. Med. Biol. 63 095019 \n')
fprintf('****************************************************************************\n\n');

%% User selected options
pixelz = dicom_headers(1).SliceThickness;
delta_x = options.delta_x; % Number of sagittal slices to average
over_factor = options.over_factor_Z; % Oversampling factor for LSF_Z analysis
thr = options.thr_Z; % Threshold to make the ESF_Z symmetrical. If set to zero, the ESF remains unchanged
use_antialias = options.use_antialias_Z;
window = options.use_window_Z;
pad_mm = options.pad_Z;
manual_flag = strcmp(options.mode, 'manual');

%% Determine the axial analysis range and flip de image if necessary
mip = mipss(im3D);  % creates a MIP of the cylindrical phantom
z_profile = squeeze(sum(sum(im3D,1),2));

bkg = z_profile - z_profile;
bkg( z_profile < max(z_profile(:)) * 0.1 )= 1;

bkg_indexes = find(bkg);
bkg_indexes2 = bkg_indexes(2:end)-bkg_indexes(1:end-1);
slice = find(bkg_indexes2-1);

f1 = sum( bkg( 1:slice - 1 ));
f2 = sum( bkg( slice + 1:end ));


% Determine if the image should be flipped in the axial direction or
    % not
    z_profile = squeeze(sum(sum(im3D,1),2));
    bkg = z_profile - z_profile;
    bkg( z_profile < max(z_profile(:)) * 0.1 )= 1;
    bkg_indexes = find(bkg);
    
    bkg_indexes2 = bkg_indexes(2:end)-bkg_indexes(1:end-1);
    slice = find(bkg_indexes2-1); % phantom surrounded by air on both sides (wide axial FOV scanners)
          
    if(not(isempty(slice)))

        f1 = sum( bkg( 1:slice - 1 ));
        f2 = sum( bkg( slice + 1:end ));
        
        if(f2 > f1)
            im3D = flipdim(im3D, 3);
            disp('Flipping image')
        end
    end
    
        if( z_profile(1) > z_profile(end) & isempty(slice))
            im3D = flipdim(im3D, 3);
            disp('Flipping image')
        end


z_profile = squeeze(sum(sum(im3D,1),2));
dz_profile = diff(z_profile(1:end/2)); 
maximum_derivative_index = find( abs(dz_profile) == max(abs(dz_profile(:))) );

z_range = maximum_derivative_index - round(100/pixelz) : maximum_derivative_index + round(100/pixelz);
z_range(z_range>size(im3D,3)) = [];
z_range(z_range<1) = [];

im3D = im3D(:,:,z_range);  % Crop the image to the axial range of interest

%%  Start of analysis
smooth_f = 3;
rstep = pixelz/over_factor; %  Choose bin size in mm (i.e., effective pixel size) for oversampled edge spread function (ESF)
% Pixel sizes and image size
pixelx = dicom_headers(1).PixelSpacing(2);
pixely = dicom_headers(1).PixelSpacing(1);
pixelz = dicom_headers(1).SliceThickness;
[a,b,c] = size(im3D);

% Create an average coronal image
coronal_image = squeeze(im3D(:,round(end/2)-round(delta_x/2):round(end/2)+round(delta_x/2),:));

coronal_image = squeeze(sum(coronal_image,2));
im = coronal_image;
  
% Show averaged image 
   
    f2 = figure;
    imagesc(im)
    xlabel('pixel')
    ylabel('pixel')
    colormap gray
    hold on
    axis image;
    impixelinfo;
% Determine background and foreground values to properly normalized the
% data to range [0,1]
BGvalue = 0;

if(manual_flag)  % manual definition of ROIs
  %  htext1 = text(-100,-10,'Click and drag to create a rectangular ROI representing the background.')
   % htext2 = text(-100,0,'Double click the ROI when finished.')
    title( {'Click and drag to create a rectangular ROI representing the background.' 'Double click the ROI when finished.'})


    fprintf('Click and drag to create a rectangular ROI representing the background.  Double click the ROI when finished.\n\n');
    
    h = imrect;
    ptROI = wait(h);
    
    %delete(htext1)    
    %delete(htext2)
    
    BGxROI1 = round(ptROI(1));
    BGxROI2 = BGxROI1 + round(ptROI(3));
    BGyROI1 = round(ptROI(2));
    BGyROI2 = BGyROI1 + round(ptROI(4));
    
    % force selected ROI to be within the bounds of the image
    BGxROI1 = min(BGxROI1,size(im,2));
    BGxROI2 = min(BGxROI2,size(im,2));
    BGyROI1 = min(BGyROI1,size(im,1));
    BGyROI2 = min(BGyROI2,size(im,1));
    BGxROI1 = max(BGxROI1,1);
    BGxROI2 = max(BGxROI2,1);
    BGyROI1 = max(BGyROI1,1);
    BGyROI2 = max(BGyROI2,1);
    BGvalue = mean(mean(im(BGyROI1:BGyROI2,BGxROI1:BGxROI2)));
    
else % automatic definition of ROIs
     % find the background value automatically
    
    vertical_profile = squeeze(im(:,round(end*3/4)));
    y_range = find(vertical_profile > max(vertical_profile * 0.5));

    z_profile = squeeze(sum(sum(im3D,1),2));
    bkg = z_profile-z_profile;
    bkg(z_profile<max(z_profile(:))*0.1)= 1;
    x_range = find(bkg==1);
    bkg_indexes = find(bkg)
    
    BGyROI1 = y_range(1)+5;
    BGyROI2 = y_range(end) - 5;
    BGxROI1 = x_range(1) + 5;
    BGxROI2 = x_range(end) - 5;
    
    plot( [BGxROI1  BGxROI1], [ BGyROI1  BGyROI2] );
    plot( [BGxROI1  BGxROI2], [ BGyROI1  BGyROI1] );
    plot( [BGxROI2  BGxROI2], [ BGyROI1  BGyROI2] );
    plot( [BGxROI2  BGxROI1], [ BGyROI2  BGyROI2] );
    BGvalue = mean(mean(im(BGyROI1:BGyROI2,BGxROI1:BGxROI2)));

    
end

FGvalue = 0;

if(manual_flag) % manual definition of ROIs

    title( {'Click and drag to create a rectangular ROI representing the foreground.' 'Double click the ROI when finished.'})

%    htext1 = text(-100,-10,'Click and drag to create a rectangular ROI representing the foreground.')
%    htext2 = text(-100,0,'Double click the ROI when finished.')

    fprintf('Click and drag to create a rectangular ROI representing the foreground.  Double click the ROI when finished.\n\n');
    
    h = imrect;
    ptROI = wait(h);

 %   delete(htext1)    
 %   delete(htext2)

    
    FGxROI1 = round(ptROI(1));
    FGxROI2 = FGxROI1 + round(ptROI(3));
    FGyROI1 = round(ptROI(2));
    FGyROI2 = FGyROI1 + round(ptROI(4));
    
    % force selected ROI to be within the bounds of the image
    FGxROI1 = min(FGxROI1,size(im,2));
    FGxROI2 = min(FGxROI2,size(im,2));
    FGyROI1 = min(FGyROI1,size(im,1));
    FGyROI2 = min(FGyROI2,size(im,1));
    FGxROI1 = max(FGxROI1,1);
    FGxROI2 = max(FGxROI2,1);
    FGyROI1 = max(FGyROI1,1);
    FGyROI2 = max(FGyROI2,1);
    FGvalue = mean(mean(im(FGyROI1:FGyROI2,FGxROI1:FGxROI2)));
else   % auto definition of ROIs
    vertical_profile = squeeze(im(:,round(end*3/4)));
    y_range = find(vertical_profile > max(vertical_profile * 0.5));
    x_range = round(size(im,2)*3/4) - 5 : round(size(im,2)*3/4) + 5;
    FGyROI1 = y_range(1) + 5;
    FGyROI2 = y_range(end) - 5;
    FGxROI1 = x_range(1);
    FGxROI2 = x_range(end);
    
    plot( [FGxROI1  FGxROI1], [ FGyROI1  FGyROI2] );
    plot( [FGxROI1  FGxROI2], [ FGyROI1  FGyROI1] );
    plot( [FGxROI2  FGxROI2], [ FGyROI1  FGyROI2] );
    plot( [FGxROI2  FGxROI1], [ FGyROI2  FGyROI2] );
    FGvalue = mean(mean(im(FGyROI1:FGyROI2,FGxROI1:FGxROI2)));
    
end


if(manual_flag) % manual definition of ROIs
    
    title( {'Click and drag to create a rectangular ROI representing the vertical analysis range.' 'Double click the ROI when finished.'})

    fprintf('Click and drag to create a rectangular ROI representing the analysis range.  Double click the ROI when finished.\n\n');
    %htext1 = text(-100,-10,'Click and drag to create a rectangular ROI representing the analysis range.')
    %htext2 = text(-100,0,'Double click the ROI when finished.')

    h = imrect;
    ptROI = wait(h);
    
    %delete(htext1)
    %delete(htext2)
    
    FGxROI1 = round(ptROI(1));
    FGxROI2 = FGxROI1 + round(ptROI(3));
    FGyROI1 = round(ptROI(2));
    FGyROI2 = FGyROI1 + round(ptROI(4));
    
    % force selected ROI to be within the bounds of the image
    FGxROI1 = min(FGxROI1,size(im,2));
    FGxROI2 = min(FGxROI2,size(im,2));
    FGyROI1 = min(FGyROI1,size(im,1));
    FGyROI2 = min(FGyROI2,size(im,1));
    FGxROI1 = max(FGxROI1,1);
    FGxROI2 = max(FGxROI2,1);
    FGyROI1 = max(FGyROI1,1);
    FGyROI2 = max(FGyROI2,1);

    plot([BGxROI1,BGxROI2,BGxROI2,BGxROI1,BGxROI1],[BGyROI1,BGyROI1,BGyROI2,BGyROI2,BGyROI1],'r-')
    plot([FGxROI1,FGxROI2,FGxROI2,FGxROI1,FGxROI1],[FGyROI1,FGyROI1,FGyROI2,FGyROI2,FGyROI1],'r-')
    
    close(f2)

else   % auto definition of ROIs
    
    % Find transition regions for every axial profile to determine the
    % presence of holes in the top of the phantom.
    indexes = [];
    vertical_range = 1:FGyROI2;
    
    for j = FGyROI1:FGyROI2
        axial_profile = medfilt1(double(im(j,:)),5);
        axial_profile_diff = diff(axial_profile);
        [max_value, max_index] = max(axial_profile_diff);
        indexes(j) = max_index;
    end
    
    index_mode = mode(indexes(indexes>0));
    vertical_range = vertical_range(indexes==index_mode);
    FGyROI1 = max(vertical_range(1),1);
    FGyROI2 = max(vertical_range(end),1);
    plot([FGxROI1,FGxROI2,FGxROI2,FGxROI1,FGxROI1],[FGyROI1,FGyROI1,FGyROI2,FGyROI2,FGyROI1],'r-')
    title('Axial LSF: Automatic ROI definition: background, foreground and vertical analysis range.')
    
end

% Normalize averaged image such that values are [0,1]
im = (im - BGvalue) ./ (FGvalue - BGvalue);

%% Estimate the axial profile
min_y = FGyROI1;
max_y = FGyROI2;
im = sum( im(min_y : max_y, :), 1) / (max_y - min_y + 1);
z_smooth = 1;
[a,b] = size(im);
z_axis = pixelz : pixelz : b*pixelz;
z_axis_i = pixelz : rstep : b*pixelz;
j = min_y;
z_profile = im;
esf = interp1(z_axis,z_profile,z_axis_i,'linear');

%% smooth axial ESF and find maximum_derivative_indexmum of the LSF
fz = smooth(esf,z_smooth);
d = diff(fz);
[maxx,I] = max(d);
Iref = I;
lsf = diff(esf);
filt_lsf = smooth(lsf,smooth_f);
maxlsf = max(lsf);

rcent = find(lsf == maxlsf);  % maximum index of the LSF
lsfaxis = z_axis(1 : end-1);
pad = round( pad_mm / pixelz );
    
%% crop and pad the LSF
    st = rcent-pad;
    st(st<1) = 1;
    endd = rcent+pad;
    endd(endd>length(esf)) = [];
    st:endd;
    length(lsfaxis);
    lsfaxis = lsfaxis(st:endd);
    lsf = lsf(st:endd);
    nlsf = length(lsf);
    
    if(mod(nlsf,2))
    
    else
        lsf = lsf(1:end-1);
        nlsf = length(lsf);
    end
    
   
%% Antialiasing
    
T = fftshift(fft(lsf));  % FFT of the LSF

% Calculate Hann window and apply to data if specified.

    if window == 1
        w = 0.5 .* (1 - cos(2*pi*[1:length(lsfaxis)]./length(lsfaxis)));
    else
        w = ones(length(lsfaxis),1)';
    end

lsfw = lsf .* w;
lsf = lsfw;

T = fftshift(fft(lsfw/sum(lsfw(:))));
faxis = -1/(2*pixelz):1/((nlsf-1)*pixelz):1/(2*pixelz);
MTF = abs(T);
fmax = 1/(2*pixelz);

    if(use_antialias)
        MTF(abs(faxis)>fmax) = 0;
        T(abs(faxis)>fmax) = 0;
    end

    
fprintf('MTF calculated.\n\n');

%% Normalize the lsf area to unity
lsf = real(ifft(fftshift(T)));
lsf = double(lsf/sum(lsf(:)));
lsf = fliplr(lsf);

return;
