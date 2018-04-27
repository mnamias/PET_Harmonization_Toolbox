function [faxis, MTF, lsf, esf] = LSFcalc_R(im3D, dicom_headers, options)
%% Name of code: LSFcalc_R.m
% Version: 1.0 (April 20, 2018)
% Level of code: function. This function is from the PET harmonization toolbox GUI
% Author: Mauro Namías - mnamias@gmail.com
% This code is intended to accompany the paper:
% Namías et. al
% A novel approach to quantitative harmonization in PET
% PMB (2018)
%
% function [faxis, MTF, lsf, esf] = LSFcalc_R(im3D, dicom_headers, options);
%
% Estimates the radial line spread function from a PET cylindrical phantom
% filled with 18F-FDG.
%
% Inputs: 
%   im3D: the phantom slices (i.e.: 192x192x95)
%   dicom_headers: structure with the dicom headers for each slice
%   options: MATLAB structure with processing options
% 
% Options:
%   options.over_factor_R : Oversampling factor for radial binning (default = 4).
%   options.use_window_R : Decide whether to apply a Hann window to the LSF before calculating the MTF (0 or 1, typically applied)
%   options.thr_R : Optional threshold level for the ESF [0 1]. If set to zero, no thresholding is performed. 
%   options.use_antialias_R : (0 or 1). If true, apply an antialiasing filter in frequency domain to the MTF
%   options.pad_R : zero-padding length of the LSF in mm
%   options.correct_R_trend : (0 or 1). If true, fits a 2nd order polynomial to correct residual attenuation and scatter correction problems.
%   Determine which angles of data to use in radians (-pi to pi to use all)
%   options.thlim1 : initial angle
%   options.thlim2 : final angle
%
% Outputs:
%   faxis: frequency axis in mm-1
%   MTF: estimated modulation transfer function (unitless)
%   lsf: estimated line spread function 
%   esf: estimated edge spread function 
%
% IMPORTANT!
% This code was adapted for PET imaging by Mauro Namías (2018) from the original code from
% Friedman et. al.:
% 
% S. N. Friedman, G. S. K. Fung, J. H. Siewerdsen, and B. M. W. Tsui.  "A
% simple approach to measure computed tomography (CT) modulation transfer 
% function (MTF) and noise-power spectrum (NPS) using the American College 
% of Radiology (ACR) accreditation phantom," Med. Phys. 40, 051907-1 -  
% 051907-9 (2013).
% http://dx.doi.org/10.1118/1.4800795  

% Copyright (c) 2012, Saul N. Friedman
% Copyright (c) 2010, Marian Uhercik
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
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

% User selected options
    over_factor = options.over_factor_R;  % Oversampling factor for radial binning
    window = options.use_window_R ; % Decide whether to apply a Hann window to the LSF before calculating the MTF (typically applied)
    thr = options.thr_R; % optional threshold level for the ESF [0 1]. If set to zero, no thresholding is performed. 
    pad_mm = options.pad_R; % zero-padding length of the LSF in mm
    correct_R_trend = options.correct_R_trend; % If true, fits a 2nd order polynomial to correct residual attenuation and scatter correction problems.

% Determine which angles of data to use in radians (-pi to pi to use all)
    thlim1 = options.thlim1;
    thlim2 = options.thlim2;

% pixel spacing and image sizes
[a,b,c] = size(im3D);
pixelx = dicom_headers(1).PixelSpacing(2);
pixely = dicom_headers(1).PixelSpacing(1);
smoothf = 141 * 0.1/(pixelx/5);
rstep = pixelx / over_factor; %  Choose bin size in mm (i.e., effective pixel size) for oversampled edge spread function (ESF)

% Background and foreground values:  

   % box border for background ROI in pixels
   BGxROI1 = 1;
   BGxROI2 = round(500*0.4883/pixelx);
   BGyROI1 = 1;
   BGyROI2 = round(40*0.4883/pixely);
   
   % box border for foreground ROI in pixels
   FGxROI1 = round(205*0.4883/pixelx);
   FGxROI2 = round(355*0.4883/pixelx);
   FGyROI1 = round(208*0.4883/pixely);
   FGyROI2 = round(358*0.4883/pixely);

% 1 = prompted to draw rectangle(s) to select ROI(s), 0 = use values below 



%% Determine analysis range

mip = mipss(im3D);  % creates a MIP of the cylindrical phantom
z_profile = squeeze(sum(sum(im3D,1),2));

bkg = z_profile - z_profile;
bkg( z_profile < max(z_profile(:)) * 0.1 )= 1;

bkg_indexes = find(bkg);
bkg_indexes2 = bkg_indexes(2:end)-bkg_indexes(1:end-1);
slice = find(bkg_indexes2-1)

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



    mip = mipss(im3D);
    f1 = figure;
    imagesc(-mip);
    colormap gray
    hold on;
    axis image;
    impixelinfo;
    %title('Radial LSF: axial range definition')
    %title( {'Click and drag to create a rectangular ROI representing the axial range.' 'Double click the ROI when finished.'})
    
    manual_flag = strcmp(options.mode, 'manual');
    
    if(manual_flag) % manual definition of ROIs
        title( {'Click and drag to create a rectangular ROI representing the axial range.' 'Double click the ROI when finished.'})

  

        fprintf('Click and drag to create a rectangular ROI representing the axial analysis range.  Double click the ROI when finished.\n\n');
        h = imrect;
        ptROI = wait(h);
        
        BGxROI1 = round(ptROI(1));
        BGxROI2 = BGxROI1 + round(ptROI(3));
        BGyROI1 = round(ptROI(2));
        BGyROI2 = BGyROI1 + round(ptROI(4));

        % force selected ROI to be within the bounds of the image
        BGxROI1 = min(BGxROI1,size(mip,2));
        BGxROI2 = min(BGxROI2,size(mip,2));
        BGyROI1 = min(BGyROI1,size(mip,1));
        BGyROI2 = min(BGyROI2,size(mip,1));
        BGxROI1 = max(BGxROI1,1)
        BGxROI2 = max(BGxROI2,1)
        BGyROI1 = max(BGyROI1,1);
        BGyROI2 = max(BGyROI2,1);

        im3D = im3D(:,:,  BGxROI1: BGxROI2);  % crop the image to the selected analysis range
        z_range = BGxROI1: BGxROI2;

        if(manual_flag)
            close(f1)
        end
        
    
    else   % Automatic ROIs
        z_profile = squeeze(sum(sum(im3D,1),2));
        dz_profile = diff(z_profile(1:end/2)); 
        maximum_derivative_index = find( abs(dz_profile) == max(abs(dz_profile(:))) );
        z_range = maximum_derivative_index + round(30/pixelx) : maximum_derivative_index + round(80/pixelx);
        z_range(z_range>size(im3D,3)) = [];
        z_range(z_range<1) = [];
        
        plot( [min(z_range) min(z_range)], [size(im3D,1)*0.2 size(im3D,1)*0.8], 'g --');
        plot( [max(z_range) max(z_range)], [size(im3D,1)*0.2 size(im3D,1)*0.8], 'g --');
        title('Radial LSF estimation: Automatic ROI definition (axial analysis range)');
        
        
    end
    
    
% Average the axial analysis range to make a 2D slice  
    imavg = sum(im3D,3);
    im = imavg;
 
% Show averaged image 
    imaxisx = [1:size(im,2)] .* pixelx;
    imaxisy = [1:size(im,1)] .* pixely;
   
    f2 = figure
    imagesc(-imavg, [-max(imavg(:))*1.5 0])
    xlabel('pixel')
    ylabel('pixel')
    colormap gray
    hold on
    axis image;
    impixelinfo;
    title('Averaged 2D Slice for LSFr processing')

% Determine background and foreground values to properly normalized the
% data to range [0,1]

    if manual_flag
        
        title( {'Click and drag to create a rectangular ROI representing the background.' 'Double click the ROI when finished.'})

        fprintf('Click and drag to create a rectangular ROI representing the background.  Double click the ROI when finished.\n\n');
        h = imrect;
        ptROI = wait(h);
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
    
        BGvalue = mean(mean(imavg(BGyROI1:BGyROI2,BGxROI1:BGxROI2)));

        title( {'Click and drag to create a rectangular ROI representing the foreground.' 'Double click the ROI when finished.'})

        fprintf('Click and drag to create a rectangular ROI representing the foreground.  Double click the ROI when finished.\n\n');
        h = imrect;
        ptROI = wait(h);

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
    
        FGvalue = mean(mean(imavg(FGyROI1:FGyROI2,FGxROI1:FGxROI2)));

        plot([BGxROI1,BGxROI2,BGxROI2,BGxROI1,BGxROI1],[BGyROI1,BGyROI1,BGyROI2,BGyROI2,BGyROI1],'r-')
        plot([FGxROI1,FGxROI2,FGxROI2,FGxROI1,FGxROI1],[FGyROI1,FGyROI1,FGyROI2,FGyROI2,FGyROI1],'r-')
        hold off
        
        close(f2)
        
        imavg = (imavg - BGvalue) ./ (FGvalue - BGvalue);     % Normalize averaged image such that values are [0,1]

    else  % Automatic ROI 
        h = fspecial('gaussian',5,3);
        fimavg = imfilter(imavg,h);

        s = regionprops(fimavg>max(fimavg(:))*0.1, fimavg, {'Centroid','WeightedCentroid'});
        % Calculate the center of the cylinder

        centroid.y = s.Centroid(2);
        centroid.x = s.Centroid(1);

        xaxis = ([1:size(imavg,2)] - centroid.x) .* pixelx;
        yaxis = ([1:size(imavg,1)] - centroid.y) .* pixely;

        % Convert all image coords to polar coords with origin at center of the cylinder
        r = zeros(size(imavg));
        th = zeros(size(imavg));

        for i = 1 :size(imavg,2)
            for j = 1 : size(imavg,1)
                [th(j,i) r(j,i)] = cart2pol(xaxis(i),yaxis(j));
            end
        end

        % Create ESF
        rmax = max(max(r));
        esf = zeros(ceil(rmax/rstep)+1,1);
        nsamp = zeros(length(esf),1);
        rbin = 0:pixelx:rmax+rstep;

        % data with radius falling within a given bin are averaged together for a
        % low noise approximation of the ESF at the given radius

        thlim1 = min(thlim1,thlim2);
        thlim2 = max(thlim1,thlim2);
    
        h2 = msgbox('Estimating radial ESF, please wait...');

        for i = 1:length(rbin)
 %           i/length(rbin)*100
            R1 = find(r >= rbin(i));
            R2 = find(r < rbin(i) + rstep);
            R3 = intersect(R1,R2);
            R4 = find(th >= thlim1);
            R5 = find(th < thlim2);
            R6 = intersect(R5,R4);
            R7 = intersect(R3,R6);
            R8 = R7;
            R=R8(R8~=0);

            [X Y] = ind2sub(size(imavg),R);

            esf(i) = sum(imavg(R));
            nsamp(i) = length(R);
        end
        close(h2);
        
        i1 = find(nsamp,1,'first');
        i2 = find(nsamp,1,'last');
        nsamp = nsamp(i1:i2);
        esf = esf(i1:i2);
        rbin = rbin(i1:i2);

        I = find(nsamp > 0);
        esf(I) = esf(I)./nsamp(I);
        %esf(esf(rbin<80)<0.8) = mode(esf(esf>max(esf(:))*0.5));  %remove outliers
        esf_diff = diff(esf);
        [max_value max_index] = max(abs(esf_diff));
        
        FGvalue = mean( esf( round(max_index*0.2) : round(max_index*0.7) ));
        BGvalue = mean( esf( round(max_index*1.1) : round(max_index*1.5) ));
        imavg = (imavg - BGvalue) ./ (FGvalue - BGvalue);     % Normalize averaged image such that values are [0,1]

        end
    
   
    f3 = figure
    imagesc([imaxisy(1) imaxisy(end)],[imaxisx(1) imaxisx(end)], -imavg, [-1.5 0])
    xlabel(' mm')
    ylabel(' mm')
    title('Radial LSF processing: Averaged and normalized [0, 1] image')
    colormap gray
    hold on
    axis image;
    impixelinfo;

% Find phantom centroid

    h = fspecial('gaussian',5,3);
    fimavg = imfilter(imavg,h);

    s = regionprops(fimavg>max(fimavg(:))*0.1, fimavg, {'Centroid','WeightedCentroid'});
    imlog = imavg;

% Calculate the center of the cylinder
    centroid.y = s.Centroid(2);
    centroid.x = s.Centroid(1);

    xaxis = ([1:size(imavg,2)] - centroid.x) .* pixelx;
    yaxis = ([1:size(imavg,1)] - centroid.y) .* pixely;

% Convert all image coords to polar coords with origin at center of the
% cylinder

    r = zeros(size(imavg));
    th = zeros(size(imavg));

    for i = 1 :size(imavg,2)
        for j = 1 : size(imavg,1)
            [th(j,i) r(j,i)] = cart2pol(xaxis(i),yaxis(j));
        end
    end

% Create oversampled ESF excluding data specified in mask

    rmax = max(max(r));

    esf = zeros(ceil(rmax/rstep)+1,1);
    nsamp = zeros(length(esf),1);

    rbin = 0:rstep:rmax+rstep;


% data with radius falling within a given bin are averaged together for a
% low noise approximation of the ESF at the given radius

    fprintf('Calculating the MTF.\n');
    fprintf('The may take several minutes. \n');
    fprintf('Please be patient.\n\n');

    thlim1 = min(thlim1,thlim2);
    thlim2 = max(thlim1,thlim2);

    h2 = msgbox('Estimating radial ESF, please wait...');
    for i = 1:length(rbin)
        %i/length(rbin)*100
        R1 = find(r >= rbin(i));
        R2 = find(r < rbin(i) + rstep);
        R3 = intersect(R1,R2);
        R4 = find(th >= thlim1);
        R5 = find(th < thlim2);
        R6 = intersect(R5,R4);
        R7 = intersect(R3,R6);
        R8 = R7;
        R=R8(R8~=0);

        [X Y] = ind2sub(size(imavg),R);    

        esf(i) = sum(imavg(R));
        nsamp(i) = length(R);
    end
    
    close(h2)

    i1 = find(nsamp,1,'first');
    i2 = find(nsamp,1,'last');
    nsamp = nsamp(i1:i2);
    esf = esf(i1:i2);
    rbin = rbin(i1:i2);

    I = find(nsamp > 0);
    esf(I) = esf(I)./nsamp(I);
    esf(esf(rbin<80)<0.8) = 1;  %remove outliers

    if(manual_flag)
        close(f3)
    end
    
    % Plot oversampled ESF

    f4 = figure
    plot(rbin,smooth(esf,10))
    xlabel('Radius (mm)')
    ylabel('Normalized ESF')
    grid on
    axis([rbin(1) rbin(end) -1.2*min(esf) 1.1*max(esf)])

    sesf = esf;

    if(correct_R_trend)
        % if desired, fit a 2nd order polynomial to the ESF to correct residual
        % attenuation and scatter correction errors

        esf_lf = sesf(round(15/rstep:75/rstep));
        rbin_lf = rbin(round(15/rstep:75/rstep));
        p = polyfit(rbin_lf,medfilt1(esf_lf',5),2);
        esf_fit = polyval(p,rbin);
        plot(rbin,sesf,'k -', rbin,esf_fit,'g -' )
        xlabel('radius (mm)')
        ylabel('normalized + oversampled ESF')
        grid on
        axis([1 150 -0.25 1.25])
        title('Radial ESF: Low Frequency Trend Correction');    
        esf = esf./esf_fit';

        close(f4)
        
        f5 = figure;
        plot(rbin,esf,'k -')
        xlabel('Radius (mm)')
        ylabel('Normalized ESF')
        grid on
        axis([1 150 -0.25 1.25])
        title('Radial ESF: Low Frequency Trend Correction');    
    end
    
    if(manual_flag)
        close(f5)
    end
    

% Calculate the LSF from the ESF
    lsf = diff(esf);
    
    % artifacts caused by empty bins likely have a value of -/+ 1, so try to 
    % remove them.
    I = find(abs(lsf)> 0.9);
    lsf(I) = 0;

    lsf(rbin<10) = 0;

    slsf = smooth(lsf,5);

    if max(slsf) < max(abs(slsf))
        lsf = -lsf;
    end


    %% Centering and padding of the LSF
        if(thr >=0)
            st = 70/rstep;
            endd = 130/rstep;
            maxlsf = max(abs(lsf(st:endd)));
        else
            maxlsf = max(abs(lsf));
        end


    rcent = find(abs(lsf) == maxlsf);
    pad = round(pad_mm/rstep);
    st = rcent-pad;
    endd = rcent+pad;
        if(st<1)
            st = 1;
        end
        
        if(endd>length(lsf))
            endd = length(lsf);
        end

    lsf = lsf(st:endd);
    lsfaxis = -length(lsf)/2*rstep:rstep:length(lsf)/2*rstep;  

    if(length(lsfaxis)>length(lsf))
        lsfaxis = lsfaxis(1:length(lsf));
    end

    slsf = smooth(lsf,5);
    lsf_norm = abs(slsf)/max(abs(slsf(:)));
    nlsf = length(lsf);


    %calculate Hann window and apply to data if specified.
    if window == 1
        w = 0.5 .* (1 - cos(2*pi*[1:length(lsf)]./length(lsf)));
        w = w';
    else
        w = ones(length(lsf),1);
    end

    if( thr > 0 && thr < 2)
        lsfw = lsf .* w';
    else        
        if(thr ==0)
            lsfw = lsf .* w;
        end
        
        if(thr<0)
            w = ones(length(lsf),1);
            lsfw = lsf .* w';
        end
    end

    lsf = lsfw;

% Plot windowed LSF 

    
%     f6 = figure
%     whos lsfaxis
%     whos lsf
%     plot(lsfaxis, lsf, 'k -')
%     xlabel('Position (mm)')
%     ylabel('Normalized and Windowed LSF')
%     axis([lsfaxis(1) lsfaxis(end) 1.1*min(lsfw) 1.1*max(lsfw)])
%     grid on


    % Calculate the MTF from the LSF
    lsfw = lsfw/sum(lsfw(:));
    T = fftshift(fft(lsfw));
    faxis = -1/(2*rstep):1/((nlsf-1)*rstep):1/(2*rstep);
    MTF = abs(T);

    if(thr == 2)
        MTF(MTF>1) = (MTF(MTF>1)-1)*2+1;
    end

    fmax = 1/(2*dicom_headers(1).PixelSpacing(1));

    % Optional antialiasing
    if(options.use_antialias_R && thr>0)
        filter = faxis-faxis+1;
        filter(abs(faxis)>fmax*1)= 0;
        length(filter)
        filter = smooth(filter,round(length(filter)*fmax/4))';
        MTF = MTF.*filter;
        T = T.*filter;
    end

    if(options.use_antialias_R && thr==0)
        filter = faxis-faxis+1;
        filter(abs(faxis)>fmax*1)= 0;        
        MTF = MTF.*filter';
        T = T.*filter';    
    end

% Plot the MTF

%     f7 = figure;
%     semilogx(faxis,abs(T),'k.-');   
%     xlabel('Spatial Frequency (mm^{-1})');
%     ylabel('Radial MTF')
%     grid on
%     fprintf('MTF calculated.\n\n');

    % Downsample lsf if needed
    ratio = rstep/pixelx;
    lsf = real(ifft(fftshift(T)));

    maxi = find(abs(lsf)==max(abs(lsf(:))));
    pad = round(20*over_factor/pixelx);
    lsf = lsf(maxi-pad:maxi+pad);

    lsf = double(lsf/sum(lsf(:)));
  
 
return;