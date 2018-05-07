function [Wavg, Savg, rbin, npsr, npsr_norm, avg_std, mean_bkg, psf, xaxis, yaxis, zaxis] = NPScalc( image , dicom_headers, options )
% Name of code: NPScalc.m
% Version: 1.0 (April 20, 2018)
% Name of author: Mauro Namías
%
% function [Wavg, Savg, rbin, npsr, npsr_norm, avg_std, mean_bkg, psf, xaxis, yaxis, zaxis] = NPScalc( image , dicom_headers )
%
% Overview and purpose of the code: This code is intended to accompany the paper:
% M. Namías, T. Bradshaw, V.O. Menezes, M.A.D. Machado and R. Jeraj. "A novel
% approach for quantitative harmonization in PET". Physics in Medicine and
% Biology, 2018. 
% In particular, NPScalc is a function to estimate the 3D Noise Power Spectrum
% (NPS) from a cylindrical PET phantom.
%
% Inputs and range of inputs:
%   image: a 3D matrix with the full PET phantom volume (i.e.: 192x192x100) 
%   dicom_headers: a MATLAB structure of dicom headers for each slice
% 
% Outputs:
%   Wavg: average 3D Noise Power Spectrum matrix
%   Savg: averag Amplitude Spectrum
%   rbin: radial bins in mm
%   npsr: radial average of the NPS
%   npsr_norm: normalized radial average of the NPS
%   avg_std: average standard deviation of the phantom volume
%   mean_bkg: mean value of the active phantom volume
%   psf: Noise PSF
%   xaxis: x frequency axis 
%   yaxis: y frequency axis 
%   zaxis: z frequency axis 
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


image = double(image);
manual_flag = strcmp(options.mode, 'manual');  % manual or automatic processing

%% Determine the phantom analysis range


%%
    mip = mipss(image); % show a mip of the phantom image    
    f1 = figure;
    imagesc(-mip, [-max(mip(:))*1.5 0])
    colormap gray
    hold on
    axis image;
    impixelinfo;
    
    
    if(manual_flag)
    
    title( {'Click and drag to create a rectangular ROI representing the axial range.' 'Double click the ROI when finished.'})    
    fprintf('Click and drag to create a rectangular ROI representing the analysis range.  Double click the ROI when finished.\n\n');
     
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
    BGxROI1 = max(BGxROI1,1);
    BGxROI2 = max(BGxROI2,1);
    BGyROI1 = max(BGyROI1,1);
    BGyROI2 = max(BGyROI2,1);
    
    else % auto processing
        
        z_profile = squeeze(sum(sum(image,1),2));
     
        z_range = find(z_profile > max(z_profile(:)) * 0.8);
        z_range( z_range > size(image, 3)) = [];
        z_range(z_range < 1) = [];
        
        plot( [min(z_range) min(z_range)], [size(image,1)*0.2 size(image,1)*0.8], 'g --');
        plot( [max(z_range) max(z_range)], [size(image,1)*0.2 size(image,1)*0.8], 'g --');
        title('Noise Power Spectrum estimation: Automatic ROI definition (axial analysis range)');
        BGxROI1 = z_range(1);
        BGxROI2 = z_range(end);
        
    end
    
    image = image(:,:,  BGxROI1: BGxROI2);

    imavg = sum(image,3);  % get an average slice by summing in the Z-direction
    [a,b,c] = size(image);



% Filter the average image to reduce noise and find the centroid
    h = fspecial('gaussian',[31 31], 27);
    fimavg = imfilter(imavg,h);
    s = regionprops(fimavg>max(fimavg(:))*0.2, fimavg, {'Centroid','WeightedCentroid'});
    cx = s.WeightedCentroid(1);  %x coordinate of the centroid
    cy = s.WeightedCentroid(2);  %y coordinate of the centroid

% Center the image 
    dx = round(b/2-cx);
    dy = round(a/2-cy);
    image = circshift(image,dy,1);  
    image = circshift(image,dx,2);

% Processing options 
    calczero = 0;  % when set to 0, the NPS(0,0,0) value is estimated based on the surrounding 26 voxel values to try remove any zero-frequency artifacts
    use_window = 0; % use a Kaiser window to smooth the spectrum;
    rstep = 0.025;  % This will be the effective bin size after forming the radially averaged NPS.
    calcrstep = 1;

% Image size and pixel dimensions
    imheader = dicom_headers(1);        
    pixelx = double(imheader.PixelSpacing(1));
    [Ny,Nx,Nz] = size(image);

% Find overlapping positions of the VOIs
    deltax = round(60/pixelx); % spacing of the VOIs
    if(mod(deltax,2))
    else
        deltax = deltax+1;
    end

    deltay = deltax;

    calcdeltaz = 1; % 1 = calculate size of ROI in z (slices) to be same physical size as x dimension, 0 = use value below
    deltaz = 51; %specify number of slices to use in ROI

% Specify overlap fraction for ROIs in z direction
%  N.B. value is typically [0,0.5] and represents overlap with EACH
%  neighboring ROI (i.e., at 0.5 all pixels overlap with at least one
%  neighbor)
    overlapz = 0.5;  

% Parameters pertaining to placement of ROI in the axial plane.
%  N.B. overlap of ROIs will depend on a combination of the radius chosen
%  in rcenter and the distance between the centers of ROIs chosen.
    rcenter = 45; % radius in mm along which to center ROIS, a value of 50 mm corresponds to half the radius of the phantom
    deltar = 30; % distance along the radial arc in mm between centers of ROIs in the axial plane
    NRoisXY = 12; % Number of equally spaced ROIs in the tangential direction 
    deltaangle = 2*pi/180;  % the range of angles in radians over which to look for a suitable pixel on which to center an ROI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of analysis code

        [Ny,Nx,Nz] = size(image);
       
        Nz = int16(floor(size(image,3)/2));
        im = image(:,:,round(end/2));  % middle slice
  
        Rows = imheader.Rows;        
        diameter = Rows*pixelx;
    
        % Pixel sizes
        pixelx = double(imheader.PixelSpacing(1));
        pixely = double(imheader.PixelSpacing(2));
        pixelz = double(imheader.SliceThickness);      

                
        Nx = int16(round(size(im,2)));
        Ny = int16(round(size(im,1)));
        imaxisx = [1:Nx] .* pixelx;
        imaxisy = [1:Ny] .* pixely;
        
        % Initialize sub-volumes for NPS estimation
        imavg = zeros(Ny,Nx);
        im3D1 = zeros(Ny,Nx,Nz);
        im3D2 = zeros(Ny,Nx,Nz);
  
    
    % Get subvolumes V1 and V2
    
    for z = 1:Nz
        im = double(image(:,:,z));                
        im3D1(:,:,z) = im;
        imavg = imavg + im;
    end
    
    for z = Nz+1:Nz*2
        im = double(image(:,:,z));
        im3D2(:,:,z-Nz) = im;
        imavg = imavg + im;
    end
    
    imavg = imavg /  double(Nz*2);

%% find mean phantom value
   median_im = medfilt2(imavg,[9 9]);
   mask = median_im>max(median_im(:))*0.2;
   mean_bkg = mean(imavg(mask(:)==1));

    if deltaz  > Nz
        deltaz = Nz;
        val = sprintf('Selected ROI size in z is larger than the number of slices.  Using value of %d instead.\n\n',Nz);
        fprintf(val)
    end


% Detrend the data and plot if desired

    fprintf('Detrending data.\n\n');
    im3Dzm = im3D1 - im3D2;
    

% Calculate the center of the cylinder
    fprintf('Calculating the center of the phantom and determining pixel locations.\n\n');
    
    imlog = imavg;
    if exist('bwconncomp')==2
        CC = bwconncomp(imlog,26);
        for i = 1:CC.NumObjects;
            Lobj(i) =  length(CC.PixelIdxList{1,i}) ;
        end
    else
        [CC CCn] = bwlabel(imlog,8);
        for i = 1:CCn
            Lobj(i) = length(find(CC == i));
        end   
    end

    L = find(Lobj == max(Lobj));
    S = regionprops(CC,'Centroid');

    h = fspecial('gaussian',[17 17], 17);
    fimavg = imfilter(imavg,h);  % Filtered average image
    s = regionprops(fimavg>max(fimavg(:))*0.1, fimavg, {'Centroid','WeightedCentroid'});
    centroid.y = round(s.Centroid(2));
    centroid.x = round(s.Centroid(1));

    % Plot the average image with the centroid marked
    f2 = figure;
    imagesc([imaxisx(1) imaxisx(end)], [imaxisy(1) imaxisy(end) ], -imavg, [-max(imavg(:))*1.5 0])
    xlabel('mm')
    ylabel('mm')
    title('NPS estimation: Cylinder with center and VOIs marked')
    colormap gray
    hold on
    axis image;
    impixelinfo;
 
    if exist('bwconncomp')==2
        [I J] = ind2sub(size(imavg),CC.PixelIdxList{1,L});
    else
        [I,J] = find(CC == L);
    end
    
    plot(centroid.x*pixelx,centroid.y*pixely,'*y')

    xaxis = ([1:size(imavg,2)] - centroid.x) .* pixelx;
    yaxis = ([1:size(imavg,1)] - centroid.y) .* pixely;

% Calculate polar coordinates for to place VOIs at a specified radial
% direction

    for i = 1 :size(imavg,2)
        for j = 1 : size(imavg,1)
            [th(j,i) r(j,i)] = cart2pol(xaxis(i),yaxis(j));
        end
    end

    rcalc = rcenter;
    anglestep = 2*pi/NRoisXY;

% Calculate ROI size in z direction such that it is the same physical
% dimension as the x ROI, if desired.

    if calcdeltaz == 1
        deltaz = round(deltax * pixelx ./ pixelz);
    end

    if(mod(deltaz,2))
 
    else
        deltaz=deltaz-1;
    end

    if(deltaz > floor(size(image,3)/2))
        deltaz = floor(size(image,3)/2);
    end


% Determine the number of ROIs in each direction.

    Ndeltaxy = floor( 2 * pi / anglestep );    
    final = 0;
    Ndeltaz = 0;
    while final < Nz
        Ndeltaz = Ndeltaz + 1   ;

        if(Ndeltaz == 1)
            final = final + deltaz;
        else
            final = final + deltaz * (1 - overlapz);
        end       

    end
    Ndeltaz = Ndeltaz - 1;
    

% Determine the offset needed to center the ROIs in z direction.
    offsetz = round((Nz - deltaz - round((deltaz - overlapz * deltaz)*(Ndeltaz - 1)))/2);
    if(offsetz <1) 
        offsetz = 1;
    end

    ROIavg = zeros(deltay,deltax,deltaz);
    nROI = 0;

    for i = 1:Ndeltaxy;
    
        R1 = find( abs(r - rcalc) < pixelx);
        R2 = find( th < -pi + i*anglestep + deltaangle);
        R3 = find( th > -pi + i*anglestep - deltaangle);
        R4 = intersect(R1,R2);
        R5 = intersect(R3,R4);
        R = min(R5);
        [I J] = ind2sub(size(im),R);
        y1 = J-floor(deltay/2)+1;
        y2 = J+round(deltay/2);
        x1 = I-floor(deltax/2)+1;
        x2 = I+round(deltax/2);
     
        for z = 1:Ndeltaz
            nROI = nROI + 1;
            z1 =  offsetz + (deltaz * (z-1)) * (1-overlapz);
            z2 = z1+deltaz-1;
            if(z2>Nz)
                z2 = Nz;
                z1 = z2-deltaz+1;
                disp('Warning! bad z ROI');
                ROIzm(:,:,:,nROI) =  im3Dzm(y1:y2,x1:x2,z1:z2);
            else
                ROIzm(:,:,:,nROI) =  im3Dzm(y1:y2,x1:x2,z1:z2) ;
            end
        end
    
    
        plot(I* pixelx,J* pixelx,'r.')
        plot([x1,x2,x2,x1,x1]* pixelx,[y1,y1,y2,y2,y1]* pixelx,'r-')
    
end

% Calculate power spectrum realizations from each ROI.  The ensemble
% average of these is the calculated estimate of the NPS.

    fprintf('Calculating the NPS.\n');
    fprintf('Please be patient.\n\n');

% Initialize variables
    Wavg = double(zeros(size(ROIzm,1),size(ROIzm,2),size(ROIzm,3)));
    Savg = double(zeros(size(ROIzm,1),size(ROIzm,2),size(ROIzm,3)));
    avg_std = 0;

    if(use_window ==1)
        window = kaiser3D(size(Wavg),8);
    else
        window = ones(size(Wavg));
    end

    for i = 1:nROI
        R = ROIzm(:,:,:,i);
        m = mean(R(:));
        W = abs( fftn(((ROIzm(:,:,:,i)-m).*window)) ) .^2;
        W = fftshift(W);
        dummy = ROIzm(:,:,:,i);
        avg_std = avg_std + std(dummy(:))/nROI;
        Wavg = W/nROI + Wavg;
        Savg = Savg + W.^(1/2)/nROI;
    end

        avg_std = avg_std/sqrt(2);

        Wavg = Wavg /double(deltax .* deltay .* deltaz .* 2 ) .* double(pixelx .* pixely .* pixelz);
        factor = double(deltax .* deltay .* deltaz .* 2 ) .* double(pixelx .* pixely .* pixelz);
        Savg = Savg /sqrt(2);

        % Calculate the axes values in the frequency domain.
        xaxis = -(1/(2*pixelx)):(1/(2*pixelx))/(floor(size(ROIzm,2)/2)):(1/(2*pixelx));
        yaxis = -(1/(2*pixely)):(1/(2*pixely))/(floor(size(ROIzm,1)/2)):(1/(2*pixely));
        zaxis = -(1/(2*pixelz)):(1/(2*pixelz))/(floor(size(ROIzm,3)/2)):(1/(2*pixelz));
        
% The calculated zero-frequency value often contains an artifact.  The value
% can be approximated based on an average of the surrounding values.
    if calczero == 0;
        J = find(abs(yaxis)==min(abs(yaxis)));
        I = find(abs(xaxis)==min(abs(xaxis)));
    
        for k = 1:size(Wavg,3)
             Wavg(J,I,k) = ( Wavg(J,I+1,k) + Wavg(J,I-1,k) + Wavg(J+1,I,k) + Wavg(J+1,I+1,k) + Wavg(J+1,I-1,k) + Wavg(J-1,I,k) + Wavg(J-1,I+1,k) + Wavg(J-1,I-1,k) + Wavg(J,I+2,k) + Wavg(J,I-2,k) + Wavg(J+2,I,k) + Wavg(J+2,I+2,k) + Wavg(J+2,I-2,k) + Wavg(J-2,I,k) + Wavg(J-2,I+2,k) + Wavg(J-2,I-2,k)) ./ 16;
        end
        
        for k = 1:size(Savg,3)
             Savg(J,I,k) = ( Savg(J,I+1,k) + Savg(J,I-1,k) + Savg(J+1,I,k) + Savg(J+1,I+1,k) + Savg(J+1,I-1,k) + Savg(J-1,I,k) + Savg(J-1,I+1,k) + Savg(J-1,I-1,k) + Savg(J,I+2,k) + Savg(J,I-2,k) + Savg(J+2,I,k) + Savg(J+2,I+2,k) + Savg(J+2,I-2,k) + Savg(J-2,I,k) + Savg(J-2,I+2,k) + Savg(J-2,I-2,k)) ./ 16;
        end

    end

% Plot a 1D slightly off-centered NPS profile in each direction.
    if(mod(deltaz,2))
        eps = 0;
    else
    eps = 1;
    end
    
   
        f3 = figure;
        subplot(3,1,1);
        plot(xaxis,Wavg(:,ceil(end/2),ceil(end/2)),'k.-');
        xlabel('x-dir (mm^{-1})');
        ylabel('NPS_x ((Bq/ml)^2 mm^3)');
        subplot(3,1,2);
        plot(yaxis,Wavg(ceil(end/2),:,ceil(end/2)),'k.-');
        xlabel('y-dir (mm^{-1})');
        ylabel('NPS_y ((Bq/ml)^2 mm^3)');
        subplot(3,1,3);
        zvalue = squeeze(Wavg(ceil(end/2),ceil(end/2),:));
        plot(zaxis,zvalue,'k.-');
        xlabel('z-dir (mm^{-1})')
        ylabel('NPS_z ((Bq/ml)^2 mm^3)');
   


% Calculate radially averaged NPS profile for xy plane

    % Switch to polar coords
    imslice = Wavg(:,:,floor(deltaz/2)+eps);
    disp(floor(deltaz/2)+eps);
    r = zeros(size(imslice));
    th = zeros(size(imslice));

    for i = 1 :size(imslice,2)
        for j = 1 : size(imslice,1)
            [th(j,i) r(j,i)] = cart2pol(xaxis(i),yaxis(j));
        end
    end

% Create oversampled radial NPS profile
    rmax = 1/(2*pixelx);

    if(calcrstep)
        rstep =   1/(2*pixelx)/double(size(ROIzm,1)/2);
    end

    npsr = zeros(ceil(rmax/rstep)+1,1);
    nsamp = zeros(length(npsr),1);
    rbin = 0:rstep:rmax+rstep;


% data with radius falling within a given bin are averaged together for a
% low noise approximation of the NPS at the given radius

    for i = 1:length(rbin)
        R1 = find(r >= rbin(i));
        R2 = find(r < rbin(i) + rstep);
        R = intersect(R1,R2);
        
        [X Y] = ind2sub(size(imslice),R);
        
        npsr(i) = sum(imslice(R));
        nsamp(i) = length(R);
    end

    i1 = find(nsamp,1,'first');
    i2 = find(nsamp,1,'last');
    nsamp = nsamp(i1:i2);
    npsr = npsr(i1:i2);
    rbin = rbin(i1:i2);

    I = find(nsamp > 0);
    npsr(I) = npsr(I)./nsamp(I);
    npsr_norm = npsr/avg_std.^2;

% Plot radially averaged NPS
    
    f4 = figure;
    semilogy(rbin,npsr,'.-')
    xlabel('spatial frequency (mm^{-1})')
    ylabel('NPS_r (Bq/ml^2 mm^3)')
    grid on
    axis([rbin(1) rbin(end) -1.2*min(npsr) 1.1*max(npsr)])
    title('Average Radial Noise Power Spectrum')
    
    fprintf('NPS calculated.\n\n');

% Estimate the Noise PSF
    psf = real(fftshift(ifftn(circshift(ifftshift(double(Savg)),[0 0 0]))));
    psf = double(psf);
    psf = psf/sum(psf(:));    
    spectrum = double(sqrt(Wavg));
    
    psf = real(ifftshift(ifftn(fftshift(spectrum))));
    psf = psf/sum(psf(:));

    if(manual_flag)
    close(f1)
    close(f2)
    close(f3)
    close(f4)
    end
    
    
return;