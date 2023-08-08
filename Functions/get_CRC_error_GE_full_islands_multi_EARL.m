function [error, fails, results] = get_CRC_error_GE_full_islands_multi_EARL(imagenes, data,  maskx1, maskx2, masky1, masky2, I, bkg_voi, pixel_sizes, FWHM, zfilter, brand, EARL_mode, half_life);

[a,b,c,d] = size(imagenes)

if(FWHM>15)
    FWHM = 15-randn;
end


if(FWHM>0)

if(strcmp(brand,'GE'))
    
for i = 1:d
%im_filt = filtra_3D(imagen, pixel_sizes, [FWHMs(1) FWHMs(1) FWHMs(2)]);
im_filt = filtra_3D(imagenes(:,:,:,i), pixel_sizes, [FWHM FWHM 0]);
imagenes(:,:,:,i) = im_filt;
end

    else
  
    
for i = 1:d
%im_filt = filtra_3D(imagen, pixel_sizes, [FWHMs(1) FWHMs(1) FWHMs(2)]);
im_filt = filtra_3D(imagenes(:,:,:,i), pixel_sizes, [FWHM FWHM FWHM]);
imagenes(:,:,:,i) = im_filt;
end

    end

end
 





if(zfilter >0)
z_kernel = zeros(1,1,3);
z_kernel(:,:,1) = 1;
z_kernel(:,:,2) = zfilter;
z_kernel(:,:,3) = 1;

z_kernel = z_kernel/sum(z_kernel(:));

%im_filt = imfilter(im_filt,z_kernel);

for i = 1:d
%im_filt = imfilter(im_filt,z_kernel);
imagenes(:,:,:,i) = imfilter(imagenes(:,:,:,i) ,z_kernel);
end


end

results = [];

error = 0;
fails = 0;

for i = 1:d;

    bkg_voi = imagenes(masky1(7):masky2(7),maskx1(7):maskx2(7),I-2:I+2,i);

     data(i).sphere_A = data(1).sphere_A;
   data(i).sphere_vol  = data(1).sphere_vol;
    data(i).sphere_reference_time = data(1).sphere_reference_time;
    
%[CRCmax, CRCmean, SUVmax, SUVmean, CV, error, fail] = get_CRCs_islands(im_filt, ratio, maskx1, maskx2, masky1, masky2, I, bkg_voi);
[results(i).CRCmax, results(i).CRCmean, results(i).SUVmax, results(i).SUVmean, results(i).CV, results(i).error, results(i).fail] = get_CRCs_islands_EARL(imagenes(:,:,:,i), data(i), maskx1, maskx2, masky1, masky2, I, bkg_voi, EARL_mode, half_life);
error = error + results(i).error;
fails = fails + results(i).fail;
end


% if(CV>15/100)
%     fail = 999;
% end

% disp(FWHMs(1))
% disp(FWHMs(2))
 disp(['FWHM: ' num2str(FWHM)]);
 disp(['Fails: ' num2str(fails)]);
% 
% disp(error)
% disp(CRCmax)
% disp(CRCmean)
% disp(CV)
