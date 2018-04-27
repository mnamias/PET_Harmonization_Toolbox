function im_filt = filtra_3D(image, pixel_sizes, FWHMs);
% function im_filt = filtra_3D(image, pixel_sizes, FWHMs);
% Copyright 2018 Mauro Namías
% function im_filt = filtra_3D(image, pixel_sizes, FWHMs);
% filters the input image with the corresponding FWHMs (1x3 vector)

sigmas = FWHMs/2.35;

sigmas_pix = sigmas./pixel_sizes;

sizes = ceil(6*sigmas_pix);
sizes = sizes + not(mod(sizes,2));

gaussKernel = fast_3D_Gaussian(sigmas_pix,sizes);
gaussKernel = padarray(gaussKernel, [1 1 1], 0, 'post');

im_filt = imfilter(image,gaussKernel);



