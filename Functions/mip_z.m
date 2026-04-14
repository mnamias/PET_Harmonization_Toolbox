function mip = mip_z(image);
% function mip = mipss(image);
% generates a MIP of the image

mip = max(image,[],3);
mip = squeeze(mip);

