function mip = mipss(image);
% function mip = mipss(image);
% generates a MIP of the image

mip = max(image,[],1);
mip = squeeze(mip);

