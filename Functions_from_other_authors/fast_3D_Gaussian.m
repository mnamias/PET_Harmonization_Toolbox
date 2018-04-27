function gaussKernel = fast_3D_Gaussian(sigmas,size)
% Creates a 3D Gaussian Kernel

if(size(3)>1)

%tic
for x = 1:size(1)
for y=1:size(2)
for z=1:size(3)
gaussKernel(y, x, z) = exp(- ( (x-( (size(1)+1) /2))^2/(2*sigmas(1)^2) + (y-( (size(2)+1) /2))^2/(2*sigmas(2)^2) + (z-( (size(3)+1) /2))^2/(2*sigmas(3)^2) ) );
end
end
end

else
    
z = 1;
for x = 1:size(1)
for y=1:size(2)
%gaussKernel(y, x, z) = exp(- ( (x-( (size(1)) /2))^2/(2*sigmas(1)^2) + (y-( (size(2)) /2))^2/(2*sigmas(2)^2) + (z-( (size(3)) /2))^2/(2*sigmas(3)^2) ) );
gaussKernel(y, x, z) = exp(- ( (x-( (size(1)+1) /2))^2/(2*sigmas(1)^2) + (y-( (size(2)+1) /2))^2/(2*sigmas(2)^2)  ));
end
end

end
gaussKernel = gaussKernel/sum(gaussKernel(:));

