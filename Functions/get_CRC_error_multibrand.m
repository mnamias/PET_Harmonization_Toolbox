function [error results spheres] = get_CRC_error_multibrand(spheres, sphere_signal, mean_bkg, pixel_sizes, FWHM, zfilter, brand, showfig);
% Name of code: get_CRC_error_multibrand.m
% Version: 1.0 (April 20, 2018)
% Level of code: function. This function is from the PET harmonization toolbox GUI
% Copyright (c) 2018, Mauro Namías mnamias@gmail.com
% All rights reserved.
% This code is intended to accompany the paper:
% Namías et. al
% A novel approach to quantitative harmonization in PET
% PMB (2018)
%
% [error results spheres] = get_CRC_error_multibrand(spheres, sphere_signal, mean_bkg, pixel_sizes, FWHM, zfilter, brand, showfig);
%
% Estimates error between the simulated spheres CRC´s and the target CRC values. 
%
% Inputs: 
%   spheres: (Nmax x 6) MATLAB structure with the simulated spheres 
%   sphere_signal: ideal sphere uptake (default: 9.75 for EARL harmonization).
%   mean_bkg: background signal (default: 1)
%   pixel_sizes: 1x3 array with pixel sizes in mm
%   FWHM: FWHM of the Gaussian filter 
%   zfilter: weight of the axial moving average filter (0:none 2:heavy 4:std or 6:light)
%   brand: filter type: 'GE' --> use 2D filter + axial filter. 'SIEMENS'--> 3D isotropic filter
%   showfig: if equal to 1, will show figures. 
%
% Outputs:
% error: RMS error between measured CRCs and target (EARL) CRCs. 
% results: structure with CRC statistics
% spheres: (Nmax x 6) structure with filtered spheres.
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

[a,b] = size(spheres)

if(FWHM>15)
    FWHM = 15-randn;
end

if(strcmp(brand,'GE'))
    
    if(FWHM>0)

        for i = 1:a
            for j = 1:b
                spheres(i,j).sphere = filtra_3D(spheres(i,j).sphere, pixel_sizes, [FWHM FWHM 0]);
            end
        end
        
        spheres(1).bkg_noise = filtra_3D(spheres(1).bkg_noise, pixel_sizes, [FWHM FWHM 0]);

    end
    

else
      
    if(FWHM>0)
        for i = 1:a
            for j = 1:b
                spheres(i,j).sphere = filtra_3D(spheres(i,j).sphere, pixel_sizes, [FWHM FWHM FWHM]);
            end
        end
        
        spheres(1).bkg_noise = filtra_3D(spheres(1).bkg_noise, pixel_sizes, [FWHM FWHM FWHM]);

    end

end


if(zfilter >0)
    z_kernel = zeros(1,1,3);
    z_kernel(:,:,1) = 1;
    z_kernel(:,:,2) = zfilter;
    z_kernel(:,:,3) = 1;
    z_kernel = z_kernel/sum(z_kernel(:));

    for i = 1:a
        for j = 1:b
            spheres(i,j).sphere = imfilter(spheres(i,j).sphere ,z_kernel);
        end
    end
    spheres(1).bkg_noise = imfilter(spheres(1).bkg_noise ,z_kernel);
end

results = [];

   
[results.CRCmax, results.CRCmean, results.error, results.fail] =  get_CRCs_quick(spheres, mean_bkg, sphere_signal, showfig);

error = results.error;

 disp(['FWHM: ' num2str(FWHM)]);
 disp(['Fails: ' num2str(results.fail)]);

