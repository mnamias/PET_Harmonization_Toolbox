function [sims tabulated_results] = simulate_spheres(lsf_r, lsf_a, psf_noise, options, Nmax, CV, dicom_headers)
%% Name of code: simulate_spheres.m
% Version: 1.0 (April 20, 2018)
% Level of code: function. This function is from the PET harmonization toolbox GUI
% Copyright (c) 2018, Mauro Namías mnamias@gmail.com
% All rights reserved.
% This code is intended to accompany the paper:
% Namías et. al
% A novel approach to quantitative harmonization in PET
% PMB (2018)
%
% function [sims tabulated_results] = simulate_spheres(lsf_r, lsf_a, psf_noise, options, Nmax, CV, dicom_headers)
%
% Simulates the spheres from the NEMA phantom based on resolution and noise
% measurements.
%
% Inputs: 
%   lsf_r: radial line spread function
%   lsf_a: axial line spread function
%   psf_noise: noise psf
%   options: simulation options (MATLAB structure)
%   dicom_headers: MATLAB structure with the dicom headers for each slice
%   Nmax: number of noise realizations
%   CV: desired background CV of the NEMA phantom%   
%
% Options: processing options
%   options.over_factor_R : Oversampling factor used for LSF_R analysis
%
% Outputs:
%   sims: (Nmax x 6) MATLAB structure containing the simulated spheres
%   tabulated results: CRC summery statistics in tabular form
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
load spheres.mat
mean_bkg = 1;
over_factor = options.over_factor_R;
lsf_z = double(lsf_a);
lsf = double(lsf_r);

pixel_sizes = [dicom_headers(1).PixelSpacing' dicom_headers(1).SliceThickness];

dummy = zeros(1,1,length(lsf_z));
dummy(1,1,:) = lsf_z;
lsf_z = dummy;

FWHMs = [1 1 0];
CRCmax_a = []
CRCmean_a = []
CVs = []
avg_std = mean_bkg*(CV/100);
sims = [];
progressbar()


    for ii = 1:Nmax
        progressbar(ii/Nmax)
        %% Model Multiplicative Noise 
        for s = 1:6
            dx = round((rand*pixel_sizes(1))-pixel_sizes(1)/2);
            dy = round((rand*pixel_sizes(1))-pixel_sizes(1)/2);
            dz = round((rand*pixel_sizes(3))-pixel_sizes(3)/2);
    
            x_factor = pixel_sizes(1)/over_factor;
            z_factor = pixel_sizes(3);

            spheres(s).sphere_res = resample3Dimage_aniso(circshift(spheres(s).sphere_f,[dx dy dz]), x_factor, x_factor, z_factor); %% resampled phantom        
            spheres(s).sphere_blur =  imfilter(spheres(s).sphere_res, fliplr(lsf));
            spheres(s).sphere_blur =  imfilter(spheres(s).sphere_blur,lsf');

            dx = round((rand*over_factor)-over_factor/2);
            dy = round((rand*over_factor)-over_factor/2);
        
            spheres(s).sphere_blur_rp =  imresize(circshift(spheres(s).sphere_blur,[0 0 0]),1/over_factor,'Antialiasing',1);
            spheres(s).sphere_blur_rp =  imfilter(spheres(s).sphere_blur_rp, lsf_z);
    
            noise = single(randn(size(spheres(s).sphere_blur_rp)));
            noise = convolution3D_FFTdomain(single(noise),single(psf_noise));
            
            p = avg_std^2/mean_bkg;
            nn_vst = real(noise./std(noise(:)).* (spheres(s).sphere_blur_rp*p).^(1/2));

            spheres(s).sphere_sim = spheres(s).sphere_blur_rp-nn_vst;
            sims(ii,s).sphere = spheres(s).sphere_sim;
            
        end
    end
    
    
            sims(1).bkg_noise = noise./std(noise(:)) * p^(1/2);

progressbar(1)

filtro_z = 0;
FWHM = 0;
sphere_6 = spheres(6).sphere_f;
sphere_signal = max(sphere_6(:));

showfig = 1;

[error results] = get_CRC_error_multibrand(sims, sphere_signal, mean_bkg, pixel_sizes, FWHM, filtro_z,'GE', showfig);

 %[cell2mat({results.CRCmax'}) ;cell2mat({results.CRCmean'})]
 
 tabulated_results = [ [mean(results.CRCmax,1)';mean(results.CRCmean,1)'], [min(results.CRCmax,[],1)';min(results.CRCmean,[],1)'] , [max(results.CRCmax,[],1)';max(results.CRCmean,[],1)'],[std(results.CRCmax,[],1)';std(results.CRCmean,[],1)'] ]
