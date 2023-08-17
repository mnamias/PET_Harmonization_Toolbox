%%
% Name of the code: harmonization_menu Toolbox
%
% Author: Mauro Nam�as (mnamias@gmail.com)
% Copyright (c) 2018-2023, Mauro Nam�as
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
%
%   Distributed under the terms of the "New BSD License."  Please see
%   license.txt.
%
% Overview: This is a GUI for the paper "A novel
% approach for quantitative harmonization_menu in PET" M Nam�as et al 2018 Phys. Med. Biol. 63 095019
%
% Inputs and ranges of inputs: the only input is a set of DICOM files of a
% cylindrical phantom.
%
% Outputs: Radial and axial LSFs, NPS and estimated CRCmax and CRCmean
% values for a NEMA phantom with a 10:1 sphere to background ratio.
%
% Version: 1.0.3-beta
%
% built-in functions:
%
% Level of the code: main function (GUI)
%
% Other codes related to the code: see callbacks
%
% Date of creation: April 17, 2018
%
% Modification history: 
% 1.0.0-beta: initial version.
% 1.0.1-beta: Added missing dependencies.
% 1.0.2-beta: Added missing dependencies.
% 1.0.2-beta: Changed main figure size to fit small displays.
% 1.0.3-beta: Updated user manual and copyright. 






function varargout = harmonization_toolbox(varargin)
% HARMONIZATION_TOOLBOX MATLAB code for harmonization_toolbox.fig
%      HARMONIZATION_TOOLBOX, by itself, creates a new HARMONIZATION_TOOLBOX or raises the existing
%      singleton*.
%
%      H = HARMONIZATION_TOOLBOX returns the handle to a new HARMONIZATION_TOOLBOX or the handle to
%      the existing singleton*.
%
%      HARMONIZATION_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HARMONIZATION_TOOLBOX.M with the given input arguments.
%
%      HARMONIZATION_TOOLBOX('Property','Value',...) creates a new HARMONIZATION_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before harmonization_toolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to harmonization_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help harmonization_toolbox

% Last Modified by GUIDE v2.5 08-Aug-2023 11:47:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @harmonization_toolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @harmonization_toolbox_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before harmonization_toolbox is made visible.
function harmonization_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to harmonization_toolbox (see VARARGIN)

% Choose default command line output for harmonization_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes harmonization_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);
disclaimer = 'THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.';

h = warndlg( {'PET Harmonization Toolbox v1.1 (August 16, 2023).' ' ' 'Copyright 2018-2023 Mauro Nam�as (mnamias@gmail.com)' ' ' 'This toolbox is to accompany the method published in PMB 2018 "A novel approach for quantitative harmonization in PET" by M. Nam�as et. al.' ' ' disclaimer} ); 
 
addpath([pwd filesep 'Functions']);
addpath([pwd filesep 'Functions_from_other_authors']);

% --- Outputs from this function are returned to the command line.
function varargout = harmonization_toolbox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function load_cyl_Callback(hObject, eventdata, handles)
% hObject    handle to load_cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image
global dicom_headers

% reads the desired directory containing dicom files
set(handles.text1,'String','Reading files...')
[image, dicom_headers] = read_directory();

thk = dicom_headers(1).SliceThickness;
thk2 = abs(dicom_headers(2).ImagePositionPatient(3)-dicom_headers(1).ImagePositionPatient(3));

dicom_headers(1).SliceThickness = min(thk,thk2);

showfig = 0;
results = check_phantom(image, dicom_headers, showfig);


mip = mipss(image);
axes(handles.axes1) 
imshow(-squeeze(image(:,round(end/2),:)),[-max(image(:))*1.5 0])
set(handles.text1,'String','Reading files... done!')
set(handles.text2,'String',['Date: ' dicom_headers(1).SeriesDate])
set(handles.text3,'String',['Time: ' dicom_headers(1).SeriesTime])
set(handles.text4,'String',['Series description: ' dicom_headers(1).SeriesDescription])

% Enable processing menu
set(handles.proc_menu,'Enable','on')



%set(h, 'position', [500 400 500 400]); %makes box bigger
h = msgbox({'Phantom Acquisition QA:' ['X-tilt: ' num2str(results.x_tilt) ' degrees']...
    ['X-offset: ' num2str(results.x_offset) ' mm']...
    ['Y-tilt: ' num2str(results.y_tilt) ' degrees']...
    ['Y-offset: ' num2str(results.y_offset) ' mm']...
    ['Z-offset: ' num2str(results.z_offset) ' mm']...
    });

ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set(ch, 'FontName', 'Century Gothic No11 L');
set(ch, 'FontSize', 10);

% --------------------------------------------------------------------
function proc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to proc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function axial_lsf_Callback(hObject, eventdata, handles)
% hObject    handle to axial_lsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global faxis_a
global MTF_a
global lsf_a
global esf_a
global image
global dicom_headers
global flag_lsf_a
global flag_lsf_r
global flag_NPS


load harmonization_options.mat options
options.mode = 'manual';

[faxis_a, MTF_a, lsf_a, esf_a] = LSFcalc_Z(image, dicom_headers, options);

% plot MTF_a and LSF_a
axes(handles.axes4);
loglog(faxis_a(faxis_a>=0),MTF_a(faxis_a>=0),'k-.')
title('Axial MTF')
xlabel('Spatial frequency [mm-1]')
ylabel('Amplitude');

axes(handles.axes5);
pix_size = dicom_headers(1).SliceThickness;
saxis = 1:length(lsf_a);
saxis = saxis*pix_size;
plot(saxis, lsf_a);
title('Axial LSF')
xlabel('Spatial position [mm]')
ylabel('Amplitude');

flag_lsf_a = 1;

if(flag_lsf_a & flag_lsf_r & flag_NPS)    
    set(handles.harmonization_menu,'Enable','on');
end



% --------------------------------------------------------------------
function radial_lsf_Callback(hObject, eventdata, handles)
% hObject    handle to radial_lsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global faxis_r
global MTF_r
global lsf_r
global esf_r
global image
global dicom_headers
global flag_lsf_a
global flag_lsf_r
global flag_NPS

load harmonization_options.mat options
options.mode = 'manual';
[faxis_r, MTF_r, lsf_r, esf_r] = LSFcalc_R(image, dicom_headers, options);

% plot MTF_r and LSF_r
axes(handles.axes2);
loglog(faxis_r(faxis_r>=0),MTF_r(faxis_r>=0),'k-.')
title('Radial MTF')
xlabel('Spatial frequency [mm-1]')
ylabel('Amplitude');

axes(handles.axes3);
pix_size = dicom_headers(1).PixelSpacing(1)/options.over_factor_R;
saxis = 1:length(lsf_r);
saxis = saxis*pix_size;
plot(saxis, lsf_r);
title('Radial LSF')
xlabel('Spatial position [mm]')
ylabel('Amplitude');

flag_lsf_r = 1;

if(flag_lsf_a & flag_lsf_r & flag_NPS)    
    set(handles.harmonization_menu,'Enable','on');
end

% --------------------------------------------------------------------
function nps_Callback(hObject, eventdata, handles)
% hObject    handle to nps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image
global dicom_headers
global Wavg
global psf_noise
global avg_std
global mean_bkg
global flag_lsf_a
global flag_lsf_r
global flag_NPS

% Estimate the NPS and the noise psf
options.mode = 'manual';
[Wavg, Savg, rbin, npsr, npsr_norm, avg_std, mean_bkg, psf_noise, xaxis, yaxis, zaxis] = NPScalc( image , dicom_headers , options);

            
        axes(handles.axes13);
        subplot(1,3,1);
        plot(xaxis,Wavg(:,ceil(end/2),ceil(end/2)),'k.-');
        xlabel('x-dir (mm^{-1})');
        ylabel('NPS_x ((kBq/ml)^2 mm^3)');
        subplot(1,3,2);
        plot(yaxis,Wavg(ceil(end/2),:,ceil(end/2)),'k.-');
        xlabel('y-dir (mm^{-1})');
        ylabel('NPS_y ((kBq/ml)^2 mm^3)');
        subplot(1,3,3);
        zvalue = squeeze(Wavg(ceil(end/2),ceil(end/2),:));
        plot(zaxis,zvalue,'k.-');
        xlabel('z-dir (mm^{-1})')
        ylabel('NPS_z ((kBq/ml)^2 mm^3)');
        
        flag_NPS = 1;

if(flag_lsf_a & flag_lsf_r & flag_NPS)    
    set(handles.harmonization_menu,'Enable','on');
end



% --- Executes when uipanel3 is resized.
function uipanel3_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disclaimer = 'THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.';

h = msgbox( {'PET Harmonization Toolbox v1.1 (August 16th, 2023).' ' ' 'Copyright 2018-2023 Mauro Nam�as (mnamias@gmail.com)' ' ' 'This toolbox is to accompany the method published in PMB 2018 "A novel approach for quantitative harmonization in PET" M Nam�as et al 2018 Phys. Med. Biol. 63 095019.' ' ' disclaimer} ); 
 

% --------------------------------------------------------------------
function harmonization_menu_Callback(hObject, eventdata, handles)
% hObject    handle to harmonization_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function sim_Callback(hObject, eventdata, handles)
% hObject    handle to sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global psf_noise
global lsf_a
global lsf_r
global sims
global tabulated_results

simulation



% --------------------------------------------------------------------
function CRCs_Callback(hObject, eventdata, ~)
% hObject    handle to CRCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CRCs


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Auto_Callback(hObject, eventdata, handles)
% hObject    handle to Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global faxis_a
global MTF_a
global lsf_a
global esf_a
global image
global dicom_headers
global faxis_r
global MTF_r
global lsf_r
global esf_r
global Wavg
global psf_noise
global avg_std
global mean_bkg
global sim_flag



load harmonization_options.mat options
options.mode = 'auto';

[Wavg, Savg, rbin, npsr, npsr_norm, avg_std, mean_bkg, psf_noise, xaxis, yaxis, zaxis] = NPScalc( image , dicom_headers, options );
[faxis_a, MTF_a, lsf_a, esf_a] = LSFcalc_Z(image, dicom_headers, options);
[faxis_r, MTF_r, lsf_r, esf_r] = LSFcalc_R(image, dicom_headers, options);

save temp.mat faxis_r MTF_r  lsf_r  esf_r faxis_a MTF_a lsf_a esf_a options dicom_headers

% plot MTF_r and LSF_r
    axes(handles.axes2);
    loglog(faxis_r(faxis_r>=0),MTF_r(faxis_r>=0),'k-.')
    title('Radial MTF')
    xlabel('Spatial frequency [mm-1]')
    ylabel('Amplitude');

    axes(handles.axes3);
    pix_size = dicom_headers(1).PixelSpacing(1)/options.over_factor_R;
    saxis = 1:length(lsf_r);
    saxis = saxis*pix_size;
    plot(saxis, lsf_r);
    
    xq = min(saxis):0.1:max(saxis)
    vq = interp1(saxis,lsf_r,xq)
    maxi = max(vq)
    indexes = find(vq>=maxi/2)
    FWHM_radial = length(indexes)*0.1
    
    title(['Radial LSF, FWHM = ' num2str(FWHM_radial)  ' mm' ])
    xlabel('Spatial position [mm]')
    ylabel('Amplitude');

  
    
% plot MTF_a and LSF_a
    axes(handles.axes4);
    loglog(faxis_a(faxis_a>=0),MTF_a(faxis_a>=0),'k-.')
    title('Axial MTF')
    xlabel('Spatial frequency [mm-1]')
    ylabel('Amplitude');

    
    axes(handles.axes5);
    dz = abs(dicom_headers(1).ImagePositionPatient(3) - dicom_headers(2).ImagePositionPatient(3));    
    pix_size = dz;
    
    saxis = 1:length(lsf_a);
    saxis = saxis*pix_size;
    plot(saxis, lsf_a);
    
    xq = min(saxis):0.1:max(saxis)
    vq = interp1(saxis,lsf_a,xq)
    maxi = max(vq)
    indexes = find(vq>=maxi/2)
    FWHM_axial = length(indexes)*0.1
    
    
    title(['Axial LSF, FWHM = ' num2str(FWHM_axial) ' mm' ])   
    xlabel('Spatial position [mm]')
    ylabel('Amplitude');

  % plot the NPS result        
        
        axes(handles.axes13);
        subplot(1,3,1);
        plot(xaxis,Wavg(:,ceil(end/2),ceil(end/2)),'k.-');
        xlabel('x-dir (mm^{-1})');
        ylabel('NPS_x ((Bq/ml)^2 mm^3)');
        subplot(1,3,2);
        plot(yaxis,Wavg(ceil(end/2),:,ceil(end/2)),'k.-');
        xlabel('y-dir (mm^{-1})');
        ylabel('NPS_y ((Bq/ml)^2 mm^3)');
        subplot(1,3,3);
        zvalue = squeeze(Wavg(ceil(end/2),ceil(end/2),:));
        plot(zaxis,zvalue,'k.-');
        xlabel('z-dir (mm^{-1})')
        ylabel('NPS_z ((Bq/ml)^2 mm^3)');

        sim_flag = 1;
        set(handles.harmonization_menu,'Enable','on')


% --------------------------------------------------------------------
function UNIF_Callback(hObject, eventdata, handles)
% hObject    handle to UNIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image
global dicom_headers

[NU, CV, SUV, slice_flags, positions] = cyl_nema_nu(image, dicom_headers)


% --------------------------------------------------------------------
function load_NEMA_Callback(hObject, eventdata, handles)
% hObject    handle to load_NEMA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image
global dicom_headers

% reads the desired directory containing dicom files
set(handles.text1,'String','Reading files...')
[image, dicom_headers] = read_directory();

thk = dicom_headers(1).SliceThickness;
thk2 = abs(dicom_headers(2).ImagePositionPatient(3)-dicom_headers(1).ImagePositionPatient(3));

dicom_headers(1).SliceThickness = min(thk,thk2);

showfig = 0;
results = check_phantom(image, dicom_headers, showfig);


mip = single(mip_z(image));
axes(handles.axes1) 
imshow(-mip,[-max(mip(:))*1.5 0])
set(handles.text1,'String','Reading files... done!')
set(handles.text2,'String',['Date: ' dicom_headers(1).SeriesDate])
set(handles.text3,'String',['Time: ' dicom_headers(1).SeriesTime])
set(handles.text4,'String',['Series description: ' dicom_headers(1).SeriesDescription])

% Enable processing menu
set(handles.proc_menu,'Enable','on')



% --------------------------------------------------------------------
function NEMA_harm_Callback(hObject, eventdata, handles)
% hObject    handle to NEMA_harm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global image
global dicom_headers
global CRCmax
global CRCmean
global CV
global EARL_mode 
global data
global maskx1
global maskx2
global masky1
global masky2
global I
global bkg_voi

ratio = 9.75;
data = struct()


data(1).sphere_A = dicom_headers(1).RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose/1e6; %% Sphere activity in MBq
data(1).sphere_vol = 1000; %% Sphere dilution volume in ml
data(1).sphere_reference_time = dicom_headers(1).RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime(1:6);
data(1).img_time = dicom_headers(1).SeriesTime;



prompt = {'Enter the volume of the stock solution used to filled the spheeres [ml]:'};
dlg_title = 'Volume';
num_lines = 1;
def = {num2str(data(1).sphere_vol)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
data(1).sphere_vol = str2num(cell2mat(answer));


prompt = {'Enter the activity used to filled the spheeres [MBq]:'};
dlg_title = 'Activity';
num_lines = 1;
def = {num2str(data(1).sphere_A)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
data(1).sphere_A = str2num(cell2mat(answer));

prompt = {'Enter the activity reference time used to filled the spheeres [hhmmss]:'};
dlg_title = 'Activity measurement time';
num_lines = 1;
def = {data(1).sphere_reference_time};
answer = inputdlg(prompt,dlg_title,num_lines,def);
data(1).sphere_reference_time = cell2mat(answer);

prompt = {'Enter the start of image acquisition [hhmmss]:'};
dlg_title = 'Image start time';
num_lines = 1;
def = {data(1).img_time};
answer = inputdlg(prompt,dlg_title,num_lines,def);
data(1).img_time  = cell2mat(answer);

pixel_sizes = [dicom_headers(1).PixelSpacing' dicom_headers(1).SliceThickness]
options = optimset('TolX',0.01, 'TolFun',0.00001,'Display','iter')
   
CRCmax_min = [0.31 0.59 0.73 0.83 0.91 0.95];
CRCmax_max = [0.49 0.85 1.01 1.09 1.13 1.16];
   
CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
CRCmean_max = [0.38 0.60 0.73 0.78 0.85 0.89];
diameters = [10 13 17 22 28 37];
   
       
im_per = permute(image(:,:,5:end-5),[3 2 1]);
im = mipss(im_per)';
im2 = image;
im2(:,:,1:5) = 0;
im2(:,:,end-15:end) = 0;
zp = squeeze(sum(sum(im2,1),2));
[M,I] = max(zp);

figure
imagesc(-(im-max(im(:))))
xlabel('pixel')
ylabel('pixel')
colormap gray
hold on
axis image;
impixelinfo;
            
fprintf('Select background ROI .\n\n');
title('Select Background ROI')

h = imrect;
ptmask = wait(h);
            
ptmask = wait(h);

maskx1 = round(ptmask(1));
maskx2 = maskx1 + round(ptmask(3));
masky1 = round(ptmask(2));
masky2 = masky1 + round(ptmask(4));

% force selected ROI to be within the bounds of the image
maskx1 = min(maskx1,size(im,2));
maskx2 = min(maskx2,size(im,2));
masky1 = min(masky1,size(im,1));
masky2 = min(masky2,size(im,1));
maskx1(7) = max(maskx1,1);
maskx2(7) = max(maskx2,1);
masky1(7) = max(masky1,1);
masky2(7) = max(masky2,1);

mask(masky1(7):masky2(7),maskx1(7):maskx2(7)) = 0;

plot([maskx1(7),maskx2(7),maskx2(7),maskx1(7),maskx1(7)],[masky1(7),masky1(7),masky2(7),masky2(7),masky1(7)],'r-')
            
bkg_voi = image(masky1(7):masky2(7),maskx1(7):maskx2(7),I-1:I+1);
            

            

SUVmax = zeros(6,1);
SUVmean = zeros(6,1);
CRCmax = zeros(6,1);
CRCmean = zeros(6,1);

%% Estimate Sphere CRCs
for i = 1:6
title(['Select Sphere: ' num2str(i)  '.\n\n'])
fprintf(['Select Sphere: ' num2str(i)  '.\n\n']);
h = imrect;
ptmask = wait(h);            
%ptmask = wait(h);            
maskx1(i) = round(ptmask(1));
maskx2(i) = maskx1(i) + round(ptmask(3));
masky1(i) = round(ptmask(2));
masky2(i) = masky1(i) + round(ptmask(4));

voi =  image(masky1(i):masky2(i),maskx1(i):maskx2(i),I-10:I+10);
SUVmax(i) = max(voi(:));

end

title('Estimating CRCs...');

t1 = data(1).sphere_reference_time;
t2 = data(1).img_time;
h1 = str2num(t1(1:2));
h2 = str2num(t2(1:2));
m1 = str2num(t1(3:4));
m2 = str2num(t2(3:4));
s1 = str2num(t1(5:6));
s2 = str2num(t2(5:6));

dt = (h2-h1)*60+(m2-m1)+(s2-s1)/60; %% Elapsed time in minutes
Fd = exp(-log(2)*dt/109.8); %% F18 decay factor
true_sphere_A = data(1).sphere_A/data(1).sphere_vol * Fd *1e6; %% expected sphere activity concentration, in Bq/ml
            
                        
bkg = mean(bkg_voi(:));
for i = 1:6
    voi =  image(masky1(i):masky2(i),maskx1(i):maskx2(i),I-10:I+10);
    maxi = SUVmax(i);
    thr = (maxi+bkg)/2;
    SUVmean(i) = mean(voi(voi>=thr));
    CRCmax(i) = SUVmax(i)/true_sphere_A;
    CRCmean(i) = SUVmean(i)/true_sphere_A;                  
end

CV = std(bkg_voi(:))/mean(bkg_voi(:));
    

save temp.mat CRCmax CRCmean CV dicom_headers data
% 
% light = 6;
% standard = 4;
% heavy = 2;
% none = 0;
% 
% filtro_z = 4;
% FWHM = 0;
% %filtro_z = none;
% %[FWHM fval] = fminsearch(@(FWHM) get_CRC_error_GE(image, ratio, maskx1, maskx2, masky1, masky2, I, bkg_voi, pixel_sizes, FWHM, filtro_z), 7.0,options)
% 
% brand = 'Siemens'
% half_life = 109.9; % F18
% 
% 
% EARL_mode = 1
%                         
% options.TolFun = 1e-2
% options.TolX = 0.05
% [FWHM fval] = fminsearch(@(FWHM) get_CRC_error_GE_full_islands_multi_EARL(imagees, data, maskx1, maskx2, masky1, masky2, I, bkg_voi, pixel_sizes, FWHM, filtro_z, brand, EARL_mode, half_life), 5.0,options)
% 
% [error, fails, results] = get_CRC_error_multibrand_multi(imagees, data, maskx1, maskx2, masky1, masky2, I, bkg_voi, pixel_sizes, FWHM, filtro_z,brand, EARL_mode, half_life);
%                         
% true_sphere_A/mean(bkg_voi(:))
%                         
% [cell2mat({results.fail}); cell2mat({results.CRCmax}) ;cell2mat({results.CRCmean});cell2mat({results.CV})]
%                        
% [cell2mat({results.fail}); cell2mat({results.SUVmax}) ;cell2mat({results.SUVmean});cell2mat({results.CV})]
% 
EARL_mode = 1
CRCs_NEMA

                       
