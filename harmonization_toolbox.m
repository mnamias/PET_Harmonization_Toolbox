%%
% Name of the code: harmonization_menu Toolbox
%
% Author: Mauro Namías (mnamias@gmail.com)
% Copyright (c) 2018, Mauro Namías
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
% approach for quantitative harmonization_menu in PET" by M. Namías et. al., PMB, May 2018
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

% Last Modified by GUIDE v2.5 23-Apr-2018 11:02:09

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

h = warndlg( {'PET Harmonization Toolbox v1.0 (April 23, 2018).' ' ' 'Copyright 2018 Mauro Namías (mnamias@gmail.com)' ' ' 'This toolbox is to accompany the method published in PMB 2018 "A novel approach for quantitative harmonization in PET" by M. Namías et. al.' ' ' disclaimer} ); 
 
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

        % plot the result
        axes(handles.axes16);
        psf_mip = mipss(psf_noise);
        imshow(-psf_mip,[-max(psf_mip(:)) 0])
        
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

h = msgbox( {'PET Harmonization Toolbox v1.0 (April 23, 2018).' ' ' 'Copyright 2018 Mauro Namías (mnamias@gmail.com)' ' ' 'This toolbox is to accompany the method published in PMB 2018 "A novel approach for quantitative harmonization in PET" by M. Namías et. al.' ' ' disclaimer} ); 
 

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

  % plot the NPS result
        axes(handles.axes16);
        psf_mip = mipss(psf_noise);
        imshow(-psf_mip,[-max(psf_mip(:)) 0])
        
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
