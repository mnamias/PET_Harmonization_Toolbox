function varargout = CRCs_NEMA(varargin)
%%
% Name of the code: CRCs_NEMA.m
%
% Author: Mauro Nam�as (mnamias@gmail.com)
% Copyright (c) 2023, Mauro Nam�as
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
% Overview: This is a sub-figure GUI for the paper "A novel
% approach for quantitative harmonization_menu in PET" by M. Nam�as et. al., PMB, April 2018



% CRCS_NEMA MATLAB code for CRCs_NEMA.fig
%      CRCS_NEMA, by itself, creates a new CRCS_NEMA or raises the existing
%      singleton*.
%
%      H = CRCS_NEMA returns the handle to a new CRCS_NEMA or the handle to
%      the existing singleton*.
%
%      CRCS_NEMA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRCS_NEMA.M with the given input arguments.
%
%      CRCS_NEMA('Property','Value',...) creates a new CRCS_NEMA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CRCs_NEMA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CRCs_NEMA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CRCs_NEMA

% Last Modified by GUIDE v2.5 08-Aug-2023 14:09:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CRCs_NEMA_OpeningFcn, ...
                   'gui_OutputFcn',  @CRCs_NEMA_OutputFcn, ...
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


% --- Executes just before CRCs_NEMA is made visible.
function CRCs_NEMA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CRCs_NEMA (see VARARGIN)

% Choose default command line output for CRCs_NEMA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CRCs_NEMA wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global CRCmax
global CRCmean
global CV
global image

diameters = [10 13 17 22 28 37]
data = [diameters', CRCmax, CRCmean]


set(handles.uitable3,'data',data);

set(handles.edit9,'String',num2str(CV*100));

mip = single(mip_z(image));
axes(handles.axes1) 
imshow(-mip,[-max(mip(:))*1.5 0])



CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
CRCmax_max = [0.57 0.85 1.01 1.09 1.13 1.16];
CRCmax_target = (CRCmax_min+CRCmax_max)/2;

CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
CRCmean_target = (CRCmean_min+CRCmean_max)/2;

axes(handles.axes3) 
plot(diameters, CRCmax, 'k o')
hold on
plot(diameters(1:6),CRCmax_max,'k --')
hold on
plot(diameters(1:6),CRCmax_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmax')



axes(handles.axes4) 
plot(diameters, CRCmean, 'k o')
hold on
plot(diameters(1:6),CRCmean_max,'k --')
hold on
plot(diameters(1:6),CRCmean_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmean')

% axes(handles.axes1);
% imshow(-collage, [-15 1]);


% 
% tableData_max = get(handles.uitable3, 'data');
% 
% tableData_max(:,2:end) = tabulated_results(1:6,:);
% 
% tableData_mean = get(handles.uitable4, 'data');
% 
% tableData_mean(:,2:end) = tabulated_results(7:end,:);
% 
% set(handles.uitable3,'data',tableData_max);
% set(handles.uitable4,'data',tableData_mean);
% 
% n_sims = size(sims,1);
% set(handles.slider1,'min',1);
% set(handles.slider1,'max',n_sims);
% set(handles.slider1,'SliderStep',[1/(n_sims-1) ,1/(n_sims-1)*5]);
% set(handles.slider1,'Value',1);
% 
% 
% sphere_size = size(sims(1,6).sphere);
% collage = zeros(sphere_size(1) * 2, sphere_size(1) * 3);
% 
% 
% collage(1:end/2, 1:end/3) = squeeze(sims(1,1).sphere(:,:,round(end/2)));
% collage(1:end/2, end/3+1:end*2/3) = squeeze(sims(1,2).sphere(:,:,round(end/2)));
% collage(1:end/2, end*2/3+1:end) = squeeze(sims(1,3).sphere(:,:,round(end/2)));
% 
% collage(end/2 + 1:end, 1:end/3) = squeeze(sims(1,4).sphere(:,:,round(end/2)));
% collage(end/2 + 1:end, end/3 + 1:end*2/3) = squeeze(sims(1,5).sphere(:,:,round(end/2)));
% collage(end/2 + 1:end, end*2/3 + 1:end) = squeeze(sims(1,6).sphere(:,:,round(end/2)));
% 
% axes(handles.axes1);
% imshow(-collage, [-15 1]);
% 
% bkg_CV = std(sims(1).bkg_noise(:));
% set(handles.edit4, 'String', num2str(bkg_CV*100));
% 
% set(handles.edit7, 'String', 'N/A');
% 





% --- Outputs from this function are returned to the command line.
function varargout = CRCs_NEMA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
set(handles.radiobutton2,'Value',1)
set(handles.radiobutton3,'Value',0)
set(handles.radiobutton4,'Enable','on')
set(handles.radiobutton5,'Enable','on')
set(handles.radiobutton6,'Enable','on')
set(handles.radiobutton7,'Enable','on')

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image
global dicom_headers
global CRCmax
global CRCmean
global CV
global data
global EARL_mode 
global maskx1
global maskx2
global masky1
global masky2
global I
global bkg_voi


FWHM = get(handles.edit1,'String');
FWHM = str2num(FWHM)

flag2D = get(handles.radiobutton2,'Value');
z_filter = 0;

half_life = 109.8;

brand = 'Siemens'

if(flag2D)
    brand = 'GE'
    flag_light = get(handles.radiobutton5,'Value');
    if(flag_light)
        z_filter = 6;
    end
    
    flag_std = get(handles.radiobutton6,'Value');
    if(flag_std)
        z_filter = 4;
    end
    
    flag_heavy = get(handles.radiobutton7,'Value');
    if(flag_heavy)
        z_filter = 2;
    end
end

pixel_sizes = [dicom_headers(1).PixelSpacing' dicom_headers(1).SliceThickness];

showfig = 1;

imagenes = image;

[error, fails, results, imagenes] = get_CRC_error_multibrand_multi(imagenes, data, maskx1, maskx2, masky1, masky2, I, bkg_voi, pixel_sizes, FWHM,  z_filter, brand, EARL_mode, half_life);

h = msgbox(['Number of CRCs outside EARL limits: ' num2str(results.fail)]);

tabulated_results = [results.CRCmax,results.CRCmean]

CRCmax = results.CRCmax
CRCmean = results.CRCmean

CV = results.CV

diameters = [10 13 17 22 28 37]
data1 = [diameters', tabulated_results]

set(handles.uitable3,'data',data1);

%set(handles.edit9,'String',num2str(CV*100));

mip = single(mip_z(imagenes(:,:,:,1)));
axes(handles.axes1) 
imshow(-mip,[-max(mip(:))*1.5 0])




%tableData = get(handles.uitable3, 'data');
%tableData(:,2:end) = tabulated_results(1:6,:);
%set(handles.uitable3,'data',tableData);
   
filtered_flag = 1;
    
save temp.mat results

set(handles.edit9, 'String', num2str(CV*100));




diameters = [10 13 17 22 28 37]

 if(EARL_mode == 1)
   CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
   CRCmax_max = [0.57 0.85 1.01 1.09 1.13 1.16];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
   CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
    elseif (EARL_mode == 2)
        
   CRCmax_min = [0.52 0.85 1.00 1.01 1.01 1.05];
   CRCmax_max = [0.88 1.22 1.38 1.32 1.26 1.29];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.39 0.63 0.76 0.80 0.82 0.85];
   CRCmean_max = [0.61 0.86 0.97 0.99 0.97 1.00];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
 end
    
axes(handles.axes3) 
cla
plot(diameters, CRCmax, 'k o')
hold on
plot(diameters(1:6),CRCmax_max,'k --')
hold on
plot(diameters(1:6),CRCmax_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmax')
axis([5 40 min(CRCmax_min)*0.9 max(CRCmax_max)*1.1])

axes(handles.axes4) 
cla
plot(diameters, CRCmean, 'k o')
hold on
plot(diameters(1:6),CRCmean_max,'k --')
hold on
plot(diameters(1:6),CRCmean_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmean')
axis([5 40 min(CRCmean_min)*0.9 max(CRCmean_max)*1.1])
   



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global image
global dicom_headers
global CRCmax
global CRCmean
global CV
global data
global EARL_mode 
global maskx1
global maskx2
global masky1
global masky2
global I
global bkg_voi

flag_2 = get(handles.radiobutton15,'Value')

if(flag_2)
    EARL_mode = 2
end


FWHM = get(handles.edit1,'String');
FWHM = str2num(FWHM)

flag2D = get(handles.radiobutton2,'Value');
z_filter = 0;

half_life = 109.8;

brand = 'Siemens'

if(flag2D)
    brand = 'GE'
    flag_light = get(handles.radiobutton5,'Value');
    if(flag_light)
        z_filter = 6;
    end
    
    flag_std = get(handles.radiobutton6,'Value');
    if(flag_std)
        z_filter = 4;
    end
    
    flag_heavy = get(handles.radiobutton7,'Value');
    if(flag_heavy)
        z_filter = 2;
    end
end

pixel_sizes = [dicom_headers(1).PixelSpacing' dicom_headers(1).SliceThickness];

showfig = 1;

imagenes = image;

options.TolFun = 1e-2
options.TolX = 0.05

[FWHM fval] = fminsearch(@(FWHM) get_CRC_error_GE_full_islands_multi_EARL(imagenes, data, maskx1, maskx2, masky1, masky2, I, bkg_voi, pixel_sizes, FWHM, z_filter, brand, EARL_mode, half_life), 5.0,options)
                  
set(handles.edit1,'String', num2str(FWHM));


[error, fails, results, imagenes] = get_CRC_error_multibrand_multi(imagenes, data, maskx1, maskx2, masky1, masky2, I, bkg_voi, pixel_sizes, FWHM,  z_filter, brand, EARL_mode, half_life);

h = msgbox(['Number of CRCs outside EARL limits: ' num2str(results.fail)]);

tabulated_results = [results.CRCmax,results.CRCmean]

CRCmax = results.CRCmax
CRCmean = results.CRCmean

CV = results.CV

diameters = [10 13 17 22 28 37]
data1 = [diameters', tabulated_results]

set(handles.uitable3,'data',data1);

%set(handles.edit9,'String',num2str(CV*100));

mip = single(mip_z(imagenes(:,:,:,1)));
axes(handles.axes1) 
imshow(-mip,[-max(mip(:))*1.5 0])




%tableData = get(handles.uitable3, 'data');
%tableData(:,2:end) = tabulated_results(1:6,:);
%set(handles.uitable3,'data',tableData);
   
filtered_flag = 1;
    
save temp.mat results

set(handles.edit9, 'String', num2str(CV*100));



diameters = [10 13 17 22 28 37]

 if(EARL_mode == 1)
   CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
   CRCmax_max = [0.57 0.85 1.01 1.09 1.13 1.16];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
   CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
    elseif (EARL_mode == 2)
        
   CRCmax_min = [0.52 0.85 1.00 1.01 1.01 1.05];
   CRCmax_max = [0.88 1.22 1.38 1.32 1.26 1.29];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.39 0.63 0.76 0.80 0.82 0.85];
   CRCmean_max = [0.61 0.86 0.97 0.99 0.97 1.00];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
 end
    
axes(handles.axes3) 
cla
plot(diameters, CRCmax, 'k o')
hold on
plot(diameters(1:6),CRCmax_max,'k --')
hold on
plot(diameters(1:6),CRCmax_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmax')
axis([5 40 min(CRCmax_min)*0.9 max(CRCmax_max)*1.1])

axes(handles.axes4) 
cla
plot(diameters, CRCmean, 'k o')
hold on
plot(diameters(1:6),CRCmean_max,'k --')
hold on
plot(diameters(1:6),CRCmean_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmean')
axis([5 40 min(CRCmean_min)*0.9 max(CRCmean_max)*1.1])
   




% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
set(handles.radiobutton2,'Value',0)
set(handles.radiobutton3,'Value',1)

set(handles.radiobutton4,'Enable','off')
set(handles.radiobutton5,'Enable','off')
set(handles.radiobutton6,'Enable','off')
set(handles.radiobutton7,'Enable','off')




% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
set(handles.radiobutton4,'Value',1)
set(handles.radiobutton5,'Value',0)
set(handles.radiobutton6,'Value',0)
set(handles.radiobutton7,'Value',0)

% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
set(handles.radiobutton4,'Value',0)
set(handles.radiobutton5,'Value',1)
set(handles.radiobutton6,'Value',0)
set(handles.radiobutton7,'Value',0)

% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
set(handles.radiobutton4,'Value',0)
set(handles.radiobutton5,'Value',0)
set(handles.radiobutton6,'Value',1)
set(handles.radiobutton7,'Value',0)

% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7
set(handles.radiobutton5,'Value',0)
set(handles.radiobutton6,'Value',0)
set(handles.radiobutton4,'Value',0)
set(handles.radiobutton7,'Value',1)

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global sims
global filtered_sims
global filtered_flag

% Display selected noise realization

index = round(get(handles.slider1,'Value'));

if(filtered_flag)
    selected_sims = filtered_sims;
else
    selected_sims = sims;
end


sphere_size = size(sims(1,6).sphere);
collage = zeros(sphere_size(1) * 2, sphere_size(1) * 3);

collage(1:end/2, 1:end/3) = squeeze(selected_sims(index,1).sphere(:,:,round(end/2)));
collage(1:end/2, end/3+1:end*2/3) = squeeze(selected_sims(index,2).sphere(:,:,round(end/2)));
collage(1:end/2, end*2/3+1:end) = squeeze(selected_sims(index,3).sphere(:,:,round(end/2)));

collage(end/2 + 1:end, 1:end/3) = squeeze(selected_sims(index,4).sphere(:,:,round(end/2)));
collage(end/2 + 1:end, end/3 + 1:end*2/3) = squeeze(selected_sims(index,5).sphere(:,:,round(end/2)));
collage(end/2 + 1:end, end*2/3 + 1:end) = squeeze(selected_sims(index,6).sphere(:,:,round(end/2)));

axes(handles.axes1);
imshow(-collage, [-15 1]);

set(handles.edit3,'String', num2str(index));

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
global sims
global filtered_sims
global filtered_flag

if(filtered_flag)
    selected_sims = filtered_sims;
else
    selected_sims = sims;
end

current_index = get(handles.slider1,'Value')
index = round(str2num(get(handles.edit3,'String')))


if(index <= size(sims,1) & index > 1)
    set(handles.slider1,'Value',index);
    
    sphere_size = size(sims(1,6).sphere);
    collage = zeros(sphere_size(1) * 2, sphere_size(1) * 3);

    collage(1:end/2, 1:end/3) = squeeze(selected_sims(index,1).sphere(:,:,round(end/2)));
    collage(1:end/2, end/3+1:end*2/3) = squeeze(selected_sims(index,2).sphere(:,:,round(end/2)));
    collage(1:end/2, end*2/3+1:end) = squeeze(selected_sims(index,3).sphere(:,:,round(end/2)));

    collage(end/2 + 1:end, 1:end/3) = squeeze(selected_sims(index,4).sphere(:,:,round(end/2)));
    collage(end/2 + 1:end, end/3 + 1:end*2/3) = squeeze(selected_sims(index,5).sphere(:,:,round(end/2)));
    collage(end/2 + 1:end, end*2/3 + 1:end) = squeeze(selected_sims(index,6).sphere(:,:,round(end/2)));

    axes(handles.axes1);
    imshow(-collage, [-15 1]);
    
else
    set(handles.edit3,'String',num2str(current_index));
end

    




% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton14.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton14
global CRCmax
global CRCmean
global EARL_mode

EARL_mode = 1


set(handles.radiobutton14,'Value',1)
set(handles.radiobutton15,'Value',0)


diameters = [10 13 17 22 28 37]

 if(EARL_mode == 1)
   CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
   CRCmax_max = [0.57 0.85 1.01 1.09 1.13 1.16];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
   CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
    elseif (EARL_mode == 2)
        
   CRCmax_min = [0.52 0.85 1.00 1.01 1.01 1.05];
   CRCmax_max = [0.88 1.22 1.38 1.32 1.26 1.29];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.39 0.63 0.76 0.80 0.82 0.85];
   CRCmean_max = [0.61 0.86 0.97 0.99 0.97 1.00];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
 end
    
axes(handles.axes3) 
cla
plot(diameters, CRCmax, 'k o')
hold on
plot(diameters(1:6),CRCmax_max,'k --')
hold on
plot(diameters(1:6),CRCmax_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmax')
axis([5 40 min(CRCmax_min)*0.9 max(CRCmax_max)*1.1])

axes(handles.axes4) 
cla
plot(diameters, CRCmean, 'k o')
hold on
plot(diameters(1:6),CRCmean_max,'k --')
hold on
plot(diameters(1:6),CRCmean_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmean')
axis([5 40 min(CRCmean_min)*0.9 max(CRCmean_max)*1.1])
   
   

% --- Executes on button press in radiobutton15.
function radiobutton15_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CRCmax
global CRCmean
global EARL_mode

EARL_mode = 2


set(handles.radiobutton14,'Value',0)
set(handles.radiobutton15,'Value',1)

diameters = [10 13 17 22 28 37]

 if(EARL_mode == 1)
   CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
   CRCmax_max = [0.57 0.85 1.01 1.09 1.13 1.16];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
   CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
    elseif (EARL_mode == 2)
        
   CRCmax_min = [0.52 0.85 1.00 1.01 1.01 1.05];
   CRCmax_max = [0.88 1.22 1.38 1.32 1.26 1.29];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.39 0.63 0.76 0.80 0.82 0.85];
   CRCmean_max = [0.61 0.86 0.97 0.99 0.97 1.00];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
 end
    
axes(handles.axes3) 
cla
plot(diameters, CRCmax, 'k o')
hold on
plot(diameters(1:6),CRCmax_max,'k --')
hold on
plot(diameters(1:6),CRCmax_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmax')
axis([5 40 min(CRCmax_min)*0.9 max(CRCmax_max)*1.1])

axes(handles.axes4) 
cla
plot(diameters, CRCmean, 'k o')
hold on
plot(diameters(1:6),CRCmean_max,'k --')
hold on
plot(diameters(1:6),CRCmean_min,'k --')
xlabel('Sphere diameter [mm]')
%ylabel('CRCmean')
axis([5 40 min(CRCmean_min)*0.9 max(CRCmean_max)*1.1])
   
    

% Hint: get(hObject,'Value') returns toggle state of radiobutton15
