function varargout = CRCs(varargin)
%%
% Name of the code: CRCs.m
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
% Overview: This is a sub-figure GUI for the paper "A novel
% approach for quantitative harmonization_menu in PET" by M. Namías et. al., PMB, April 2018



% CRCS MATLAB code for CRCs.fig
%      CRCS, by itself, creates a new CRCS or raises the existing
%      singleton*.
%
%      H = CRCS returns the handle to a new CRCS or the handle to
%      the existing singleton*.
%
%      CRCS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CRCS.M with the given input arguments.
%
%      CRCS('Property','Value',...) creates a new CRCS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CRCs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CRCs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CRCs

% Last Modified by GUIDE v2.5 23-Apr-2018 10:41:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CRCs_OpeningFcn, ...
                   'gui_OutputFcn',  @CRCs_OutputFcn, ...
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


% --- Executes just before CRCs is made visible.
function CRCs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CRCs (see VARARGIN)

% Choose default command line output for CRCs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CRCs wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global sims
global tabulated_results

tableData_max = get(handles.uitable3, 'data')

tableData_max(:,2:end) = tabulated_results(1:6,:);

tableData_mean = get(handles.uitable4, 'data');

tableData_mean(:,2:end) = tabulated_results(7:end,:);

set(handles.uitable3,'data',tableData_max);
set(handles.uitable4,'data',tableData_mean);

n_sims = size(sims,1)
set(handles.slider1,'min',1);
set(handles.slider1,'max',n_sims);
set(handles.slider1,'SliderStep',[1/(n_sims-1) ,1/(n_sims-1)*5]);
set(handles.slider1,'Value',1);


sphere_size = size(sims(1,6).sphere);
collage = zeros(sphere_size(1) * 2, sphere_size(1) * 3);


collage(1:end/2, 1:end/3) = squeeze(sims(1,1).sphere(:,:,round(end/2)));
collage(1:end/2, end/3+1:end*2/3) = squeeze(sims(1,2).sphere(:,:,round(end/2)));
collage(1:end/2, end*2/3+1:end) = squeeze(sims(1,3).sphere(:,:,round(end/2)));

collage(end/2 + 1:end, 1:end/3) = squeeze(sims(1,4).sphere(:,:,round(end/2)));
collage(end/2 + 1:end, end/3 + 1:end*2/3) = squeeze(sims(1,5).sphere(:,:,round(end/2)));
collage(end/2 + 1:end, end*2/3 + 1:end) = squeeze(sims(1,6).sphere(:,:,round(end/2)));

axes(handles.axes1);
imshow(-collage, [-15 1]);

bkg_CV = std(sims(1).bkg_noise(:));
set(handles.edit4, 'String', num2str(bkg_CV*100));

set(handles.edit7, 'String', 'N/A');






% --- Outputs from this function are returned to the command line.
function varargout = CRCs_OutputFcn(hObject, eventdata, handles) 
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
global sims
global filtered_sims
global tabulated_results
global dicom_headers
global filtered_flag



FWHM = get(handles.edit1,'String');
FWHM = str2num(FWHM)

flag2D = get(handles.radiobutton2,'Value');
z_filter = 0;


if(flag2D)
    
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

load spheres.mat
sphere_6 = spheres(6).sphere_f;
sphere_signal = max(sphere_6(:));
mean_bkg = 1;

showfig = 1;

if(flag2D)
    [error results filtered_sims] = get_CRC_error_multibrand(sims, sphere_signal, mean_bkg, pixel_sizes, FWHM, z_filter,'GE', showfig);    
else    
    [error results filtered_sims] = get_CRC_error_multibrand(sims, sphere_signal, mean_bkg, pixel_sizes, FWHM, z_filter,'Siemens', showfig);
end

 h = msgbox(['Number of CRCs outside EARL limits: ' num2str(results.fail)])

 set(handles.edit7,'String', num2str(results.fail));
 
tabulated_results = [ [mean(results.CRCmax,1)';mean(results.CRCmean,1)'], [min(results.CRCmax,[],1)';min(results.CRCmean,[],1)'] , [max(results.CRCmax,[],1)';max(results.CRCmean,[],1)'],[std(results.CRCmax,[],1)';std(results.CRCmean,[],1)'] ];

tableData_max = get(handles.uitable3, 'data');
tableData_max(:,2:end) = tabulated_results(1:6,:);
tableData_mean = get(handles.uitable4, 'data');
tableData_mean(:,2:end) = tabulated_results(7:end,:);
set(handles.uitable3,'data',tableData_max);
set(handles.uitable4,'data',tableData_mean);
    
filtered_flag = 1;
    
bkg_CV = std(filtered_sims(1).bkg_noise(:));

set(handles.edit4, 'String', num2str(bkg_CV*100));

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sims
global filtered_sims
global tabulated_results
global dicom_headers
global filtered_flag

flag2D = get(handles.radiobutton2,'Value');

z_filter = 0;
brand = 'Siemens';

if(flag2D)
    brand = 'GE';
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

load spheres.mat
sphere_6 = spheres(6).sphere_f;
sphere_signal = max(sphere_6(:));
mean_bkg = 1;

options_opt = optimset('TolX',0.01, 'TolFun',0.00001,'Display','iter')
showfig = 0;

h = msgbox('Finding optimal filter, please wait... ');
[FWHM fval] = fminsearch(@(FWHM) get_CRC_error_multibrand(sims, sphere_signal, mean_bkg, pixel_sizes, FWHM, z_filter, brand, showfig), 6.0, options_opt);
close(h)

set(handles.edit1,'String',num2str(FWHM));

showfig = 1;
[error results filtered_sims] = get_CRC_error_multibrand(sims, sphere_signal, mean_bkg, pixel_sizes, FWHM, z_filter,brand, showfig);    

 h = msgbox(['Number of CRCs outside EARL limits: ' num2str(results.fail)])
 set(handles.edit7,'String', num2str(results.fail));
 
 
tabulated_results = [ [mean(results.CRCmax,1)';mean(results.CRCmean,1)'], [min(results.CRCmax,[],1)';min(results.CRCmean,[],1)'] , [max(results.CRCmax,[],1)';max(results.CRCmean,[],1)'],[std(results.CRCmax,[],1)';std(results.CRCmean,[],1)'] ];

tableData_max = get(handles.uitable3, 'data');
tableData_max(:,2:end) = tabulated_results(1:6,:);
tableData_mean = get(handles.uitable4, 'data');
tableData_mean(:,2:end) = tabulated_results(7:end,:);
set(handles.uitable3,'data',tableData_max);
set(handles.uitable4,'data',tableData_mean);
    
filtered_flag = 1;

    
bkg_CV = std(filtered_sims(1).bkg_noise(:));

set(handles.edit4, 'String', num2str(bkg_CV*100));

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
