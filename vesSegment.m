function varargout = vesSegment(varargin)
% VESSEGMENT MATLAB code for vesSegment.fig
%      VESSEGMENT, by itself, creates a new VESSEGMENT or raises the existing
%      singleton*.
%
%      H = VESSEGMENT returns the handle to a new VESSEGMENT or the handle to
%      the existing singleton*.
%
%      VESSEGMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VESSEGMENT.M with the given input arguments.
%
%      VESSEGMENT('Property','Value',...) creates a new VESSEGMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vesSegment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vesSegment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vesSegment



% Last Modified by GUIDE v2.5 10-Aug-2017 14:57:48


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vesSegment_OpeningFcn, ...
                   'gui_OutputFcn',  @vesSegment_OutputFcn, ...
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


% --- Executes just before vesSegment is made visible.
function vesSegment_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vesSegment (see VARARGIN)

% Choose default command line output for vesSegment
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Data
if exist('Data')
    if isfield(Data,'Graph')
        Data = rmfield(Data,'Graph');
    end
end

% --- Outputs from this function are returned to the command line.
function varargout = vesSegment_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to File_loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
% clear all
[filename,pathname] = uigetfile({'*.mat;*.tiff;*.tif'},'Please select the Angiogram Data');
h = waitbar(0,'Please wait... loading the data');
[~,~,ext] = fileparts(filename);
if strcmp(ext,'.mat')
    load([pathname filename]);
    if exist('Output','var')
        if isfield(Data,'angio')
            if ~strcmp(Data.rawdatapath,Output.rawdatapath)
                error('Output raw data path and currently loaded raw data path did not match');
            end
        else
            [~,~,ext] = fileparts(Output.rawdatapath);
            if strcmp(ext,'.mat')
                if ~exist(Output.rawdatapath,'file')
                    error('Raw data was moved from Original location');
                end
                temp = load(Output.rawdatapath);
                fn = fieldnames(temp);
                Data.angio = temp.(fn{1});
            elseif strcmp(ext,'.tiff') || strcmp(ext,'.tif')
                info = imfinfo(Output.rawdatapath);
                for u = 1:length(info)
                    if u == 1
                        temp = imread(Output.rawdatapath,1);
                        angio = zeros([length(info) size(temp)]);
                        angio(u,:,:) = temp;
                    else
                        angio(u,:,:) = imread(Output.rawdatapath,u);
                    end
                end
                Data.angio = angio;
            end
            Data.rawdatapath = Output.rawdatapath;
        end
        if isfield(Output,'procSteps')
            Data.procSteps = Output.procSteps;
        end
        if isfield(Output,'angioF')
            Data.angioF = Output.angioF;
        end
        if isfield(Output,'angioT')
            Data.angioT = Output.angioT;
        end
        if isfield(Output,'segangio')
            Data.segangio = Output.segangio;
        end
        if isfield(Output,'fv')
            Data.fv = Output.fv;
        end
    else
        temp = load([pathname filename]);
        fn = fieldnames(temp);
        Data.angio = temp.(fn{1});
        Data.rawdatapath = [pathname filename];
    end
elseif strcmp(ext,'.tiff') || strcmp(ext,'.tif')
    info = imfinfo([pathname filename]);
    for u = 1:length(info)
        if u == 1
            temp = imread([pathname filename],1);
            angio = zeros([length(info) size(temp)]);
            angio(u,:,:) = temp;
        else
            angio(u,:,:) = imread([pathname filename],u);
        end
    end
    Data.angio = angio;
    Data.rawdatapath = [pathname filename];
end
[z,x,y] = size(Data.angio);
set(handles.edit_Zstartframe,'String',num2str(1));
set(handles.edit_ZMIP,'String',num2str(z));
% set(handles.edit_XcenterZoom,'String',num2str(1));
% set(handles.edit_XwidthZoom,'String',num2str(x));
% set(handles.edit_YcenterZoom,'String',num2str(1));
% set(handles.edit_YwidthZoom,'String',num2str(y));

set(handles.edit_XcenterZoom,'String',num2str(mean([1 size(Data.angio,2)-1])))
set(handles.edit_YcenterZoom,'String',num2str(mean([1 size(Data.angio,3)-1])))
set(handles.edit_XwidthZoom,'String',num2str(size(Data.angio,2)));
set(handles.edit_YwidthZoom,'String',num2str(size(Data.angio,2)));

% set(handles.edit_imageInfo,'String','Image info');
% set(handles.edit_imageInfo,'String',[num2str(x) 'X' num2str(y) 'X' num2str(z)]);
str = sprintf('%s\n%s','Image info',[num2str(x) 'X' num2str(y) 'X' num2str(z)]);
set(handles.edit_imageInfo,'String',str);
waitbar(1);
close(h);

draw(hObject, eventdata, handles);


function draw(hObject, eventdata, handles)

global Data

if (isfield(Data,'angioF') || isfield(Data,'angioT')) && (get(handles.checkbox_showFiltData,'Value') == 1)
    if isfield(Data,'angioT')
        I = Data.angioT;
    else
        I = Data.angioF;
    end
elseif isfield(Data,'segangio') && (get(handles.checkbox_showFiltData,'Value') == 0) && (get(handles.checkbox_showRawData,'Value') == 0) && (get(handles.checkbox_showSeg,'Value') == 1)
    I = Data.segangio*0;
else
    I = Data.angio;
end

[Sz,Sx,Sy] = size(I);

Zstartframe = str2double(get(handles.edit_Zstartframe,'String'));
Zstartframe = min(max(Zstartframe,1),Sz);
ZMIP = str2double(get(handles.edit_ZMIP,'String'));
Zendframe = min(max(Zstartframe+ZMIP-1,1),Sz);
Xstartframe = str2double(get(handles.edit_XcenterZoom,'String'));
Xstartframe = min(max(Xstartframe,1),Sx);
XMIP = str2double(get(handles.edit_XwidthZoom,'String'));
Xendframe = min(max(Xstartframe+XMIP-1,1),Sx);
Ystartframe = str2double(get(handles.edit_YcenterZoom,'String'));
Ystartframe = min(max(Ystartframe,1),Sy);
YMIP = str2double(get(handles.edit_YwidthZoom,'String'));
Yendframe = min(max(Ystartframe+YMIP-1,1),Sy);

Zimg = squeeze(max(I(Zstartframe:Zendframe,:,:),[],1));
%Ximg = squeeze(max(I(:,Xstartframe:Xendframe,:),[],2));
%Yimg = squeeze(max(I(:,:,Ystartframe:Yendframe),[],3));
%ZMIPimg = squeeze(max(I,[],1));

if isfield(Data,'segangio') && (get(handles.checkbox_showSeg,'Value') == 1)
    ZimgS = squeeze(max(Data.segangio(Zstartframe:Zendframe,:,:),[],1));
%    XimgS = squeeze(max(Data.segangio(:,Xstartframe:Xendframe,:),[],2));
%    YimgS = squeeze(max(Data.segangio(:,:,Ystartframe:Yendframe),[],3));
%    imgS = squeeze(max(Data.segangio,[],1));
end

axes(handles.axes1)
if (get(handles.checkbox_showFiltData,'Value') == 1) || (get(handles.checkbox_showRawData,'Value') == 1)
colormap('gray')
imagesc(Zimg)
if isfield(Data,'segangio') && (get(handles.checkbox_showSeg,'Value') == 1)
    hold on
    img = double(ZimgS);
    green = cat(3, zeros(size(img)),ones(size(img)), zeros(size(img)));
    hold on
    h = imagesc(green);
    hold off
    set(h, 'AlphaData', img*0.25)
else
    if get(handles.checkbox_showColorbar,'Value') == 1
        colorbar;
    end
end
axis image
axis on
if isfield(Data,'ZoomXrange') && isfield(Data,'ZoomYrange')
    xlim(Data.ZoomXrange);
    ylim(Data.ZoomYrange);
end
elseif get(handles.checkbox_showSeg,'Value') == 1 && isfield(Data,'segangio')
%     C = [0 0 0; 0 0.25 0; 0 0.5 0; 0 1 0];
    img = double(ZimgS);
    green = cat(3, zeros(size(img)),img*0.25, zeros(size(img)));
    imagesc(green);
%     colormap(C);
%     colorbar; 
    axis image; 
    axis on
    if isfield(Data,'ZoomXrange') && isfield(Data,'ZoomYrange')
        xlim(Data.ZoomXrange);
        ylim(Data.ZoomYrange);
    end
%     set(h, 'AlphaData', img*0.25)
end

% Display Graph
if get(handles.checkboxDisplayGraph,'value')==1
    nodes = Data.Graph.nodes;
    edges = Data.Graph.edges;
    lst = find(nodes(:,1)>=Data.ZoomXrange(1) & nodes(:,1)<=Data.ZoomXrange(2) & ...
               nodes(:,2)>=Data.ZoomYrange(1) & nodes(:,2)<=Data.ZoomYrange(2) & ...
               nodes(:,3)>=Zstartframe & nodes(:,3)<=Zendframe );
    hold on
    h=plot(nodes(lst,1),nodes(lst,2),'m.');
    set(h,'markersize',12)
    
    for ii=1:length(lst)
        lst2 = find(edges(:,1)==lst(ii));
        plot(nodes(edges(lst2,:),1), nodes(edges(lst2,:),2), 'm-' );
    end
    hold off
end

% Display procSteps in text box
if isfield(Data,'procSteps')
    for u = 1:size(Data.procSteps,1)
        if u == 1
            str = sprintf('%s',string(Data.procSteps(u,1)));
        else
            str = sprintf('%s\n%s',str,string(Data.procSteps(u,1)));
        end
    end
    set(handles.edit_dispProcSteps,'String',str);
else
    str = '';
    set(handles.edit_dispProcSteps,'String',str);
end

% axes(handles.axes2)
% colormap('gray')
% imagesc(Ximg)
% if isfield(Data,'angioT') && (get(handles.checkbox_showSeg,'Value') == 1)
%     hold on
%     img = double(XimgS);
%     green = cat(3, zeros(size(img)),ones(size(img)), zeros(size(img)));
%     hold on
%     h = imshow(green);
%     hold off
%     set(h, 'AlphaData', img*0.25)
% end
% axis image
% 
% axes(handles.axes3)
% colormap('gray')
% imagesc(Yimg)
% if isfield(Data,'angioT') && (get(handles.checkbox_showSeg,'Value') == 1)
%     hold on
%     img = double(YimgS);
%     green = cat(3, zeros(size(img)),ones(size(img)), zeros(size(img)));
%     hold on
%     h = imshow(green);
%     hold off
%     set(h, 'AlphaData', img*0.25)
% end
% axis image
% 
% axes(handles.axes4)
% colormap('gray')
% imagesc(ZMIPimg)
% if isfield(Data,'angioT') && (get(handles.checkbox_showSeg,'Value') == 1)
%     hold on
%     img = double(imgS);
%     green = cat(3, zeros(size(img)),ones(size(img)), zeros(size(img)));
%     hold on
%     h = imshow(green);
%     hold off
%     set(h, 'AlphaData', img*0.25)
% end
% axis image


function edit_Zstartframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Zstartframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Zstartframe as text
%        str2double(get(hObject,'String')) returns contents of edit_Zstartframe as a double
global Data

[z,x,y] = size(Data.angio);
ii = get(handles.edit_Zstartframe,'Value');
ii = min(max(ii,1),z);
set(handles.edit_Zstartframe,'Value',ii);
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_Zstartframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Zstartframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_XcenterZoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XcenterZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XcenterZoom as text
%        str2double(get(hObject,'String')) returns contents of edit_XcenterZoom as a double
global Data
Xcenter = str2double(get(handles.edit_XcenterZoom,'String'));
XwidthZoom = str2double(get(handles.edit_XwidthZoom,'String'));
Data.ZoomXrange = [max((Xcenter-XwidthZoom/2),1) min((Xcenter+XwidthZoom/2),size(Data.angio,2))];
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_XcenterZoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_XcenterZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_YcenterZoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_YcenterZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_YcenterZoom as text
%        str2double(get(hObject,'String')) returns contents of edit_YcenterZoom as a double
global Data
Ycenter = str2double(get(handles.edit_YcenterZoom,'String'));
YwidthZoom = str2double(get(handles.edit_YwidthZoom,'String'));
Data.ZoomYrange = [max((Ycenter-YwidthZoom/2),1) min((Ycenter+YwidthZoom/2),size(Data.angio,2))];
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_YcenterZoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_YcenterZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ZMIP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ZMIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ZMIP as text
%        str2double(get(hObject,'String')) returns contents of edit_ZMIP as a double
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_ZMIP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ZMIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_XwidthZoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XwidthZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XwidthZoom as text
%        str2double(get(hObject,'String')) returns contents of edit_XwidthZoom as a double
global Data
Xcenter = str2double(get(handles.edit_XcenterZoom,'String'));
XwidthZoom = str2double(get(handles.edit_XwidthZoom,'String'));
Data.ZoomXrange = [max((Xcenter-XwidthZoom/2),1) min((Xcenter+XwidthZoom/2),size(Data.angio,2))];

draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_XwidthZoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_XwidthZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_YwidthZoom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_YwidthZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_YwidthZoom as text
%        str2double(get(hObject,'String')) returns contents of edit_YwidthZoom as a double
global Data
Ycenter = str2double(get(handles.edit_YcenterZoom,'String'));
YwidthZoom = str2double(get(handles.edit_YwidthZoom,'String'));
Data.ZoomYrange = [max((Ycenter-YwidthZoom/2),1) min((Ycenter+YwidthZoom/2),size(Data.angio,2))];

draw(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit_YwidthZoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_YwidthZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Zmoveleft.
function pushbutton_Zmoveleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Zmoveleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ii = str2double(get(handles.edit_Zstartframe,'String'));
ii = ii-1;
set(handles.edit_Zstartframe,'String',num2str(ii));
draw(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_Xmoveleft.
function pushbutton_Xmoveleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Xmoveleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
ii = str2double(get(handles.edit_XcenterZoom,'String'));
XwidthZoom = str2double(get(handles.edit_XwidthZoom,'String'));
ii = ii-round(XwidthZoom*.10);
ii = min(max(ii,round(XwidthZoom/2)),size(Data.angio,2)-round(XwidthZoom/2));
set(handles.edit_XcenterZoom,'String',num2str(ii));
Xcenter = str2double(get(handles.edit_XcenterZoom,'String'));
XwidthZoom = str2double(get(handles.edit_XwidthZoom,'String'));
Data.ZoomXrange = [max((Xcenter-XwidthZoom/2),1) min((Xcenter+XwidthZoom/2),size(Data.angio,2))];
draw(hObject, eventdata, handles);

% --- Executes on button press in pushbutton_Ymoveleft.
function pushbutton_Ymoveleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Ymoveleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
ii = str2double(get(handles.edit_YcenterZoom,'String'));
YwidthZoom = str2double(get(handles.edit_YwidthZoom,'String'));
ii = ii-round(YwidthZoom*.10);
ii = min(max(ii,round(YwidthZoom/2)),size(Data.angio,3)-round(YwidthZoom/2));
set(handles.edit_YcenterZoom,'String',num2str(ii));
Ycenter = str2double(get(handles.edit_YcenterZoom,'String'));
Data.ZoomYrange = [max((Ycenter-YwidthZoom/2),1) min((Ycenter+YwidthZoom/2),size(Data.angio,2))];
draw(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_Zmoveright.
function pushbutton_Zmoveright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Zmoveright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ii = str2double(get(handles.edit_Zstartframe,'String'));
ii = ii+1;
set(handles.edit_Zstartframe,'String',num2str(ii));
draw(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_Xmoveright.
function pushbutton_Xmoveright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Xmoveright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
ii = str2double(get(handles.edit_XcenterZoom,'String'));
XwidthZoom = str2double(get(handles.edit_XwidthZoom,'String'));
ii = ii+round(XwidthZoom*.10);
ii = min(max(ii,round(XwidthZoom/2)),size(Data.angio,2)-round(XwidthZoom/2));
set(handles.edit_XcenterZoom,'String',num2str(ii));
Xcenter = str2double(get(handles.edit_XcenterZoom,'String'));
XwidthZoom = str2double(get(handles.edit_XwidthZoom,'String'));
Data.ZoomXrange = [max((Xcenter-XwidthZoom/2),1) min((Xcenter+XwidthZoom/2),size(Data.angio,2))];
draw(hObject, eventdata, handles);

% --- Executes on button press in pushbutton_Ymoveright.
function pushbutton_Ymoveright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Ymoveright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
ii = str2double(get(handles.edit_YcenterZoom,'String'));
YwidthZoom = str2double(get(handles.edit_YwidthZoom,'String'));
ii = ii+round(YwidthZoom*.10);
ii = min(max(ii,round(YwidthZoom/2)),size(Data.angio,3)-round(YwidthZoom/2));
set(handles.edit_YcenterZoom,'String',num2str(ii));
Ycenter = str2double(get(handles.edit_YcenterZoom,'String'));
Data.ZoomYrange = [max((Ycenter-YwidthZoom/2),1) min((Ycenter+YwidthZoom/2),size(Data.angio,2))];
draw(hObject, eventdata, handles);


% --------------------------------------------------------------------
function Filters_Callback(hObject, eventdata, handles)
% hObject    handle to Filters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Filters_GausFilter_Callback(hObject, eventdata, handles)
% hObject    handle to Filters_GausFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
prompt = {'Please enter sigma for Gaussian Filter'};
defaultans = {'2'};
x = inputdlg(prompt,'Gaussian Filter',1,defaultans);
sigma = str2double(x{1});
if isfield(Data,'angioF')
    I = Data.angioF;
else
    I = Data.angio;
end
h = waitbar(0,'Please wait... applying gaussian filter');
Data.angioF = imgaussfilt3(I,sigma);
if isfield(Data,'procSteps')
    Data.procSteps(end+1,:) =  {{'Gaussian Filter'},{'Sigma'},{sigma}};
else
    Data.procSteps = {{'Gaussian Filter'},{'Sigma'},{sigma}};
end
waitbar(1);
close(h);
draw(hObject, eventdata, handles);

% --------------------------------------------------------------------
function Filters_MedFilter_Callback(hObject, eventdata, handles)
% hObject    handle to Filters_MedFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
prompt = {'Enter filter size in Z :','Enter filter size in X :','Enter filter size in Y :'};
defaultans = {'3','3','3'};
x = inputdlg(prompt,'Filter Size',1,defaultans );
sz = str2double(x{1});
sx = str2double(x{2});
sy = str2double(x{3});
if isfield(Data,'angioF')
    I = Data.angioF;
else
    I = Data.angio;
end
h = waitbar(0,'Please wait... applying Median filter');
Data.angioF = medfilt3(I,[sz sx sy]);
if isfield(Data,'procSteps')
    Data.procSteps(end+1,:) =  {{'Median Filter'},{'Size'},{[sz sx sy]}};
else
    Data.procSteps =  {{'Median Filter'},{'Size'},{[sz sx sy]}};
end
waitbar(1);
close(h);
draw(hObject, eventdata, handles);


% --------------------------------------------------------------------
function Filters_TubenessFilter_Callback(hObject, eventdata, handles)
% hObject    handle to Filters_TubenessFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
h = waitbar(0,'Please wait... applying Tubeness filter');

if isfield(Data,'angioF')
    I = Data.angioF;
else
    I = Data.angio;
end
I = double(I);
% Z = size(I,3);
% beta = 100;
% c = 500;
[k,l,m] = size(I);
% [nz,nx,ny] = size(I);
% L   = zeros(k,l,m);
% Vs  = zeros(k,l,m);
alpha = 0.25;
gamma12 = 0.5;
gamma23 = 0.5;
T = zeros(k,l,m);
prompt = {'Enter Gaussian filter start value :','Enter Gaussian filter end value :','Enter Gaussian filter step size :'};
defaultans = {'2','3','1'};
x = inputdlg(prompt,'Tubeness filter parameters',1,defaultans );
sigma = str2double(x{1}):str2double(x{3}):str2double(x{2});
for i = 1:length(sigma)
    
    waitbar(i-1/length(sigma));
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigma(i));
    
%     Normalizing the Hessian Matrix
    Dxx = sigma(i)^2*Dxx; Dyy = sigma(i)^2*Dyy;  Dzz = sigma(i)^2*Dzz; Dxy = sigma(i)^2*Dxy;  Dxz = sigma(i)^2*Dxz; Dyz = sigma(i)^2*Dyz;
    
    
    
    [Lambda1,Lambda2,Lambda3,~,~,~] = eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);

    
    SortL = sort([Lambda1(:)'; Lambda2(:)'; Lambda3(:)'],'descend');
    Lambda1 = reshape(SortL(1,:),size(Lambda1));
    Lambda2 = reshape(SortL(2,:),size(Lambda2));
    Lambda3 = reshape(SortL(3,:),size(Lambda3));
    
    idx = find(Lambda3 < 0 & Lambda2 < 0 & Lambda1 < 0);
    T(idx ) = abs(Lambda3(idx)).*(Lambda2(idx)./Lambda3(idx)).^gamma23.*(1+Lambda1(idx)./abs(Lambda2(idx))).^gamma12;
    idx = find(Lambda3 < 0 & Lambda2 < 0 & Lambda1 > 0 & Lambda1 < abs(Lambda2)/alpha);
    T(idx ) = abs(Lambda3(idx)).*(Lambda2(idx)./Lambda3(idx)).^gamma23.*(1-alpha*Lambda1(idx)./abs(Lambda2(idx))).^gamma12;
    
    %         L1 = (2*Lambda1-Lambda2-Lambda3)./(2*sqrt(Lambda1.^2+Lambda2.^2+Lambda3.^2-Lambda1.*Lambda2-Lambda1.*Lambda3-Lambda2.*Lambda3));
    %         L1 = exp(-alpha*(L1-1).^2);
    %         L1(abs(Lambda1)> abs(Lambda2)) = 0;
    %         L1(Lambda2>0 | Lambda3>0) = 0;
    %         L1 = -L1.*Lambda2;
    %         L = max(L,L1);
    %
    %         Ra = abs(Lambda2./Lambda3);
    %         Rb = abs(Lambda1./(Lambda2.*Lambda3));
    %         s = sqrt(Lambda1.^2+Lambda2.^2+Lambda3.^2);
    %         Vs1 = 1-exp(-Ra.^2/(2*alpha)).*exp(-Rb.^2/(2*beta)).*(1-exp(-s.^2/(2*c)));
    %         Vs1(Lambda2>0 | Lambda3>0) = 0;
    %         Vs1(abs(Lambda1) > abs(Lambda2)) = 0;
    %         Vs = max(Vs,Vs1);
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(T,sigma(i));
    % Normalizing the Hessian Matrix
     Dxx = sigma(i)^2*Dxx; Dyy = sigma(i)^2*Dyy;  Dzz = sigma(i)^2*Dzz; Dxy = sigma(i)^2*Dxy;  Dxz = sigma(i)^2*Dxz; Dyz = sigma(i)^2*Dyz;
    
    [Lambda1,Lambda2,Lambda3,V1,V2,V3] = eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    
    SortL = sort([Lambda1(:)'; Lambda2(:)'; Lambda3(:)'],'ascend');
    Lambda1 = reshape(SortL(1,:),size(Lambda1));
    Lambda2 = reshape(SortL(2,:),size(Lambda2));
    Lambda3 = reshape(SortL(3,:),size(Lambda3));
    
    E = -sigma(i)^2.*Lambda2;
    E(E<0) = 0;
    if i == 1
        Emax = E;
    else
        Emax = max(E,Emax);
    end
end
%     T = L;
T = Emax;

T = (T-min(T(:)))/(max(T(:))-min(T(:)));
Data.angioT = T;
if isfield(Data,'procSteps')
    Data.procSteps(end+1,:) =  {{'Tubeness Filter'},{'Sigma'},{sigma}};
else
    Data.procSteps =  {{'Tubeness Filter'},{'Sigma'},{sigma}};
end
close(h);
draw(hObject, eventdata, handles);

% % --------------------------------------------------------------------
% function Filters_thresholding_Callback(hObject, eventdata, handles)
% % hObject    handle to Filters_thresholding (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% global Data
% % Check if angioT exists before segmentation
% if isfield(Data,'angioT')
%     prompt = {'Please enter threshold value for segmentation'};
%     defaultans = {'0.075'};
%     x = inputdlg(prompt,'Segmenatation',1,defaultans);
%     threshold = str2double(x{1});
%     T = Data.angioT;
%     T_seg = zeros(size(T));
%     idx = find(T > threshold);
%     T_seg(idx) = 1;
%     CC = bwconncomp(T_seg);
%     T_segM = T_seg;
%     for uuu = 1:length(CC.PixelIdxList)
%          if length(CC.PixelIdxList{uuu}) < 100
%              T_segM(CC.PixelIdxList{uuu}) = 0;
%          end
%     end 
%     Data.segangio = T_segM;
% end


% --- Executes on button press in checkbox_showRawData.
function checkbox_showRawData_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showRawData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showRawData

ii = get(handles.checkbox_showRawData,'Value');
if ii == 1
    set(handles.checkbox_showFiltData,'Value',0);
% else
%      set(handles.checkbox_showFiltData,'Value',1);
end
draw(hObject, eventdata, handles);

% --- Executes on button press in checkbox_showFiltData.
function checkbox_showFiltData_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showFiltData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showFiltData

ii = get(handles.checkbox_showFiltData,'Value');
if ii == 1
    set(handles.checkbox_showRawData,'Value',0);
% else
%     set(handles.checkbox_showRawData,'Value',1);
end
draw(hObject, eventdata, handles);


% --- Executes on button press in checkbox_showSeg.
function checkbox_showSeg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showSeg
draw(hObject, eventdata, handles);


% --------------------------------------------------------------------
function Filters_ResetFilter_Callback(hObject, eventdata, handles)
% hObject    handle to Filters_ResetFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data

if isfield(Data,'angioF')
    Data = rmfield(Data,'angioF');
end

if isfield(Data,'angioT')
    Data = rmfield(Data,'angioT');
end

if isfield(Data,'segangio')
    Data = rmfield(Data,'segangio');
end

if isfield(Data,'procSteps')
    Data = rmfield(Data,'procSteps');
end
draw(hObject, eventdata, handles);

% --------------------------------------------------------------------
function File_savedata_Callback(hObject, eventdata, handles)
% hObject    handle to File_savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = waitbar(0,'Please wait... saving the data');
global Data
if isfield(Data,'angioF')|| isfield(Data,'angioT') || isfield(Data,'segangio')
    if isfield(Data,'rawdatapath')
        Output.rawdatapath = Data.rawdatapath;
    end
    if isfield(Data,'procSteps')
        Output.procSteps = Data.procSteps;
    end
    if isfield(Data,'angioF')
        Output.angioF = Data.angioF;
    end
    if isfield(Data,'angioT')
        Output.angioT = Data.angioT;
    end
    if isfield(Data,'segangio')
        Output.segangio = Data.segangio;
    end
    if isfield(Data,'fv')
        Output.fv = Data.fv;
    end
    [FileName,PathName] = uiputfile('*.mat');
    save([PathName FileName],'Output');
end
waitbar(1);
close(h);
    



% --------------------------------------------------------------------
function Segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to Segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function segmentation_thresholding_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation_thresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = waitbar(0,'Please wait... saving the data');
global Data
% Check if angioT exists before segmentation
if isfield(Data,'angioT')
    prompt = {'Please enter threshold value for segmentation'};
    defaultans = {'0.05'};
    x = inputdlg(prompt,'Segmenatation',1,defaultans);
    threshold = str2double(x{1});
    T = Data.angioT;
    T_seg = zeros(size(T));
    idx = find(T > threshold);
    T_seg(idx) = 1;
    CC = bwconncomp(T_seg);
    T_segM = T_seg;
    for uuu = 1:length(CC.PixelIdxList)
         if length(CC.PixelIdxList{uuu}) < 100
             T_segM(CC.PixelIdxList{uuu}) = 0;
         end
    end 
    Data.segangio = T_segM;
    if isfield(Data,'procSteps')
        Data.procSteps(end+1,:) =  {{'Thresholding on tubeness filter'},{'Threshold value'},{threshold}};
    else
        Data.procSteps =  {{'Thresholding on tubeness filter'},{'Threshold value'},{threshold}};
    end
end
waitbar(1);
close(h);
draw(hObject, eventdata, handles);



function edit_dispProcSteps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dispProcSteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dispProcSteps as text
%        str2double(get(hObject,'String')) returns contents of edit_dispProcSteps as a double


% --- Executes during object creation, after setting all properties.
function edit_dispProcSteps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dispProcSteps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function segmentation_SeedBasedSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation_SeedBasedSegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
addpath(genpath([pwd '\seed_based_segmentation']));
if isfield(Data,'segangio')
    options.fg_seed_vol = Data.segangio;
    [seg_vol, seg_prob, fg_seed_vol, bg_seed_vol] = segment_vessels_random_walker(Data.angio, options)
    
end



% --- Executes on button press in pushbutton_Zoomin.
function pushbutton_Zoomin_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
axes(handles.axes1)
title('\fontsize{16}\color{red}ZOOM IN');
k = waitforbuttonpress;
if k==0
    point1 = get(handles.axes1,'CurrentPoint');     % button down detected
    finalRect = rbbox;                              % return figure units
    point2 = get(handles.axes1,'CurrentPoint');     % button up detected
    point1 = round(point1(1,1:2));                          % extract x and y
    point2 = round(point2(1,1:2));
  
%     if isfield(cna,'ROIs')
%         cna.ROIs(end+1,:,:,:) = [point1; point2; [zrange(1) zrange(end)]];
%     else
%         cna.ROIs(1,:,:,:) = [point1; point2;[zrange(1) zrange(end)]];
%     end
%     draw(hObject, eventdata, handles);
Data.ZoomYrange = [point1(2) point2(2)];
Data.ZoomXrange = [point1(1) point2(1)];
set(handles.edit_XcenterZoom,'String',num2str(mean([point1(1) point2(1)-1])));
set(handles.edit_YcenterZoom,'String',num2str(mean([point1(2) point2(2)-1])));
set(handles.edit_XwidthZoom,'String',num2str(point2(1)-point1(1)+1));
set(handles.edit_YwidthZoom,'String',num2str(point2(2)-point1(2)+1));
end
% rect_pos = rbbox;
draw(hObject, eventdata, handles);



% --- Executes on button press in pushbutton_ZoomOut.
function pushbutton_ZoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
Data.ZoomYrange = [1 size(Data.angio,3)];
Data.ZoomXrange = [1 size(Data.angio,2)];
set(handles.edit_XcenterZoom,'String',num2str(mean([1 size(Data.angio,2)-1])))
set(handles.edit_YcenterZoom,'String',num2str(mean([1 size(Data.angio,3)-1])))
set(handles.edit_XwidthZoom,'String',num2str(size(Data.angio,2)));
set(handles.edit_YwidthZoom,'String',num2str(size(Data.angio,2)));
draw(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_displayMesh.
function pushbutton_displayMesh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_displayMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
wait_h = waitbar(0,'Please wait... calculating the mesh');
if isfield(Data,'segangio')
    Mask = permute(Data.segangio,[2 3 1]);
    if isfield(Data,'fv')
        fv2 = Data.fv;
    else
        fv = isosurface(Mask);
        fv2 = reducepatch(fv,200000); 
    end
    
%     fv2 = fv;
    f = fv2.faces;
    v = fv2.vertices;
    figure(2);
    clf

    h=trisurf(f,v(:,1),v(:,2),v(:,3),'facecolor','red','edgecolor','none');
    daspect([1,1,1])
    view(3); axis tight
    camlight
    lighting gouraud
    xlabel('Y')
    ylabel('X')
    zlabel('Z')
    Data.fv = fv2;
end
offset = [Xstartframe,Ystartframe,Zstartframe];
save('mesh.mat','Mask','f','v','offset');
waitbar(1);
close(wait_h);



function edit_imageInfo_Callback(hObject, eventdata, handles)
% hObject    handle to edit_imageInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_imageInfo as text
%        str2double(get(hObject,'String')) returns contents of edit_imageInfo as a double


% --- Executes during object creation, after setting all properties.
function edit_imageInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_imageInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_showColorbar.
function checkbox_showColorbar_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showColorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showColorbar

draw(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_displayMeshVisible.
function pushbutton_displayMeshVisible_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_displayMeshVisible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data

wait_h = waitbar(0,'Please wait... calculating the mesh');
if isfield(Data,'segangio')
    [Sz,Sx,Sy] = size(Data.segangio);
    Zstartframe = str2double(get(handles.edit_Zstartframe,'String'));
    Zstartframe = min(max(Zstartframe,1),Sz);
    ZMIP = str2double(get(handles.edit_ZMIP,'String'));
    Zendframe = min(max(Zstartframe+ZMIP-1,1),Sz);
    if isfield(Data,'ZoomXrange')
        Xstartframe = Data.ZoomXrange(1);
        Xendframe = Data.ZoomXrange(2);
    else
        Xstartframe = 1;
        Xendframe = Sx;
    end
     if isfield(Data,'ZoomYrange')
        Ystartframe = Data.ZoomYrange(1);
        Yendframe = Data.ZoomYrange(2);
    else
        Ystartframe = 1;
        Yendframe = Sy;
    end
        Mask = permute( Data.segangio(Zstartframe:Zendframe,Ystartframe:Yendframe,Xstartframe:Xendframe), [2 3 1]);
        fv = isosurface(Mask);
        fv2 = reducepatch(fv,200000); 
    
    
%     fv2 = fv;
    f = fv2.faces;
    v = fv2.vertices;
    figure(2);
    clf

    if eventdata~=1
        h=trisurf(f,v(:,1),v(:,2),v(:,3),'facecolor','red','edgecolor','none');
        daspect([1,1,1])
        view(3); axis tight
        camlight
        lighting gouraud
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
    end
    
    offset = [Xstartframe,Ystartframe,Zstartframe];
    save('mesh.mat','Mask','f','v','offset');
end
waitbar(1);
close(wait_h);


% --- Executes on button press in pushbuttonLoadGraph.
function pushbuttonLoadGraph_Callback(hObject, eventdata, handles)

global Data

[filename,pathname] = uigetfile({'*.mat'},'Please select the Graph Data to Load');
if filename==0
    return
end

load([pathname filename]);
if ~exist('seg')
    msgbox('Selected File does not have Graph information');
    return
end

%Data.Graph.seg = seg;
if isfield(Data,'Graph')
    nNodes = size(Data.Graph.nodes,1);
    nNewNodes = size(nodes,1);
    Data.Graph.nodes(nNodes+[1:nNewNodes],1:3) = nodes;

    nEdges = size(Data.Graph.edges,1);
    nNewEdges = size(edges,1);
    Data.Graph.edges(nEdges+[1:nNewEdges],1:2) = edges + nNodes;
else
    Data.Graph.nodes = nodes;
    Data.Graph.edges = edges;
end

set(handles.checkboxDisplayGraph,'enable','on')
draw(hObject, eventdata, handles)


% --- Executes on button press in checkboxDisplayGraph.
function checkboxDisplayGraph_Callback(hObject, eventdata, handles)
draw(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonGraphMesh.
function pushbuttonGraphMesh_Callback(hObject, eventdata, handles)

load mesh.mat

graphTubularMesh( f, v, Mask, offset );


% --- Executes on button press in pushbuttonCenterXYZ.
function pushbuttonCenterXYZ_Callback(hObject, eventdata, handles)
global Data

nodesNew = centerNodesXYZ( Data.Graph.nodes, Data.Graph.edges, permute(Data.segangio,[2 3 1]) ); % permute to x,y,z

Data.Graph.nodes = nodesNew;

draw(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonCenterXY.
function pushbuttonCenterXY_Callback(hObject, eventdata, handles)
global Data

nodesNew = centerNodes( Data.Graph.nodes, Data.Graph.edges, permute(Data.segangio,[2 3 1]) ); % permute to x,y,z

Data.Graph.nodes = nodesNew;

draw(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonRegraphNodes.
function pushbuttonRegraphNodes_Callback(hObject, eventdata, handles)
global Data

[nodes, edges] = regraphNodes( Data.Graph.nodes, Data.Graph.edges);
[nodes, edges] = fillNodes( nodes, edges);

Data.Graph.nodes = nodes;
Data.Graph.edges = edges;

draw(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonStraightenNodes.
function pushbuttonStraightenNodes_Callback(hObject, eventdata, handles)
global Data

nodesNew = straightenNodes( Data.Graph.nodes, Data.Graph.edges, permute(Data.segangio,[2 3 1]), 0.5 ); % permute to x,y,z

Data.Graph.nodes = nodesNew;

draw(hObject, eventdata, handles)




% --------------------------------------------------------------------
function Filter_expTransformation_Callback(hObject, eventdata, handles)
% hObject    handle to Filter_expTransformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data

prompt = {'Enter transformation constant :'};
defaultans = {'4'};
Tc = str2double(inputdlg(prompt,'Transformation Constant',1,defaultans ));

if isfield(Data,'angioF')
    I = double(Data.angioF);
else
    I = double(Data.angio);
end

h = waitbar(0,'Please wait... applying Exponential Transform');

Data.angioF = 1-exp(-Tc*I/max(I(:)));

if isfield(Data,'procSteps')
    Data.procSteps(end+1,:) =  {{'Exponential Tranformation'},{'Transformation constant'},{Tc}};
else
    Data.procSteps =  {{'Exponential Tranformation'},{'Transformation constant'},{Tc}};
end
waitbar(1);
close(h);
draw(hObject, eventdata, handles);




% --- Executes on button press in pushbuttonMeshGraphTiles.
function pushbuttonMeshGraphTiles_Callback(hObject, eventdata, handles)
global Data

[Sz,Sx,Sy] = size(Data.segangio);

xRange = str2num(get(handles.editGraphXrange,'string'));
yRange = str2num(get(handles.editGraphYrange,'string'));
zRange = str2num(get(handles.editGraphZrange,'string'));

set(handles.edit_ZMIP,'String',num2str(zRange(4)));
iiZ=0;
for iZ=zRange(1):zRange(2):zRange(3)
    iiZ = iiZ + 1;
    set(handles.edit_Zstartframe,'String',num2str(iZ))
    draw(hObject, eventdata, handles);
    pushbutton_displayMeshVisible_Callback(hObject, 0, handles);
    drawnow
    
    load mesh.mat
    
    filenm = sprintf('graph%02d.mat',iiZ);
    graphTubularMesh( f, v, Mask, offset, filenm );

    load(filenm);
    if isfield(Data,'Graph')
        nNodes = size(Data.Graph.nodes,1);
        nNewNodes = size(nodes,1);
        Data.Graph.nodes(nNodes+[1:nNewNodes],1:3) = nodes;
        
        nEdges = size(Data.Graph.edges,1);
        nNewEdges = size(edges,1);
        Data.Graph.edges(nEdges+[1:nNewEdges],1:2) = edges + nNodes;
    else
        Data.Graph.nodes = nodes;
        Data.Graph.edges = edges;
    end
        
end

set(handles.checkboxDisplayGraph,'enable','on')
draw(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonGraphClear.
function pushbuttonGraphClear_Callback(hObject, eventdata, handles)
global Data
if exist('Data')
    if isfield(Data,'Graph')
        Data = rmfield(Data,'Graph');
    end
end
draw(hObject, eventdata, handles)



function editGraphXrange_Callback(hObject, eventdata, handles)


function editGraphYrange_Callback(hObject, eventdata, handles)


function editGraphZrange_Callback(hObject, eventdata, handles)


