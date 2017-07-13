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

% Last Modified by GUIDE v2.5 13-Jul-2017 11:33:35

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

% UIWAIT makes vesSegment wait for user response (see UIRESUME)
% uiwait(handles.figure1);


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

[filename,pathname] = uigetfile({'*.mat;*.tiff;*.tif'},'Please select the Angiogram Data');
[pathstr,name,ext] = fileparts(filename);
if strcmp(ext,'.mat')
    temp = load([pathname filename]);
    fn = fieldnames(temp);
    Data.angio = temp.(fn{1});
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
end
[z,x,y] = size(Data.angio);
set(handles.edit_Zstartframe,'String',num2str(1));
set(handles.edit_ZMIP,'String',num2str(z));
set(handles.edit_Xstartframe,'String',num2str(1));
set(handles.edit_XMIP,'String',num2str(x));
set(handles.edit_Ystartframe,'String',num2str(1));
set(handles.edit_YMIP,'String',num2str(y));

draw(hObject, eventdata, handles);


function draw(hObject, eventdata, handles)

global Data

if isfield(Data,'angioF') && (get(handles.checkbox_showFiltData,'Value') == 1)
    I = Data.angioF;
elseif isfield(Data,'segangio') && (get(handles.checkbox_showFiltData,'Value') == 0) && (get(handles.checkbox_showRawData,'Value') == 0) && (get(handles.checkbox_showSeg,'Value') == 1)
    I = Data.segangio;
else
    I = Data.angio;
end

[Sz,Sx,Sy] = size(I);

Zstartframe = str2double(get(handles.edit_Zstartframe,'String'));
Zstartframe = min(max(Zstartframe,1),Sz);
ZMIP = str2double(get(handles.edit_ZMIP,'String'));
Zendframe = min(max(Zstartframe+ZMIP-1,1),Sz);
Xstartframe = str2double(get(handles.edit_Xstartframe,'String'));
Xstartframe = min(max(Xstartframe,1),Sx);
XMIP = str2double(get(handles.edit_XMIP,'String'));
Xendframe = min(max(Xstartframe+XMIP-1,1),Sx);
Ystartframe = str2double(get(handles.edit_Ystartframe,'String'));
Ystartframe = min(max(Ystartframe,1),Sy);
YMIP = str2double(get(handles.edit_YMIP,'String'));
Yendframe = min(max(Ystartframe+YMIP-1,1),Sy);

Zimg = squeeze(max(I(Zstartframe:Zendframe,:,:),[],1));
Ximg = squeeze(max(I(:,Xstartframe:Xendframe,:),[],2));
Yimg = squeeze(max(I(:,:,Ystartframe:Yendframe),[],3));
ZMIPimg = squeeze(max(I,[],1));

if isfield(Data,'segangio') && (get(handles.checkbox_showSeg,'Value') == 1)
    ZimgS = squeeze(max(Data.segangio(Zstartframe:Zendframe,:,:),[],1));
    XimgS = squeeze(max(Data.segangio(:,Xstartframe:Xendframe,:),[],2));
    YimgS = squeeze(max(Data.segangio(:,:,Ystartframe:Yendframe),[],3));
    imgS = squeeze(max(Data.segangio,[],1));
end

axes(handles.axes1)
colormap('gray')
imagesc(Zimg)
if (get(handles.checkbox_showFiltData,'Value') == 1) || (get(handles.checkbox_showRawData,'Value') == 1)
    if isfield(Data,'angioT') && (get(handles.checkbox_showSeg,'Value') == 1)
        hold on
        img = double(ZimgS);
        green = cat(3, zeros(size(img)),ones(size(img)), zeros(size(img)));
        hold on
        h = imshow(green);
        hold off
        set(h, 'AlphaData', img*0.25)
    end
end
axis image
axis on

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



function edit_Xstartframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Xstartframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Xstartframe as text
%        str2double(get(hObject,'String')) returns contents of edit_Xstartframe as a double
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_Xstartframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Xstartframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ystartframe_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ystartframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ystartframe as text
%        str2double(get(hObject,'String')) returns contents of edit_Ystartframe as a double
draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_Ystartframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ystartframe (see GCBO)
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



function edit_XMIP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_XMIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_XMIP as text
%        str2double(get(hObject,'String')) returns contents of edit_XMIP as a double


draw(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_XMIP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_XMIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_YMIP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_YMIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_YMIP as text
%        str2double(get(hObject,'String')) returns contents of edit_YMIP as a double

draw(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit_YMIP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_YMIP (see GCBO)
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
ii = str2double(get(handles.edit_Xstartframe,'String'));
ii = ii-1;
set(handles.edit_Xstartframe,'String',num2str(ii));
draw(hObject, eventdata, handles);

% --- Executes on button press in pushbutton_Ymoveleft.
function pushbutton_Ymoveleft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Ymoveleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ii = str2double(get(handles.edit_Zstartframe,'String'));
ii = ii-1;
set(handles.edit_Ystartframe,'String',num2str(ii));
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
ii = str2double(get(handles.edit_Xstartframe,'String'));
ii = ii+1;
set(handles.edit_Xstartframe,'String',num2str(ii));
draw(hObject, eventdata, handles);

% --- Executes on button press in pushbutton_Ymoveright.
function pushbutton_Ymoveright_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Ymoveright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ii = str2double(get(handles.edit_Ystartframe,'String'));
ii = ii+1;
set(handles.edit_Ystartframe,'String',num2str(ii));
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
waitbar(1);
close(h);

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
waitbar(1);
close(h);


% --------------------------------------------------------------------
function Filters_TubenessFilter_Callback(hObject, eventdata, handles)
% hObject    handle to Filters_TubenessFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
if isfield(Data,'angioF')
    I = Data.angioF;
else
    I = Data.angio;
end
I = double(I);
Z = size(I,3);
beta = 100;
c = 500;
[k,l,m] = size(I);
[nz,nx,ny] = size(I);
L   = zeros(k,l,m);
Vs  = zeros(k,l,m);
alpha = 0.25;
gamma12 = 0.5;
gamma23 = 0.5;
T = zeros(k,l,m);
sigma = [3];
h = waitbar(0,'Please wait... applying Tubeness filter');
for i = 1:length(sigma)
    
    waitbar(i-1/length(sigma));
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigma(i));
    
%     Normalizing the Hessian Matrix
%     Dxx = i^2*Dxx; Dyy = i^2*Dyy;  Dzz = i^2*Dzz; Dxy = i^2*Dxy;  Dxz = i^2*Dxz; Dyz = i^2*Dyz;
    
    
    
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
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(T,sigma);
    
    [Lambda1,Lambda2,Lambda3,V1,V2,V3] = eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    
    SortL = sort([Lambda1(:)'; Lambda2(:)'; Lambda3(:)'],'ascend');
    Lambda1 = reshape(SortL(1,:),size(Lambda1));
    Lambda2 = reshape(SortL(2,:),size(Lambda2));
    Lambda3 = reshape(SortL(3,:),size(Lambda3));
    
    E = -sigma^2.*Lambda2;
    E(E<0) = 0;
end
%     T = L;
T = E;

T = (T-min(T(:)))/(max(T(:))-min(T(:)));
Data.angioT = T;
close(h);

% --------------------------------------------------------------------
function Filters_thresholding_Callback(hObject, eventdata, handles)
% hObject    handle to Filters_thresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data

if isfield(Data,'angioT')
    prompt = {'Please enter threshold value for segmentation'};
    defaultans = {'0.075'};
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
end


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

% --------------------------------------------------------------------
function File_savedata_Callback(hObject, eventdata, handles)
% hObject    handle to File_savedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Data
if isfield(Data,'angioF')|| isfield(Data,'angioT') || isfield(Data,'segangio')
    if isfield(Data,'angioF')
        Output.angioF = Data.angioF;
    end
    if isfield(Data,'angioT')
        Output.angioT = Data.angioT;
    end
    if isfield(Data,'segangio')
        Output.segangio = Data.segangio;
    end
    [FileName,PathName] = uiputfile('*.mat');
    save([PathName FileName],'Output');
end
    

