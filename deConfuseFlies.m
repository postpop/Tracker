function varargout = deConfuseFlies(varargin)
% DECONFUSEFLIES MATLAB code for deConfuseFlies.fig
%      DECONFUSEFLIES, by itself, creates a new DECONFUSEFLIES or raises the existing
%      singleton*.
%
%      H = DECONFUSEFLIES returns the handle to a new DECONFUSEFLIES or the handle to
%      the existing singleton*.
%
%      DECONFUSEFLIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECONFUSEFLIES.M with the given input arguments.
%
%      DECONFUSEFLIES('Property','Value',...) creates a new DECONFUSEFLIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deConfuseFlies_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deConfuseFlies_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% TODO:
% - allow manual entry of precise frame number or frame range
% - DROP: implement with VLC - NO, current version is fast enough
% - DONE: save results
% - DONE: allow shifting frame range
% - DONE: rotate frames to fix orientations
% - DONE: allow selecting frame range in lower plot
% - DONE: fix zoom out function
% - DROP: reuse previously loaded frames - NO, using pre-boxed vids is fast
% enough
% - DONE if first/last are adjacent: allow swap with first/last frame as endpoint
% - DONE: use preprocessed flyframe movie if available

% Edit the above text to modify the response to help deConfuseFlies

% Last Modified by GUIDE v2.5 03-Apr-2015 17:26:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
   'gui_Singleton',  gui_Singleton, ...
   'gui_OpeningFcn', @deConfuseFlies_OpeningFcn, ...
   'gui_OutputFcn',  @deConfuseFlies_OutputFcn, ...
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


% --- Executes just before deConfuseFlies is made visible.
function deConfuseFlies_OpeningFcn(hObject, eventdata, handles, varargin)
global p
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to deConfuseFlies (see VARARGIN)

if isempty(varargin)
   p.dirName = cd;
else
   p.dirName = fileparts(varargin{1});
end

if ~isfield(p,'init')
   p.init = 1;
else
   p.init = 0;
end
% populate file list and initialize all parameters
if p.init
   disp(['parsing ' p.dirName ' for tracking results'])
   fileList = rdir(fullfile(p.dirName, ['**' filesep '*_res*.mat']));
   p.fileNames = {fileList.name}';
   set(handles.popupmenu1, 'String', p.fileNames);
   axis(handles.axes1,'off','square')
   p.init = 0;
end
% load current fileName
p.currentFileName = p.fileNames{get(handles.popupmenu1,'Value')};
disp(['   loading results file ' p.currentFileName])
load(p.currentFileName,'vDat');
try
   p.vDat = vDat;
catch
   p.vDat = p;
end
tmp = load([p.currentFileName],'p','fp');
p.fp = tmp.fp;
try
   p.fp.tracks = tmp.p.tracks;
   p.fp.orientation = tmp.p.orientation;
catch
   p.fp.tracks = tmp.fp.tracks;  
   p.fp.orientation = tmp.fp.orientation;
end

% p.fp.orientation = tmp.p.orientation;
p.fp = fix.fixOrientations(p.fp);
oriX = p.fp.orientation(:,:,1);
oriY = p.fp.orientation(:,:,2);
p.flyAngle = cart2pol(oriX, oriY);% in RAD
p.flyAngle = p.flyAngle/2/pi*360;% in DEG
% minimize random flips and orient forwards based on speed


p.flyIdx(1:size(p.fp.tracks,1),[1 2]) = 1;
p.flyIdx(1:size(p.fp.tracks,1),2) = 2;
p.currentFly = 1; % select fly 1 per default

[currentFileDir, currentFileNam, currentFileExt] = fileparts(p.currentFileName);

if exist([currentFileDir '/' currentFileNam(1:11) '_fly.mj2'],'file');
   p.flyBoxVidExists = true;
   p.vidFileName = [currentFileDir '/' currentFileNam(1:11) '_fly.mj2'];
else
   p.flyBoxVidExists = false;
   if exist([currentFileDir '/' currentFileNam(1:11) '.mp4'], 'file')
      p.vidFileName = [currentFileDir '/' currentFileNam(1:11) '.mp4'];
   end
   if exist([currentFileDir '/' currentFileNam(1:11) '.avi'], 'file')
      p.vidFileName = [currentFileDir '/' currentFileNam(1:11) '.avi'];
   end
end

disp(['   loading ' p.vidFileName])
% try
[fileDir, fileNam, fileExt] = fileparts(p.vidFileName);
%    owd = cd;
%    cd(fileDir)
%    vr = VideoReader2([fileNam fileExt]);
%    cd(owd);
%    p.fp.vr = VideoReader(p.vidFileName);
%    %p.fp.vr.vr.getFrameAtNum(1);
%    p.fp.vr.read(1);
% catch
%    disp('falling back to built-in VideoReader.')

% if ispc
% % else
% %    p.fp.vr = VideoReaderFFMPEG(p.vidFileName);
% end

p.fp.vr = VideoReader(p.vidFileName);
p.fp.vr.read(1);
% end
p.imageChannels = size(p.fp.vr.read(1),3);
try
   p.boxW = tmp.p.flyBoxWidthIdx;%-60:61;
   p.boxH = tmp.p.flyBoxHeightIdx;%-60:61;
catch
   p.boxW = -40:41;
   p.boxH = -40:41;
end
p.NumFramesToRead = 2^nextpow2(128);% always make sure that this is square

% selected frame in differen coordinates:
p.selFrameSub = zeros(2,1);   % ?
p.selFrameIdx = zeros(2,2);   % ?
p.selFrame = zeros(1,2);      % ?
p.framesToRead = round(linspace(p.fp.initFrame,p.fp.vr.NumberOfFrames,p.NumFramesToRead));
p.frames = fix.extractFlyFrames(p);
plotFrames(handles);


set(handles.edit1,'String',int2str(p.fp.initFrame));
p.skipSize = str2double(get(handles.edit3, 'String'));

% Choose default command line output for deConfuseFlies
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes deConfuseFlies wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = deConfuseFlies_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
global p
% reset gui, load selected fle
disp('   selecting new file.')
deConfuseFlies_OpeningFcn(hObject, [], handles, p.dirName)


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end

% ---- CALLED WHEN CLICKED ON A DISPLAYED FRAME
function ImageClickCallback ( objectHandle , eventData, handles )
global p
axesHandle  = get(objectHandle,'Parent');
coordinates = get(axesHandle,'CurrentPoint');
coordinates = coordinates(1,1:2);
% convert mouse click coords to frame idx (linear and x,y)
p.selFrameIdx(1,:) = p.selFrameIdx(2,:);% a poor man's FIFO queue
p.selFrameIdx(2,:) = ceil(coordinates./length(p.boxW));% convert to ??
p.selFrameSub(1) = p.selFrameSub(2);% a poor man's FIFO queue
p.selFrameSub(2) = (p.selFrameIdx(2,2)-1)*ceil(sqrt(p.NumFramesToRead))+p.selFrameIdx(2,1);
p.selFrame(1) = p.selFrame(2);
if p.selFrameSub(2)>0
   p.selFrame(2) = p.framesToRead(p.selFrameSub(2));% a poor man's FIFO queue
else
   p.selFrame(2) = 0;
end
disp(p.selFrame)
% plot frames
p.imh = plotFrames(handles);


% ---- CALLED WHEN CLICKED ON A TRACKS
function ImageClickCallbackTrack ( objectHandle , eventData, handles )
global p
coordinates = get(objectHandle,'CurrentPoint');
% swap frames and convert mouse click coords to frame idx (linear and x,y)
p.selFrame(1) = p.selFrame(2);
p.selFrame(2) = (coordinates(1));
p.selFrameIdx(1) = p.selFrameIdx(2);
p.selFrameIdx(2) = 0;
p.selFrameSub(1) = p.selFrameSub(2);
p.selFrameSub(2) = 0;
disp(p.selFrame)
% plot frames
plotFrames(handles);



function imh = plotFrames(handles)
global p
axes(handles.axes1)
%imh = montage(p.frames(:,(p.currentFly-1)*p.fp.nFlies + (1:length(p.boxH)),:,:));
imh = montage(p.frames);% plot frames
set(imh,'ButtonDownFcn',{@ImageClickCallback handles});

% plot boxes around selected frames
if p.selFrameSub(1)>0
   rh = rectangle('Position',[(p.selFrameIdx(1,1)-1)*length(p.boxW),(p.selFrameIdx(1,2)-1)*length(p.boxH), length(p.boxW), length(p.boxH)]);
   set(rh,'LineWidth',2,'EdgeColor','r');
end
if p.selFrameSub(2)>0
   rh = rectangle('Position',[(p.selFrameIdx(2,1)-1)*length(p.boxW),(p.selFrameIdx(2,2)-1)*length(p.boxH), length(p.boxW), length(p.boxH)]);
   set(rh,'LineWidth',2,'EdgeColor','b');
end
% plot displayed frame ticks and flyidx below
axes(handles.axes2)
cla
plot([p.framesToRead;p.framesToRead], repmat(1:p.fp.nFlies,length(p.framesToRead),1)','k')
set(gca,'XLim',[p.fp.initFrame p.fp.vr.NumberOfFrames]);%,'YTick',1:p.fp.nFlies,'box','off')
hold on
if min(p.selFrame)>0
   plot(min(p.selFrame)*[1;1], [1; p.fp.nFlies], 'LineWidth',2,'Color','r')
end
if max(p.selFrame)>0
   plot(max(p.selFrame)*[1;1], [1; p.fp.nFlies], 'LineWidth',2,'Color','b')
end
plot(p.fp.initFrame + repmat(1:length(p.flyIdx),2,1)',  p.flyIdx)
hold off
set(gca,'ButtonDownFcn',{@ImageClickCallbackTrack handles});


% --- ZOOM IN
function pushbutton1_Callback(hObject, eventdata, handles)
global p
% grab p.NumFramesToRead frames between the selected points
%p.framesToRead = round(linspace(p.framesToRead(min(p.selFrameSub)), p.framesToRead(max(p.selFrameSub)), p.NumFramesToRead));
p.framesToRead = round(linspace(min(p.selFrame), max(p.selFrame), p.NumFramesToRead));
p.frames = fix.extractFlyFrames(p);
% set selected frames to beginning and start of this new seq
p.selFrameSub = [1 p.NumFramesToRead];% linear frame index
[p.selFrameIdx(:,1), p.selFrameIdx(:,2)] = ind2sub([4 , 4],p.selFrameSub);% [x, y] frame index
% plot frames
p.imh = plotFrames(handles);


% --- ZOOM RESET
function pushbutton2_Callback(hObject, eventdata, handles)
global p
p.framesToRead = round(linspace(p.fp.initFrame, p.fp.vr.NumberOfFrames, p.NumFramesToRead));
p.frames = fix.extractFlyFrames(p);
p.selFrameSub = [1 p.NumFramesToRead];
p.selFrame = [p.fp.initFrame p.fp.vr.NumberOfFrames];
plotFrames(handles);



% --- SWAP FLIES
function pushbutton3_Callback(hObject, eventdata, handles)
global p
% get swap points
swapPoints = p.framesToRead(sort(p.selFrameSub)) - p.fp.initFrame + 1;% convert frameidx to array indices by correcting for initFrame
% swap
if diff(swapPoints)==1
   p.flyIdx(swapPoints(1)+1:end,:) = fliplr(p.flyIdx(swapPoints(1)+1:end,:));
else
   p.flyIdx(swapPoints(1):swapPoints(2),:) = fliplr(p.flyIdx(swapPoints(1):swapPoints(2),:));
end
%plot new traces
p.frames = fix.extractFlyFrames(p);
p.imh = plotFrames(handles);


% --- ZOOM OUT
function pushbutton4_Callback(hObject, eventdata, handles)
global p
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dFrame = round((p.framesToRead(end) - p.framesToRead(1))*p.NumFramesToRead);
dFrame = limit(dFrame,1,floor((p.fp.vr.NumberOfFrames-p.fp.initFrame)/p.NumFramesToRead));% make sure step size is such that p.NumFramesToRead steps fit in number of fraes
p.framesToRead = round(mean(p.framesToRead)) + (-floor(p.NumFramesToRead/2):ceil(p.NumFramesToRead/2))*dFrame;

% if p.framesToRead(end)>p.fp.vr.NumberOfFrames% hit upper end - go back from end
%    p.framesToRead = p.fp.vr.NumberOfFrames - (p.NumFramesToRead:-1:1)*dFrame;
% elseif p.framesToRead(end)<p.fp.initFrame% hit upper end - go back from end
%    p.framesToRead = p.fp.initFrame + (1:p.NumFramesToRead)*dFrame;
% end
p.framesToRead = limit(p.framesToRead, p.fp.initFrame,p.fp.vr.NumberOfFrames);
p.frames = fix.extractFlyFrames(p);
p.selFrameSub = [floor(p.NumFramesToRead/2) ceil(p.NumFramesToRead/2)];
p.selFrame = p.framesToRead(p.selFrameSub);
plotFrames(handles);

function edit1_Callback(hObject, eventdata, handles)  % START FRAME
global p
p.selFrame(1) = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)        % SKIP FRAMES
global p
p.skipSize = str2double(get(handles.edit3,'String'));

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)       % SKIP FRAMES
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
%
% function edit2_Callback(hObject, eventdata, handles)
% p.selFrame(2) = str2double(get(hObject,'String'));
%
% % --- Executes during object creation, after setting all properties.
% function edit2_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
global p
if get(handles.radiobutton1, 'Value')
   p.currentFly = 1;
else
   p.currentFly = 2;
end
disp(['   selected fly ' int2str(p.currentFly) '.'])


% --- Executes on button press in pushbutton5. SAVE BUTTON
function pushbutton5_Callback(hObject, eventdata, handles)
global p
saveFileName = [p.currentFileName(1:end-4) '_fixed.mat'];
disp(['   saving to ' saveFileName])
fp = p.fp;
save(saveFileName, 'p','fp')
disp(['   done.'])


% --- Executes on button press in pushbutton6. FWD
function pushbutton6_Callback(hObject, eventdata, handles)
global p
p.framesToRead = round(max(p.framesToRead)+(0:p.skipSize:p.skipSize*p.NumFramesToRead));
p.framesToRead = limit(p.framesToRead, p.fp.initFrame,p.fp.vr.NumberOfFrames);
p.frames = fix.extractFlyFrames(p);
p.selFrameSub = [1 p.NumFramesToRead];
p.selFrame = [p.fp.initFrame p.fp.vr.NumberOfFrames];
plotFrames(handles);



% --- Executes on button press in pushbutton7. REV
function pushbutton7_Callback(hObject, eventdata, handles)
global p
p.framesToRead = round(min(p.framesToRead)-p.skipSize*p.NumFramesToRead + (0:p.skipSize:p.skipSize*p.NumFramesToRead));
p.framesToRead = limit(p.framesToRead, p.fp.initFrame,p.fp.vr.NumberOfFrames);
p.frames = fix.extractFlyFrames(p);
p.selFrameSub = [1 p.NumFramesToRead];
p.selFrame = [p.fp.initFrame p.fp.vr.NumberOfFrames];
plotFrames(handles);


% --- Executes on button press in pushbutton8.           INIT FRAME VIEW
function pushbutton8_Callback(hObject, eventdata, handles)
global p
p.framesToRead = round(p.fp.initFrame:p.skipSize:p.skipSize*p.NumFramesToRead);
p.framesToRead = limit(p.framesToRead, p.fp.initFrame,p.fp.vr.NumberOfFrames);
p.frames = fix.extractFlyFrames(p);
p.selFrameSub = [1 p.NumFramesToRead];
p.selFrame = [p.fp.initFrame p.fp.vr.NumberOfFrames];
plotFrames(handles);


%  ROTATE
function pushbutton9_Callback(hObject, eventdata, handles)
global p

% get rotation points
swapPoints = p.framesToRead(sort(p.selFrameSub));
% rotate CW by adding 180 to all angles beyond rotation point
if diff(swapPoints)<=1
   endIdx = size(p.flyAngle,1);
else
   endIdx = swapPoints(2);
end
thisFlyIdx = p.flyIdx(swapPoints(1)+1:endIdx,p.currentFly);
linearInd = sub2ind(size(p.flyAngle), swapPoints(1)+1:endIdx, thisFlyIdx');
p.flyAngle(linearInd) = p.flyAngle(linearInd) + 180;
%plot new traces
p.frames = fix.extractFlyFrames(p);
p.imh = plotFrames(handles);
