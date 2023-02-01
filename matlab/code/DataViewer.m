function varargout = DataViewer(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DataViewer_OpeningFcn, ...
    'gui_OutputFcn',  @DataViewer_OutputFcn, ...
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


%UPDATE TIME SERIES
function handles2give = UpdateTimeSeries(handles2give)
%% General parameters
LineWidth = 1.0;
FontSize = 12;

data = handles2give.data;
SelectedTrial = strcmp(data.Mouse_Name,handles2give.MouseName) & ...
    data.Cell_Counter == handles2give.CellNumber & ...
    data.Sweep_Counter == handles2give.TrialNumber;

%% Update Information Display
set(handles2give.InformationDisplay_MouseName_Tag, 'String', handles2give.MouseName);
set(handles2give.InformationDisplay_CellNumber_Tag, 'String', handles2give.CellNumber);
set(handles2give.InformationDisplay_TrialNumber_Tag, 'String', handles2give.TrialNumber);
ExtractMouseGenotype = data.Mouse_Genotype(SelectedTrial);
set(handles2give.InformationDisplay_MouseGenotype_Tag, 'String', ExtractMouseGenotype{1,1});
ExtractTrialStartTime = datestr(datetime(data.Sweep_StartTime(SelectedTrial,:)));
set(handles2give.InformationDisplay_TrialStartTime_Tag, 'String', ExtractTrialStartTime);
ExtracttdTomato = data.Cell_tdTomatoExpressing(SelectedTrial);
if strcmp(ExtracttdTomato{1,1}, 'True')
    set(handles2give.InformationDisplay_CelltdTomato_Tag, 'String', 'tdTomato+');
else
    set(handles2give.InformationDisplay_CelltdTomato_Tag, 'String', 'tdTomato-');
end
ExtractCellDepth = data.Cell_Depth(SelectedTrial);
set(handles2give.InformationDisplay_CellDepth_Tag, 'String', num2str(ExtractCellDepth));
ExtractTrialType = data.Sweep_Type(SelectedTrial);
set(handles2give.InformationDisplay_TrialType_Tag, 'String', ExtractTrialType{1,1});
ExtractCellType = data.Cell_Type(SelectedTrial);
set(handles2give.InformationDisplay_CellType_Tag, 'String', ExtractCellType{1,1});


%% Data to plot

% Membrane Potential
ExtractMembranePotential = data.Sweep_MembranePotential(SelectedTrial);
handles2give.MembranePotential = ExtractMembranePotential{1,1};
handles2give.MembranePotentialTimeVector = 1:(length(handles2give.MembranePotential));
ExtractMembranePotential_SamplingRate = data.Sweep_MembranePotential_SamplingRate(SelectedTrial);
MembranePotential_SamplingRate = ExtractMembranePotential_SamplingRate;
handles2give.MembranePotentialTimeVector = handles2give.MembranePotentialTimeVector * (1000 / MembranePotential_SamplingRate);

% Whisker Angle
ExtractWhiskerAngle = data.Sweep_WhiskerAngle(SelectedTrial);
handles2give.WhiskerAngle = ExtractWhiskerAngle{1,1};
handles2give.WhiskerAngleTimeVector = 1:(length(handles2give.WhiskerAngle));
WhiskerAngle_SamplingRate = data.Sweep_WhiskerAngle_SamplingRate(SelectedTrial);
handles2give.WhiskerAngleTimeVector = handles2give.WhiskerAngleTimeVector * (1000 / WhiskerAngle_SamplingRate);

% Contact Time
ExtractContactTime = data.Sweep_ActiveContactTimes(SelectedTrial);
handles2give.ContactTime = ExtractContactTime{1,1};
% set the same time vector as the membrane potential
handles2give.ContactTimeTimeVector = handles2give.MembranePotentialTimeVector; 
if isnan(handles2give.ContactTime)    % check if there is any contact time recorded
    handles2give.ContactTime = zeros(length(handles2give.ContactTimeTimeVector),1); 
else
    % ContactTime has two columns (beginning of the contact and the end)
    % make a vector of ContactTime and convert in ms
    handles2give.ContactTime = reshape(handles2give.ContactTime,(2*length(handles2give.ContactTime(:,1))),1)*1000;
    % sort to have beginning and ending of contact in successive manner
    handles2give.ContactTime = sort(handles2give.ContactTime);
    % round to have ms precision
    handles2give.ContactTime = round(handles2give.ContactTime);
    % initialize vector containing y-values at 0
    temp_ContactTime = zeros(length(handles2give.ContactTimeTimeVector),1); 
    for i=1:2:length(handles2give.ContactTime)-1
        position1 = find(handles2give.ContactTimeTimeVector == handles2give.ContactTime(i), 1, 'first');
        position2 = find(handles2give.ContactTimeTimeVector == handles2give.ContactTime(i+1), 1, 'last');
        for j=position1:position2
            temp_ContactTime(j) = 1;
        end
    end
     handles2give.ContactTime = temp_ContactTime;
end

% Whisking Time
ExtractWhiskingTime = data.Sweep_WhiskingTimes(SelectedTrial);
handles2give.WhiskingTime = ExtractWhiskingTime{1,1};
% set the same time vector as the membrane potential
handles2give.WhiskingTimeTimeVector = handles2give.MembranePotentialTimeVector; 
if isempty(handles2give.WhiskingTime)    % check if there is any whisking time recorded
    handles2give.WhiskingTime = zeros(length(handles2give.WhiskingTimeTimeVector),1); 
else
    % WhiskingTime has two columns (beginning of the whisking and the end)
    % make a vector of WhiskingTime and convert in ms
    handles2give.WhiskingTime = reshape(handles2give.WhiskingTime,(2*length(handles2give.WhiskingTime(:,1))),1)*1000;
    % sort to have beginning and ending of whisking in successive manner
    handles2give.WhiskingTime = sort(handles2give.WhiskingTime);
    % round to have ms precision
    handles2give.WhiskingTime = round(handles2give.WhiskingTime);
    % initialize vector containing y-values at 0
    temp_WhiskingTime = zeros(length(handles2give.WhiskingTimeTimeVector),1); 
    for i=1:2:length(handles2give.WhiskingTime)-1
        position1 = find(handles2give.WhiskingTimeTimeVector == handles2give.WhiskingTime(i), 1, 'first');
        position2 = find(handles2give.WhiskingTimeTimeVector == handles2give.WhiskingTime(i+1), 1, 'last');
        for j=position1:position2
            temp_WhiskingTime(j) = 1;
        end
    end
     handles2give.WhiskingTime = temp_WhiskingTime;
end

% Quiet Time
ExtractQuietTime = data.Sweep_QuietTimes(SelectedTrial);
handles2give.QuietTime = ExtractQuietTime{1,1};
% set the same time vector as the membrane potential
handles2give.QuietTimeTimeVector = handles2give.MembranePotentialTimeVector; 
if isempty(handles2give.QuietTime)    % check if there is any quiet time recorded
    handles2give.QuietTime = zeros(length(handles2give.QuietTimeTimeVector),1); 
else
    % QuietTime has two columns (beginning of the contact and the end)
    % make a vector of QuietTime and convert in ms
    handles2give.QuietTime = reshape(handles2give.QuietTime,(2*length(handles2give.QuietTime(:,1))),1)*1000;
    % sort to have beginning and ending of quiet in successive manner
    handles2give.QuietTime = sort(handles2give.QuietTime);
    % round to have ms precision
    handles2give.QuietTime = round(handles2give.QuietTime);
    % initialize vector containing y-values at 0
    temp_QuietTime = zeros(length(handles2give.QuietTimeTimeVector),1); 
    for i=1:2:length(handles2give.QuietTime)-1
        position1 = find(handles2give.QuietTimeTimeVector == handles2give.QuietTime(i), 1, 'first');
        position2 = find(handles2give.QuietTimeTimeVector == handles2give.QuietTime(i+1), 1, 'last');
        for j=position1:position2
            temp_QuietTime(j) = 1;
        end
    end
     handles2give.QuietTime = temp_QuietTime;
end

%% Plot membrane potential
axes(handles2give.MembranePotentialAxes);
plot(handles2give.MembranePotentialTimeVector, handles2give.MembranePotential,'Color','k','LineWidth',LineWidth);
set(gca,'FontSize',FontSize);
xlabel(['Time (ms)'],'FontWeight','Bold','FontSize',FontSize);
ylabel('Membrane potential (V)','FontWeight','Bold','FontSize',FontSize);
axis([handles2give.TimeAxisMin handles2give.TimeAxisMax handles2give.MembranePotentialAxisMin handles2give.MembranePotentialAxisMax])
%set(gca,'XTickLabel',[])
box(gca,'off')

%% Plot right whisker angle
axes(handles2give.RightWhiskerAngleAxes);
plot(handles2give.WhiskerAngleTimeVector, handles2give.WhiskerAngle,'Color','g','LineWidth',LineWidth);
set(gca,'FontSize',FontSize);
xlabel(['Time (ms)'],'FontWeight','Bold','FontSize',FontSize);
ylabel('Right whisker (deg)','FontWeight','Bold','FontSize',FontSize);
axis([handles2give.TimeAxisMin handles2give.TimeAxisMax handles2give.WhiskerAxisMin handles2give.WhiskerAxisMax])
%set(gca,'XTickLabel',[])
box(gca,'off')

%% Plot whisking time
axes(handles2give.WhiskingTimeAxes);
plot(handles2give.WhiskingTimeTimeVector, handles2give.WhiskingTime,'Color','g','LineWidth',LineWidth);
set(gca,'FontSize',FontSize);
xlabel(['Time (ms)'],'FontWeight','Bold','FontSize',FontSize);
ylabel('Whisking','FontWeight','Bold','FontSize',FontSize);
axis([handles2give.TimeAxisMin handles2give.TimeAxisMax handles2give.WhiskingTimeAxisMin handles2give.WhiskingTimeAxisMax])
%set(gca,'XTickLabel',[])
box(gca,'off')

%% Plot quiet time
axes(handles2give.QuietTimeAxes);
plot(handles2give.QuietTimeTimeVector, handles2give.QuietTime,'Color','b','LineWidth',LineWidth);
set(gca,'FontSize',FontSize);
xlabel(['Time (ms)'],'FontWeight','Bold','FontSize',FontSize);
ylabel('Quiet','FontWeight','Bold','FontSize',FontSize);
axis([handles2give.TimeAxisMin handles2give.TimeAxisMax handles2give.QuietTimeAxisMin handles2give.QuietTimeAxisMax])
%set(gca,'XTickLabel',[])
box(gca,'off')

%% Plot contact time
axes(handles2give.ContactTimeAxes);
plot(handles2give.ContactTimeTimeVector, handles2give.ContactTime,'Color','b','LineWidth',LineWidth);
set(gca,'FontSize',FontSize);
xlabel(['Time (ms)'],'FontWeight','Bold','FontSize',FontSize);
ylabel('Contact','FontWeight','Bold','FontSize',FontSize);
axis([handles2give.TimeAxisMin handles2give.TimeAxisMax handles2give.ContactTimeAxisMin handles2give.ContactTimeAxisMax])
%set(gca,'XTickLabel',[])
box(gca,'off')



% --- Executes just before DataViewer is made visible.
function DataViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataViewer (see VARARGIN)

% Choose default command line output for DataViewer
handles.output = hObject;

%% Set axis scaling
handles.MembranePotentialAxisMin = -0.080; set(handles.MembranePotentialAxisMinTag,'String',num2str(handles.MembranePotentialAxisMin));
handles.MembranePotentialAxisMax = 0.040; set(handles.MembranePotentialAxisMaxTag,'String',num2str(handles.MembranePotentialAxisMax));

handles.WhiskerAxisMin = 175; set(handles.WhiskerAxisMinTag,'String',num2str(handles.WhiskerAxisMin));
handles.WhiskerAxisMax = 250; set(handles.WhiskerAxisMaxTag,'String',num2str(handles.WhiskerAxisMax));

handles.WhiskingTimeAxisMin = -0.2; set(handles.WhiskingTimeAxisMinTag,'String',num2str(handles.WhiskingTimeAxisMin));
handles.WhiskingTimeAxisMax = 1.2; set(handles.WhiskingTimeAxisMaxTag,'String',num2str(handles.WhiskingTimeAxisMax));

handles.QuietTimeAxisMin = -0.2; set(handles.QuietTimeAxisMinTag,'String',num2str(handles.QuietTimeAxisMin));
handles.QuietTimeAxisMax = 1.2; set(handles.QuietTimeAxisMaxTag,'String',num2str(handles.QuietTimeAxisMax));

handles.ContactTimeAxisMin = -0.2; set(handles.ContactTimeAxisMinTag,'String',num2str(handles.ContactTimeAxisMin));
handles.ContactTimeAxisMax = 1.2; set(handles.ContactTimeAxisMaxTag,'String',num2str(handles.ContactTimeAxisMax));

handles.TimeAxisMin = 0; set(handles.TimeAxisMinTag,'String',num2str(handles.TimeAxisMin));
handles.TimeAxisMax = 3000; set(handles.TimeAxisMaxTag,'String',num2str(handles.TimeAxisMax));

%% Load data
% cd ..
currentFolder = pwd;
if ispc
    MatlabFile = fullfile([currentFolder '\Data' '\*.mat']);
    MouseFileName = dir(MatlabFile);
    load([currentFolder '\Data\' MouseFileName.name]);
elseif ismac
    MatlabFile = fullfile([currentFolder '/Data' '/*.mat']);
    MouseFileName = dir(MatlabFile);
    load([currentFolder '/Data/' MouseFileName.name]);
end
handles.data = Data;

%% Set mouse name
handles.MouseNameList = unique(Data.Mouse_Name);
set(handles.MouseName_Popupmenu,'String',handles.MouseNameList);
set(handles.MouseName_Popupmenu,'Value',1);
handles.MouseName = handles.MouseNameList{get(handles.MouseName_Popupmenu,'Value')};

%% Set cell number
handles.CellNumberList = unique(handles.data.Cell_Counter(strcmp(handles.data.Mouse_Name,handles.MouseName)));
set(handles.CellNumber_Popupmenu,'String',handles.CellNumberList);
set(handles.CellNumber_Popupmenu,'Value',1);
handles.CellNumber = handles.CellNumberList(get(handles.CellNumber_Popupmenu,'Value'));

%% Set trial number
handles.TrialNumberList = unique(handles.data.Sweep_Counter(strcmp(handles.data.Mouse_Name,handles.MouseName) & handles.data.Cell_Counter == handles.CellNumber));
set(handles.TrialNumber_Popupmenu,'String',handles.TrialNumberList);
set(handles.TrialNumber_Popupmenu,'Value',1);
handles.TrialNumber = handles.TrialNumberList(get(handles.TrialNumber_Popupmenu,'Value'));

%% Set axes
handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DataViewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in MouseName_Popupmenu.
function MouseName_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to MouseName_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MouseName_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MouseName_Popupmenu

handles.MouseName = handles.MouseNameList{get(handles.MouseName_Popupmenu,'Value')};

handles.CellNumberList = unique(handles.data.Cell_Counter(strcmp(handles.data.Mouse_Name,handles.MouseName)));
set(handles.CellNumber_Popupmenu,'String',handles.CellNumberList);
handles.CellNumber = handles.CellNumberList(1);
set(handles.CellNumber_Popupmenu,'Value',1);

handles.TrialNumberList = unique(handles.data.Sweep_Counter(strcmp(handles.data.Mouse_Name,handles.MouseName) & handles.data.Cell_Counter == handles.CellNumber));
set(handles.TrialNumber_Popupmenu,'String',handles.TrialNumberList);
handles.TrialNumber = handles.TrialNumberList(1);
set(handles.TrialNumber_Popupmenu,'Value',1);

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MouseName_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MouseName_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CellNumber_Popupmenu.
function CellNumber_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to CellNumber_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CellNumber_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CellNumber_Popupmenu

handles.CellNumber = handles.CellNumberList(get(handles.CellNumber_Popupmenu,'Value'));

handles.TrialNumberList = unique(handles.data.Sweep_Counter(strcmp(handles.data.Mouse_Name,handles.MouseName) & handles.data.Cell_Counter == handles.CellNumber));
set(handles.TrialNumber_Popupmenu,'String',handles.TrialNumberList);
handles.TrialNumber = handles.TrialNumberList(1);
set(handles.TrialNumber_Popupmenu,'Value',1);

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CellNumber_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellNumber_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TrialNumber_Popupmenu.
function TrialNumber_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to TrialNumber_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TrialNumber_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TrialNumber_Popupmenu

handles.TrialNumber = handles.TrialNumberList(get(handles.TrialNumber_Popupmenu,'Value'));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TrialNumber_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrialNumber_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% --- Change Axes Scale --- 
%
function WhiskerAxisMinTag_Callback(hObject, eventdata, handles)
% hObject    handle to WhiskerAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhiskerAxisMinTag as text
%        str2double(get(hObject,'String')) returns contents of WhiskerAxisMinTag as a double

handles.WhiskerAxisMin = round(str2double(get(handles.WhiskerAxisMinTag,'String')));

% Change the text
set(handles.WhiskerAxisMinTag,'String',num2str(handles.WhiskerAxisMin));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function WhiskerAxisMinTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiskerAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function WhiskerAxisMaxTag_Callback(hObject, eventdata, handles)
% hObject    handle to WhiskerAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhiskerAxisMaxTag as text
%        str2double(get(hObject,'String')) returns contents of WhiskerAxisMaxTag as a double

handles.WhiskerAxisMax = round(str2double(get(handles.WhiskerAxisMaxTag,'String')));

% Change the text
set(handles.WhiskerAxisMaxTag,'String',num2str(handles.WhiskerAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function WhiskerAxisMaxTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiskerAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function TimeAxisMinTag_Callback(hObject, eventdata, handles)
% hObject    handle to TimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeAxisMinTag as text
%        str2double(get(hObject,'String')) returns contents of TimeAxisMinTag as a double

handles.TimeAxisMin = round(str2double(get(handles.TimeAxisMinTag,'String')));

% Change the text
set(handles.TimeAxisMinTag,'String',num2str(handles.TimeAxisMin));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TimeAxisMinTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function TimeAxisMaxTag_Callback(hObject, eventdata, handles)
% hObject    handle to TimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeAxisMaxTag as text
%        str2double(get(hObject,'String')) returns contents of TimeAxisMaxTag as a double

handles.TimeAxisMax = round(str2double(get(handles.TimeAxisMaxTag,'String')));

% Change the text
set(handles.TimeAxisMaxTag,'String',num2str(handles.TimeAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TimeAxisMaxTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MembranePotentialAxisMinTag_Callback(hObject, eventdata, handles)
% hObject    handle to MembranePotentialAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MembranePotentialAxisMinTag as text
%        str2double(get(hObject,'String')) returns contents of MembranePotentialAxisMinTag as a double

handles.MembranePotentialAxisMin = round(str2double(get(handles.MembranePotentialAxisMinTag,'String')));

% Change the text
set(handles.MembranePotentialAxisMinTag,'String',num2str(handles.MembranePotentialAxisMin));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MembranePotentialAxisMinTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MembranePotentialAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MembranePotentialAxisMaxTag_Callback(hObject, eventdata, handles)
% hObject    handle to MembranePotentialAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MembranePotentialAxisMaxTag as text
%        str2double(get(hObject,'String')) returns contents of MembranePotentialAxisMaxTag as a double

handles.MembranePotentialAxisMax = round(str2double(get(handles.MembranePotentialAxisMaxTag,'String')));

% Change the text
set(handles.MembranePotentialAxisMaxTag,'String',num2str(handles.MembranePotentialAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MembranePotentialAxisMaxTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MembranePotentialAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ContactTimeAxisMinTag_Callback(hObject, eventdata, handles)
% hObject    handle to ContactTimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ContactTimeAxisMinTag as text
%        str2double(get(hObject,'String')) returns contents of ContactTimeAxisMinTag as a double

handles.ContactTimeAxisMin = str2double(get(handles.ContactTimeAxisMinTag,'String'));

% Change the text
set(handles.ContactTimeAxisMinTag,'String',num2str(handles.ContactTimeAxisMin));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ContactTimeAxisMinTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ContactTimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ContactTimeAxisMaxTag_Callback(hObject, eventdata, handles)
% hObject    handle to ContactTimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ContactTimeAxisMaxTag as text
%        str2double(get(hObject,'String')) returns contents of ContactTimeAxisMaxTag as a double

handles.ContactTimeAxisMax = str2double(get(handles.ContactTimeAxisMaxTag,'String'));

% Change the text
set(handles.ContactTimeAxisMaxTag,'String',num2str(handles.ContactTimeAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ContactTimeAxisMaxTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ContactTimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function WhiskingTimeAxisMinTag_Callback(hObject, eventdata, handles)
% hObject    handle to WhiskingTimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhiskingTimeAxisMinTag as text
%        str2double(get(hObject,'String')) returns contents of WhiskingTimeAxisMinTag as a double

handles.WhiskingTimeAxisMin = str2double(get(handles.WhiskingTimeAxisMinTag,'String'));

% Change the text
set(handles.WhiskingTimeAxisMinTag,'String',num2str(handles.WhiskingTimeAxisMin));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function WhiskingTimeAxisMinTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiskingTimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function WhiskingTimeAxisMaxTag_Callback(hObject, eventdata, handles)
% hObject    handle to WhiskingTimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WhiskingTimeAxisMaxTag as text
%        str2double(get(hObject,'String')) returns contents of WhiskingTimeAxisMaxTag as a double

handles.WhiskingTimeAxisMax = str2double(get(handles.WhiskingTimeAxisMaxTag,'String'));

% Change the text
set(handles.WhiskingTimeAxisMaxTag,'String',num2str(handles.WhiskingTimeAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function WhiskingTimeAxisMaxTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WhiskingTimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function QuietTimeAxisMinTag_Callback(hObject, eventdata, handles)
% hObject    handle to QuietTimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of QuietTimeAxisMinTag as text
%        str2double(get(hObject,'String')) returns contents of QuietTimeAxisMinTag as a double

handles.QuietTimeAxisMin = str2double(get(handles.QuietTimeAxisMinTag,'String'));

% Change the text
set(handles.QuietTimeAxisMinTag,'String',num2str(handles.QuietTimeAxisMin));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function QuietTimeAxisMinTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QuietTimeAxisMinTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function QuietTimeAxisMaxTag_Callback(hObject, eventdata, handles)
% hObject    handle to QuietTimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of QuietTimeAxisMaxTag as text
%        str2double(get(hObject,'String')) returns contents of QuietTimeAxisMaxTag as a double

handles.QuietTimeAxisMax = str2double(get(handles.QuietTimeAxisMaxTag,'String'));

% Change the text
set(handles.QuietTimeAxisMaxTag,'String',num2str(handles.QuietTimeAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function QuietTimeAxisMaxTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QuietTimeAxisMaxTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ResetScalingPushbutton.
function ResetScalingPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ResetScalingPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.MembranePotentialAxisMin = -0.080;
handles.MembranePotentialAxisMax = 0.040;
handles.WhiskerAxisMin = 175;
handles.WhiskerAxisMax = 250;
handles.WhiskingTimeAxisMin = -0.2;
handles.WhiskingTimeAxisMax = 1.2;
handles.QuietTimeAxisMin = -0.2;
handles.QuietTimeAxisMax = 1.2;
handles.ContactTimeAxisMin = -0.2;
handles.ContactTimeAxisMax = 1.2;
handles.TimeAxisMin = 0;
handles.TimeAxisMax = 3000;

% Change the text
set(handles.MembranePotentialAxisMinTag,'String',num2str(handles.MembranePotentialAxisMin));
set(handles.MembranePotentialAxisMaxTag,'String',num2str(handles.MembranePotentialAxisMax));
set(handles.WhiskerAxisMinTag,'String',num2str(handles.WhiskerAxisMin));
set(handles.WhiskerAxisMaxTag,'String',num2str(handles.WhiskerAxisMax));
set(handles.WhiskingTimeAxisMinTag,'String',num2str(handles.WhiskingTimeAxisMin));
set(handles.WhiskingTimeAxisMaxTag,'String',num2str(handles.WhiskingTimeAxisMax));
set(handles.QuietTimeAxisMinTag,'String',num2str(handles.QuietTimeAxisMin));
set(handles.QuietTimeAxisMaxTag,'String',num2str(handles.QuietTimeAxisMax));
set(handles.ContactTimeAxisMinTag,'String',num2str(handles.ContactTimeAxisMin));
set(handles.ContactTimeAxisMaxTag,'String',num2str(handles.ContactTimeAxisMax));
set(handles.TimeAxisMinTag,'String',num2str(handles.TimeAxisMin));
set(handles.TimeAxisMaxTag,'String',num2str(handles.TimeAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in MoveBackButton.
function MoveBackButton_Callback(hObject, eventdata, handles)
% hObject    handle to MoveBackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
move_time = (handles.TimeAxisMax - handles.TimeAxisMin) * 0.8;
handles.TimeAxisMin = handles.TimeAxisMin - move_time;
handles.TimeAxisMax = handles.TimeAxisMax - move_time;
set(handles.TimeAxisMinTag,'String',num2str(handles.TimeAxisMin));
set(handles.TimeAxisMaxTag,'String',num2str(handles.TimeAxisMax));
handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in MoveForwardsButton.
function MoveForwardsButton_Callback(hObject, eventdata, handles)
% hObject    handle to MoveForwardsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
move_time = (handles.TimeAxisMax - handles.TimeAxisMin) * 0.8;
handles.TimeAxisMin = handles.TimeAxisMin + move_time;
handles.TimeAxisMax = handles.TimeAxisMax + move_time;
set(handles.TimeAxisMinTag,'String',num2str(handles.TimeAxisMin));
set(handles.TimeAxisMaxTag,'String',num2str(handles.TimeAxisMax));
handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in AutoscaleButton.
function AutoscaleButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutoscaleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.MembranePotentialAxisMin = min(handles.MembranePotential);
handles.MembranePotentialAxisMax = max(handles.MembranePotential);
handles.WhiskerAxisMin = min(handles.WhiskerAngle);
handles.WhiskerAxisMax = max(handles.WhiskerAngle);
handles.WhiskingTimeAxisMin = -0.2;
handles.WhiskingTimeAxisMax = 1.2;
handles.QuietTimeAxisMin = -0.2;
handles.QuietTimeAxisMax = 1.2;
handles.ContactTimeAxisMin = -0.2;
handles.ContactTimeAxisMax = 1.2;
handles.TimeAxisMin = min(handles.MembranePotentialTimeVector);
handles.TimeAxisMax = max(handles.MembranePotentialTimeVector);

% Change the text
set(handles.MembranePotentialAxisMinTag,'String',num2str(handles.MembranePotentialAxisMin));
set(handles.MembranePotentialAxisMaxTag,'String',num2str(handles.MembranePotentialAxisMax));
set(handles.WhiskerAxisMinTag,'String',num2str(handles.WhiskerAxisMin));
set(handles.WhiskerAxisMaxTag,'String',num2str(handles.WhiskerAxisMax));
set(handles.WhiskingTimeAxisMinTag,'String',num2str(handles.WhiskingTimeAxisMin));
set(handles.WhiskingTimeAxisMaxTag,'String',num2str(handles.WhiskingTimeAxisMax));
set(handles.QuietTimeAxisMinTag,'String',num2str(handles.QuietTimeAxisMin));
set(handles.QuietTimeAxisMaxTag,'String',num2str(handles.QuietTimeAxisMax));
set(handles.ContactTimeAxisMinTag,'String',num2str(handles.ContactTimeAxisMin));
set(handles.ContactTimeAxisMaxTag,'String',num2str(handles.ContactTimeAxisMax));
set(handles.TimeAxisMinTag,'String',num2str(handles.TimeAxisMin));
set(handles.TimeAxisMaxTag,'String',num2str(handles.TimeAxisMax));

handles = UpdateTimeSeries(handles);

% Update handles structure
guidata(hObject, handles);
