function varargout = GUI_Cuff(varargin)
% This code calculates capacitance, charge transfer resistance, and
% solution resistance of cuff and solution.
% AccRCT and AccIMP show the Rct and capacitacne based on each frequecy
% calculation; phase and gaind calculation. Rct and Caoacitor show the
% results with an estimation using a high frequency signal.
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Cuff_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Cuff_OutputFcn, ...
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

% --- Executes just before GUI_Cuff is made visible.
function GUI_Cuff_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

set(handles.figure1, 'pointer', 'arrow')

axes(handles.axes4);
imshow(imread('CuffSchematic.png')); % Showing the schematic

% Default values
set(handles.MSR,'string','1500000'); % Maximum Frequency of DAQ Board
set(handles.MiF,'string','1'); % Minimum Frequency for measurement
set(handles.MaF,'string','10000'); % Maximum Frequency for measurement
set(handles.Re,'string','10000'); % External resistance
set(handles.LBS,'string','9'); % Log based steps
set(handles.IC0,'string','ai0'); % Input channel of DAQ board which should be connected to Output channel as follow
set(handles.IC1,'string','ai2'); % Input channel of DAQ board which should be connected to cuff signal
set(handles.OC,'string','ao0'); % Output channel of DAQ board which should be connected to cuff signal
set(handles.Amp,'string','1'); % Amplitude of signals
set(handles.Time,'string','5'); % Time duration of signals
set(handles.GF,'string','0.01'); % Gain Factor for high frequency (infinit frequency)
set(handles.CPE,'string','1'); % CPE factor

set(handles.CPEFitting, 'enable', 'off'); % CPE fitting tools

% Electrode methods
a = 'Dipole (Cuff) Magnitude & Phase';
b = 'Electrode Magnitude & Phase in Symetrical Dipole';
c = 'Dipole (Cuff) Estimated';
d = 'Dipole (Cuff) Two Equations';
e = 'Electrode Estimated in Symetrical Dipole';
s = char(a, b, c, d, e);
set(handles.Electrode_Method,'string',s)


% Fitting methods
c = 'Signle Frequency';
b = 'A Range of Frequency';
a = 'Best Fitting Factors';
s = char(a, b, c);
set(handles.popupmenu_fitting,'string',s)


set(handles.popupmenu_fitting, 'enable', 'off')


handles.FileName_String = 'Cuff.mat';
set(handles.FileName, 'string', handles.FileName_String)

handles.MSR_Value = get(handles.MSR, 'String');
handles.MiF_Value = get(handles.MiF, 'String');
handles.MaF_Value = get(handles.MaF, 'String');
handles.LBS_Value = get(handles.LBS, 'String');
handles.OC_Value = get(handles.OC, 'String');
handles.Re_Value = get(handles.Re, 'String');
handles.CPE_Value = get(handles.CPE, 'String');

handles.panel1 = get(handles.uipanel2, 'Title');
handles.panel2 = get(handles.uipanel3, 'Title');
handles.panel3 = get(handles.uipanel4, 'Title');
handles.panel4 = get(handles.uipanel6, 'Title');
handles.panel5 = get(handles.uipanel8, 'Title');
handles.panel6 = get(handles.uipanel5, 'Title');
handles.panel7 = get(handles.uipanel23, 'Title');


% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Cuff_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;




% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
set(handles.figure1, 'pointer', 'watch')
axes(handles.axes1);cla reset
axes(handles.axes2);cla reset
axes(handles.axes3);cla reset


max_sample_rate = str2double(get(handles.MSR,'String')); % Nominal max sample rate for daq board when using more than one channel


ai1 = get(handles.IC0,'String'); % ai0 to get a feedback from the generated signal
ai2 = get(handles.IC1,'String'); % ai1 for checking cuff signal
ao = get(handles.OC,'String'); % ao0 for the generating signal
ai1 = str2double(ai1(3:end));
ai2 = str2double(ai2(3:end));
ao = str2double(ao(3:end));


GainFactor = str2double(get(handles.GF,'String'));
Time = str2double(get(handles.Time,'String')); % Acquiring the signal for 5 secs
Amp = str2double(get(handles.Amp,'String')); % Setting all voltages to 10 volts
R = str2double(get(handles.Re,'String')); % This is the series resistance in circuit for cuff testing; this resistor should be 220K
CPE = str2double(get(handles.CPE,'String')); % CPE factor
MaxFre = str2double(get(handles.MaF,'String')); % Maximum frequency for acquiring signal
MinFre = str2double(get(handles.MiF,'String')); % Minimum frequency for acquiring signal
LBS = str2double(get(handles.LBS,'String')); % Log based steps


% Check if calibration is availabe or not if does not exist then consider conditions ideal
if (~isfield(handles, 'Output_Data')   ||  ~isfield(handles.Output_Data, 'Calibration'))
    handles.Output_Data.Calibration.Vdc_ch1(1) = 5;
    handles.Output_Data.Calibration.Vdc_ch2(1) = 5;
    handles.Output_Data.Calibration.Vdc(1) = 5;
    
    handles.Output_Data.Calibration.Vdc_ch1(2) = 0;
    handles.Output_Data.Calibration.Vdc_ch2(2) = 0;
    handles.Output_Data.Calibration.Vdc(2) = 0;
end
% Update the Calibration parameters
Calibration = handles.Output_Data.Calibration;


% checking which mode of operation has been asked for processing from pop
% up bottom for selecting emthod
if get(handles.Electrode_Method,'Value') == 3 % Cuff estimated
    Output = Cuff_Measurment(max_sample_rate, MaxFre, MinFre, LBS, ai1, ai2, ao, Time, Amp, R, GainFactor, CPE, 0, Calibration);
    
elseif get(handles.Electrode_Method,'Value') == 4 % Cuff Two equations
    Fre = str2double(get(handles.LBS,'String'));
    Output = Cuff_Measurment(max_sample_rate, MaxFre, MinFre, LBS, ai1, ai2, ao, Time, Amp, R, GainFactor, CPE, 1, Calibration);
    
elseif get(handles.Electrode_Method,'Value') == 5 % Electrode estimated
    Output = Electrode_Measurment(max_sample_rate, MaxFre, MinFre, LBS, ai1, ai2, ao, Time, Amp, R, GainFactor, CPE, 0, Calibration);
    
elseif get(handles.Electrode_Method,'Value') == 1 || get(handles.Electrode_Method,'Value') == 2% Just Impedance and Phase Measurement
    Output = Just_Impednace_Phase(max_sample_rate, MaxFre, MinFre, LBS, ai1, ai2, ao, Time, Amp, R, Calibration);
    
end


FREQUENCY = Output.FREQUENCY;
Capacitor = Output.Capacitor;
Rct = Output.Rct;
Rs = Output.Rs;
Measured_Imp = Output.Measured_Imp;
Measured_Phase = Output.Measured_Phase;
Calculated_Impedance = Output.Calculated_Impedance;
Calculated_Phase = Output.Calculated_Phase;
CPE_Factor = Output.CPE_Factor;
    

c ={'red', 'blue', 'green', 'yellow', 'magenta', 'cyan'};
% plotting  
if get(handles.Electrode_Method,'Value') == 1 || get(handles.Electrode_Method,'Value') == 2
    % Plot C, Impedance Magnitude, Real Impedance, and error based on current
    PLOT_IMP_PHASE(handles, FREQUENCY, Measured_Imp, Measured_Phase, c);
else
    % Plot C, Impedance Magnitude, Real Impedance, and error based on current
    PLOT_C_IMP(handles, FREQUENCY, Capacitor, Measured_Imp, Measured_Phase, Calculated_Impedance, Calculated_Phase, c);

    % Charge Transfer resistance
    PLOT_Rct(handles, FREQUENCY, Rct, c)

    % Solution Resistance
    PLOT_Rs(handles, FREQUENCY, Rs, c)

    % Phase
    PLOT_Phase(handles, FREQUENCY, Calculated_Phase, Measured_Phase, c, CPE_Factor)
end


% Put variables in handles for updating
handles.Output_Data.Mode = string(handles.Electrode_Method.String(get(handles.Electrode_Method,'Value'),:));
handles.Output_Data.FREQUENCY = FREQUENCY;
handles.Output_Data.Capacitor = Capacitor;
handles.Output_Data.Rct = Rct;
handles.Output_Data.Rs = Rs;
handles.Output_Data.Measured_Imp = Measured_Imp;
handles.Output_Data.Measured_Phase = Measured_Phase;
handles.Output_Data.Calculated_Impedance = Calculated_Impedance;
handles.Output_Data.Calculated_Phase = Calculated_Phase;
handles.Output_Data.CPE_Factor = CPE_Factor;
handles.Output_Data.Output = Output;


% put notification to save the data or reset if needed
set(handles.Run,'Backgroundcolor','r');
set(handles.DataError,'string','Data is not Saved');
set(handles.Run, 'enable', 'off');

set(handles.figure1, 'pointer', 'arrow')

guidata(hObject, handles);


% --- Executes on button press in CPEFitting.
function CPEFitting_Callback(hObject, eventdata, handles)
% hObject    handle to CPEFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 % CPE fitting tools
set(handles.figure1, 'pointer', 'watch')

if get(handles.Electrode_Method,'Value') == 1 || get(handles.Electrode_Method,'Value') == 3 || get(handles.Electrode_Method,'Value') == 4 
    Electrode_mode = 0;
else
    Electrode_mode = 1;
end

set(handles.figure1, 'pointer', 'watch')
axes(handles.axes1);cla reset
axes(handles.axes2);cla reset
axes(handles.axes3);cla reset

if get(handles.popupmenu_fitting,'Value') == 3
    Inputs = [str2double(handles.MSR.String), str2double(handles.MiF.String), str2double(handles.MaF.String)];

elseif get(handles.popupmenu_fitting,'Value') == 2
    Inputs = [str2double(handles.MSR.String), str2double(handles.MiF.String), str2double(handles.MaF.String), str2double(handles.LBS.String), str2double(handles.OC.String)];
    handles.Output_Data.Rct(:) = str2double(handles.Re.String);
    handles.Output_Data.Rs(:) = str2double(handles.CPE.String);
    
elseif get(handles.popupmenu_fitting,'Value') == 1
    Inputs = [str2double(handles.MSR.String), str2double(handles.MiF.String), str2double(handles.MaF.String), str2double(handles.LBS.String), str2double(handles.OC.String), 1];
    handles.Output_Data.Rct(:) = str2double(handles.Re.String);
    handles.Output_Data.Rs(:) = str2double(handles.CPE.String);
    
end


Output = Fitting(handles.Output_Data, Inputs, Electrode_mode);


FREQUENCY = Output.FREQUENCY;
Capacitor = Output.Capacitor;
Rct = Output.Rct;
Rs = Output.Rs;
Measured_Imp = Output.Measured_Imp;
Measured_Phase = Output.Measured_Phase;
Calculated_Impedance = Output.Calculated_Impedance;
Calculated_Phase = Output.Calculated_Phase;
CPE_Factor = Output.CPE_Factor;


c ={'red', 'blue', 'green', 'yellow', 'magenta', 'cyan'};
% C, Impedance Magnitude, Real Impedance, and error based on current
PLOT_C_IMP(handles, FREQUENCY, Capacitor, Measured_Imp, Measured_Phase, Calculated_Impedance, Calculated_Phase, c);

% Charge Transfer resistance
PLOT_Rct(handles, FREQUENCY, Rct, c)

% Solution Resistance
PLOT_Rs(handles, FREQUENCY, Rs, c)

% Phase
PLOT_Phase(handles, FREQUENCY, Calculated_Phase, Measured_Phase, c, CPE_Factor)

set(handles.DataError,'string','Data is not Saved');
set(handles.figure1, 'pointer', 'arrow')


if get(handles.Keep_Fitting,'Value') == 1
    Keep_Fitting = 1;
else
    set(handles.CPEFitting, 'enable', 'off');
    set(handles.popupmenu_fitting, 'enable', 'off')
    set(handles.MSR,'enable', 'on'); % Maximum Frequency of DAQ Board
    set(handles.MiF,'enable', 'on'); % Minimum Frequency for measurement
    set(handles.MaF,'enable', 'on'); % Maximum Frequency for measurement
    set(handles.LBS,'enable', 'on'); % Log based steps
    set(handles.IC0,'enable', 'on'); % Input channel of DAQ board which should be connected to Output channel as follow
    set(handles.IC1,'enable', 'on'); % Input channel of DAQ board which should be connected to cuff signal
    set(handles.OC,'enable', 'on'); % Output channel of DAQ board which should be connected to cuff signal
    set(handles.Amp,'enable', 'on'); % Amplitude of signals
    set(handles.Time,'enable', 'on'); % Time duration of signals
    set(handles.Re,'enable', 'on'); % Re
    set(handles.CPE,'enable', 'on'); % CPE
    set(handles.Electrode_Method,'enable', 'on'); % Method
    set(handles.GF,'enable', 'on'); % Infinite gain factor


    set(handles.MSR, 'String', handles.MSR_Value);
    set(handles.MiF, 'String', handles.MiF_Value);
    set(handles.MaF, 'String', handles.MaF_Value);
    set(handles.LBS, 'String', handles.LBS_Value);
    set(handles.OC, 'String', handles.OC_Value);
    set(handles.Re, 'String', handles.Re_Value);
    set(handles.CPE, 'String', handles.CPE_Value);

    set(handles.uipanel2, 'Title', handles.panel1);
    set(handles.uipanel3, 'Title', handles.panel2);
    set(handles.uipanel4, 'Title', handles.panel3);
    set(handles.uipanel6, 'Title', handles.panel4);
    set(handles.uipanel8, 'Title', handles.panel5);
    set(handles.uipanel5, 'Title', handles.panel6);
    set(handles.uipanel23, 'Title', handles.panel7);
end

handles.Output_Data.Mode = string(handles.Electrode_Method.String(get(handles.Electrode_Method,'Value'),:));
handles.Output_Data.FREQUENCY = FREQUENCY;
handles.Output_Data.Capacitor = Capacitor;
handles.Output_Data.Rct = Rct;
handles.Output_Data.Rs = Rs;
handles.Output_Data.Measured_Imp = Measured_Imp;
handles.Output_Data.Measured_Phase = Measured_Phase;
handles.Output_Data.Calculated_Impedance = Calculated_Impedance;
handles.Output_Data.Calculated_Phase = Calculated_Phase;
handles.Output_Data.Output = Output.Output;
handles.Output_Data.CPE_Factor = CPE_Factor;


guidata(hObject, handles);


% --- Executes on button press in LFF. load from file
function LFF_Callback(hObject, eventdata, handles)
set(handles.figure1, 'pointer', 'watch')

axes(handles.axes1);cla reset
axes(handles.axes2);cla reset
axes(handles.axes3);cla reset
set(handles.Rs_Tag,'string','');
set(handles.Rct_Tag,'string','');


[FileName,PathName] = uigetfile('*.mat','Select the MATLAB code file');
load(fullfile(PathName,FileName))
handles.FileName_String = FileName;
set(handles.FileName, 'string', handles.FileName_String)

Output = Output_Data.Output;

FREQUENCY = Output_Data.FREQUENCY;
Capacitor = Output_Data.Capacitor;
Rct = Output_Data.Rct;
Rs = Output_Data.Rs;
Measured_Imp = Output_Data.Measured_Imp;
Measured_Phase = Output_Data.Measured_Phase;
Calculated_Impedance = Output_Data.Calculated_Impedance;
Calculated_Phase = Output_Data.Calculated_Phase;
CPE_Factor = Output_Data.CPE_Factor;
    

handles.Output_Data.FREQUENCY = FREQUENCY;
handles.Output_Data.Capacitor = Capacitor;
handles.Output_Data.Rct = Rct;
handles.Output_Data.Rs = Rs;
handles.Output_Data.Measured_Imp = Measured_Imp;
handles.Output_Data.Measured_Phase = Measured_Phase;
handles.Output_Data.Calculated_Impedance = Calculated_Impedance;
handles.Output_Data.Calculated_Phase = Calculated_Phase;
handles.Output_Data.Output = Output;
handles.Output_Data.CPE_Factor= CPE_Factor;


c ={'red', 'blue', 'green', 'yellow', 'magenta', 'cyan'};
% C, Impedance Magnitude, Real Impedance, and error based on current
PLOT_C_IMP(handles, FREQUENCY, Capacitor, Measured_Imp, Measured_Phase, Calculated_Impedance, Calculated_Phase, c);

% Charge Transfer resistance
PLOT_Rct(handles, FREQUENCY, Rct, c)

% Solution Resistance
PLOT_Rs(handles, FREQUENCY, Rs, c)

% Phase
PLOT_Phase(handles, FREQUENCY, Calculated_Phase, Measured_Phase, c, CPE_Factor)

set(handles.figure1, 'pointer', 'arrow')

set(handles.popupmenu_fitting, 'enable', 'on')

guidata(hObject, handles);


% --- Executes on button press in ResetButton.
function ResetButton_Callback(hObject, eventdata, handles)
set(handles.Run,'Backgroundcolor',[0.94 0.94 0.94]);
set(handles.DataError,'string','');
set(handles.Rs_Tag,'string','');
set(handles.Rct_Tag,'string','');

axes(handles.axes1);cla reset
axes(handles.axes2);cla reset
axes(handles.axes3);cla reset
set(handles.Run, 'enable', 'on');

set(handles.figure1, 'pointer', 'arrow');

set(handles.MSR,'enable', 'on'); % Maximum Frequency of DAQ Board
set(handles.MiF,'enable', 'on'); % Minimum Frequency for measurement
set(handles.MaF,'enable', 'on'); % Maximum Frequency for measurement
set(handles.Re,'enable', 'on'); % External resistance
set(handles.LBS,'enable', 'on'); % Log based steps
set(handles.IC0,'enable', 'on'); % Input channel of DAQ board which should be connected to Output channel as follow
set(handles.IC1,'enable', 'on'); % Input channel of DAQ board which should be connected to cuff signal
set(handles.OC,'enable', 'on'); % Output channel of DAQ board which should be connected to cuff signal
set(handles.Amp,'enable', 'on'); % Amplitude of signals
set(handles.Time,'enable', 'on'); % Time duration of signals


set(handles.MSR, 'String', handles.MSR_Value);
set(handles.MiF, 'String', handles.MiF_Value);
set(handles.MaF, 'String', handles.MaF_Value);
set(handles.LBS, 'String', handles.LBS_Value);
set(handles.OC, 'String', handles.OC_Value);
set(handles.Re, 'String', handles.Re_Value);
set(handles.CPE, 'String', handles.CPE_Value);

set(handles.uipanel2, 'Title', handles.panel1);
set(handles.uipanel3, 'Title', handles.panel2);
set(handles.uipanel4, 'Title', handles.panel3);
set(handles.uipanel6, 'Title', handles.panel4);
set(handles.uipanel8, 'Title', handles.panel5);
set(handles.uipanel5, 'Title', handles.panel6);
set(handles.uipanel23, 'Title', handles.panel7);


set(handles.Cal,'Backgroundcolor',[0.94 0.94 0.94]);

set(handles.CPEFitting, 'enable', 'off'); % CPE fitting tools


% --- Executes on selection change in popupmenu_fitting.
function popupmenu_fitting_Callback(hObject, eventdata, handles)

set(handles.CPEFitting, 'enable', 'on')

set(handles.MSR,'enable', 'off'); % Maximum Frequency of DAQ Board
set(handles.MiF,'enable', 'off'); % Minimum Frequency for measurement
set(handles.MaF,'enable', 'off'); % Maximum Frequency for measurement
set(handles.LBS,'enable', 'off'); % Log based steps
set(handles.IC0,'enable', 'off'); % Input channel of DAQ board which should be connected to Output channel as follow
set(handles.IC1,'enable', 'off'); % Input channel of DAQ board which should be connected to cuff signal
set(handles.OC,'enable', 'off'); % Output channel of DAQ board which should be connected to cuff signal
set(handles.Amp,'enable', 'off'); % Amplitude of signals
set(handles.Time,'enable', 'off'); % Time duration of signals
set(handles.Re,'enable', 'off'); % Frequency for fast method
set(handles.CPE,'enable', 'off'); % Frequency for fast method
set(handles.Electrode_Method,'enable', 'on'); % Frequency for fast method
set(handles.GF,'enable', 'off'); % Frequency for fast method


set(handles.MSR,'enable', 'on'); % Acceptable Phase Error
set(handles.MiF,'enable', 'on'); % Phase Factor
set(handles.MaF,'enable', 'on'); % Magnitude Factor
    
set(handles.uipanel2, 'Title', 'Min Acceptable Phase Error');
set(handles.uipanel3, 'Title', 'Phase Error Factor');
set(handles.uipanel4, 'Title', 'Magnitude Error Factor');
set(handles.uipanel6, 'Title', 'Minimum Frequency');
set(handles.uipanel8, 'Title', 'Maximum Frequency');
set(handles.uipanel5, 'Title', 'Given Rct');
set(handles.uipanel23, 'Title', 'Given Rs');

set(handles.MSR, 'String', num2str(5));
set(handles.MiF, 'String', num2str(1));
set(handles.MaF, 'String', num2str(1));

if  isnan(nanmean(handles.Output_Data.Rct))
    set(handles.Re, 'String', num2str(100000));
    set(handles.CPE, 'String', num2str(20));
else
    
    set(handles.Re, 'String', num2str(nanmean(handles.Output_Data.Rct)));
    set(handles.CPE, 'String', num2str(nanmean(handles.Output_Data.Rs)));
end

set(handles.Re,'enable', 'on'); % Output Channel
set(handles.CPE,'enable', 'on'); % Re
    

if get(handles.popupmenu_fitting,'Value') == 2 || get(handles.popupmenu_fitting,'Value') == 1
    set(handles.LBS,'enable', 'on'); % Log based steps
    set(handles.LBS, 'String', num2str(1));
    set(handles.OC,'enable', 'on'); % Log based steps
    set(handles.OC, 'String', num2str(10000));
end
guidata(hObject, handles);


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)

FileName = get(handles.FileName,'string');
Path = uigetdir;
OutputPath = fullfile(Path,FileName);
Output_Data = handles.Output_Data;
save(OutputPath, 'Output_Data');

set(handles.Run,'Backgroundcolor',[0.94 0.94 0.94]);
set(handles.DataError,'string','');
set(handles.Run, 'enable', 'on');

axes(handles.axes1);cla reset
axes(handles.axes2);cla reset
axes(handles.axes3);cla reset


% --- Executes on button press in Cal.
function Cal_Callback(hObject, eventdata, handles)
set(handles.figure1, 'pointer', 'watch')

max_sample_rate = str2double(get(handles.MSR,'String'));
ai1 = get(handles.IC0,'String'); % ai0 to get a feedback from the generated signal
ai2 = get(handles.IC1,'String'); % ai1 for checking cuff signal
ao = get(handles.OC,'String'); % ao0 for the generating signal
ai1 = str2double(ai1(3:end));
ai2 = str2double(ai2(3:end));
ao = str2double(ao(3:end));

handles.Output_Data.Calibration = Calibration_Function(max_sample_rate, ai1, ai2, ao);

set(handles.Cal,'Backgroundcolor','g');

set(handles.figure1, 'pointer', 'arrow')

guidata(hObject, handles);


% --- Executes on selection change in Electrode_Method.
function Electrode_Method_Callback(hObject, eventdata, handles)
set(handles.MiF, 'enable', 'on');
set(handles.MaF, 'enable', 'on');
set(handles.LBS, 'enable', 'on');
axes(handles.axes4);

if get(handles.Electrode_Method,'Value') == 1 || get(handles.Electrode_Method,'Value') == 3 || get(handles.Electrode_Method,'Value') == 4
    imshow(imread('CuffSchematic.png')); % Showing the schematic
    
elseif get(handles.Electrode_Method,'Value') == 5 || get(handles.Electrode_Method,'Value') == 2
    imshow(imread('ElectrodeSchematic.png')); % Showing the schematic
end






% Buttons

function OC_Callback(hObject, eventdata, handles)
% hObject    handle to OC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function OC_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function IC0_Callback(hObject, eventdata, handles)
% hObject    handle to IC0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function IC0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LBS_Callback(hObject, eventdata, handles)
% hObject    handle to LBS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function LBS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Re_Callback(hObject, eventdata, handles)
% hObject    handle to Re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Re_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaF_Callback(hObject, eventdata, handles)
% hObject    handle to MaF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function MaF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MiF_Callback(hObject, eventdata, handles)
% hObject    handle to MiF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function MiF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MSR_Callback(hObject, eventdata, handles)
% hObject    handle to MSR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function MSR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_fitting_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Electrode_Method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Time_Callback(hObject, eventdata, handles)
% hObject    handle to Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Time_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Amp_Callback(hObject, eventdata, handles)
% hObject    handle to Amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Amp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function IC1_Callback(hObject, eventdata, handles)
% hObject    handle to IC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function IC1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FileName_Callback(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function FileName_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GF_Callback(hObject, eventdata, handles)
% hObject    handle to GF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function GF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CPE_Callback(hObject, eventdata, handles)
% hObject    handle to CPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function CPE_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Keep_Fitting_Callback(hObject, eventdata, handles)
% hObject    handle to CPE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)







% Plottings
function PLOT_C_IMP(handles, FREQUENCY, Capacitor, Measured_Imp, Measured_Phase, Calculated_Impedance, Calculated_Phase, c)
axes(handles.axes1);

% Impedance Magnitude
yyaxis left
E = plot(FREQUENCY(:,1), Calculated_Impedance(:,1), 'o-', 'MarkerSize',7);
E.Color = char(c(2));
List(1) = {'Fitted Impedance'};
E.LineWidth = 1.3;


% Real Impedance based on current
hold on
E = plot(FREQUENCY(:,1), Measured_Imp(:,1), 's-', 'MarkerSize',7); hold on
E.Color = char(c(3));
List(2) = {'Measured Impedance'};
E.LineWidth = 1.3;

yyaxis right
E = plot(FREQUENCY(:,1), Capacitor(:,1), '*-', 'MarkerSize',7); hold on
E.Color = char(c(1));
List(3) = {'CPE Value'};
E.LineWidth = 1.3;


legend(List,'FontSize',7,'fontweight','bold');
ax = gca; ax.XAxis.FontSize = 10; ax.XAxis.FontWeight = 'bold'; 
set(gca,'Color',[0.8 0.8 0.8]); xlabel('Frequency (Hz)'); 
yyaxis right; ylabel('CPE Value (sec^n * Ohm^-^1)'); yyaxis left; ylabel('Impedance Magnitude (Ohm)'); 
ax.YAxis(1).FontSize = 10; ax.YAxis(1).FontWeight = 'bold'; ax.YAxis(2).FontSize = 10; ax.YAxis(2).FontWeight = 'bold'; ax.XAxis.FontWeight = 'bold';
grid on; yyaxis left; set(ax,'XScale','log','YScale','log'); yyaxis right; set(ax,'YScale','log');
ax.YAxis(1).Color = char(c(2)); ax.YAxis(2).Color = char(c(1));


% Error
axes(handles.axes3);
List = {};
Error_Imp = abs(Calculated_Impedance - Measured_Imp)./Measured_Imp * 100;
yyaxis left
E = plot(FREQUENCY(:,1), Error_Imp(:,1), '*-', 'MarkerSize',7);
E.Color = char(c(2));
List(1) = {'Impedance Error'};
E.LineWidth = 1.3;



Error_Phase = abs(Calculated_Phase - Measured_Phase)*180/pi;
% Error
yyaxis right
E = plot(FREQUENCY(:,1), Error_Phase(:,1), 'o-', 'MarkerSize',7);
E.Color = char(c(1));
List(2) = {'Phase Error'};
E.LineWidth = 1.3;

legend(List,'FontSize',7,'fontweight','bold');
ax = gca; ax.XAxis.FontSize = 10; ax.XAxis.FontWeight = 'bold'; 
set(gca,'Color',[0.8 0.8 0.8]); xlabel('Frequency (Hz)'); ax.XAxis.FontWeight = 'bold';
yyaxis left; ylabel('Impedance Error (Percentage)'); yyaxis right; ylabel('Phase Error (Degree)'); 
ax.YAxis(1).FontSize = 10; ax.YAxis(1).FontWeight = 'bold'; ax.YAxis(2).FontSize = 10; ax.YAxis(2).FontWeight = 'bold';
grid on; yyaxis left; set(ax,'XScale','log'); yyaxis right;
ax.YAxis(1).Color = char(c(2)); ax.YAxis(2).Color = char(c(1));

function PLOT_Rct(handles, FREQUENCY, Rct, c)
% Rct
RCT = num2str(nanmean(Rct), '%.3e');
set(handles.Rct_Tag,'string',RCT);

function PLOT_Rs(handles, FREQUENCY, Rs, c)
% Rs
RS = num2str(nanmean(Rs), '%.3e');
set(handles.Rs_Tag,'string',RS);

function PLOT_Phase(handles, FREQUENCY, Calculated_Phase, Measured_Phase, c, CPE_Factor)
% Phase
axes(handles.axes2);
List = {};

yyaxis left
E = plot(FREQUENCY(:,1), Calculated_Phase(:,1)/pi*180, '*-', 'MarkerSize',6);
E.Color = char(c(2));
E.LineWidth = 1.3;
List(1) = {'Fitted Phase'};
hold on

E = plot(FREQUENCY(:,1), Measured_Phase(:,1)/pi*180, 'o-', 'MarkerSize',6);
E.Color = char(c(3));
E.LineWidth = 1.3;
List(2) = {'Measured Phase'};

yyaxis right
E = plot(FREQUENCY(:,1), CPE_Factor(:,1), 's-', 'MarkerSize',6);
E.Color = char(c(1));
E.LineWidth = 1.3;
List(3) = {'Fitted or Given CPE'};


legend(List,'FontSize',6,'fontweight','bold');
ax = gca; ax.XAxis.FontSize = 10; ax.XAxis.FontWeight = 'bold'; 
set(gca,'Color',[0.8 0.8 0.8]); xlabel('Frequency (Hz)'); ax.XAxis.FontWeight = 'bold';
yyaxis left; ylabel('Phase (Degree)'); yyaxis right; ylabel('CPE Factor'); 
ax.YAxis(1).FontSize = 10; ax.YAxis(1).FontWeight = 'bold'; ax.YAxis(2).FontSize = 10; ax.YAxis(2).FontWeight = 'bold';
grid on; set(ax,'XScale','log');
ax.YAxis(1).Color = char(c(2)); ax.YAxis(2).Color = char(c(1)); ax.YAxis(2).Limits = [0 1];

function PLOT_IMP_PHASE(handles, FREQUENCY, Measured_Imp, Measured_Phase, c)
axes(handles.axes1)

E = plot(FREQUENCY(:,1), Measured_Imp(:,1), 's-', 'MarkerSize',7); hold on
E.Color = char(c(3));
List(1) = {'Measured Impedance'};
E.LineWidth = 1.3;
legend(List,'FontSize',7,'fontweight','bold');

ax = gca; ax.XAxis.FontSize = 10; ax.XAxis.FontWeight = 'bold'; 
set(gca,'Color',[0.8 0.8 0.8]); xlabel('Frequency (Hz)'); 
ylabel('Impedance Magnitude (Ohm)'); 
ax.YAxis.FontSize = 10; ax.YAxis.FontWeight = 'bold'; ax.XAxis.FontWeight = 'bold';
grid on; set(ax,'XScale','log','YScale','log');


axes(handles.axes2)

E = plot(FREQUENCY(:,1), Measured_Phase(:,1)/pi*180, 's-', 'MarkerSize',7); hold on
E.Color = char(c(3));
List(1) = {'Measured Phase'};
E.LineWidth = 1.3;
legend(List,'FontSize',7,'fontweight','bold');

ax = gca; ax.XAxis.FontSize = 10; ax.XAxis.FontWeight = 'bold'; 
set(gca,'Color',[0.8 0.8 0.8]); xlabel('Frequency (Hz)'); 
ylabel('Phase (Degree)'); 
ax.YAxis.FontSize = 10; ax.YAxis.FontWeight = 'bold'; ax.XAxis.FontWeight = 'bold';
grid on; set(ax,'XScale','log');
