function Output = Just_Impednace_Phase(max_sample_rate, MaxFre, MinFre, LBS, ai1, ai2, ao, Time, Amp, R, Calibration)

% Get a copy of inputs for output saving
Output.Input_Parameters.max_sample_rate = max_sample_rate;
Output.Input_Parameters.MaxFre = MaxFre;
Output.Input_Parameters.MinFre = MinFre;
Output.Input_Parameters.LBS = LBS;
Output.Input_Parameters.ai1 = ai1;
Output.Input_Parameters.ai2 = ai2;
Output.Input_Parameters.ao = ao;
Output.Input_Parameters.Time = Time;
Output.Input_Parameters.Amp = Amp;
Output.Input_Parameters.R = R;

Devices = daq.getDevices;
Device_Name = Devices.ID;

%Create acquisition session
session = daq.createSession('ni');

sample_rate = max_sample_rate/5; % 3 channels are needed here
tot_samples = Time*sample_rate; % Total number of samples based on Time
%set acquisition parameters
session.NumberOfScans = tot_samples;
session.DurationInSeconds = Time;
session.Rate = sample_rate;
%set input channels
ch1 = addAnalogInputChannel(session, Device_Name, ai1, 'Voltage');
ch2 = addAnalogInputChannel(session, Device_Name, ai2, 'Voltage');
ch1.TerminalConfig = 'Differential';
ch2.TerminalConfig = 'Differential';
ch1.Range = [-Amp, Amp];
ch2.Range = [-Amp, Amp];
%set output channel
addAnalogOutputChannel(session, Device_Name, ao, 'Voltage');


Output.Frequencies = logspace(log10(MinFre), log10(MaxFre), LBS); % Frequencies to test the cuff

LEN = length(Output.Frequencies);
% Initialize the output matrix
Output.GAIN = nan(LEN,1);
Output.FREQUENCY = nan(LEN,1);
Output.Rs = nan(LEN,1);
Output.Rct = nan(LEN,1);
Output.Capacitor = nan(LEN,1);
Output.Vp_p1 = nan(LEN,1);
Output.Vp_p2 = nan(LEN,1);
Output.Measured_Imp = nan(LEN,1);
Output.Measured_Phase = nan(LEN,1);
Output.Calculated_Impedance = nan(LEN,1);
Output.Calculated_Phase = nan(LEN,1);
Output.CPE_Factor = nan(LEN,1);


for i = 1: length(Output.Frequencies)
    % Data acquisition for ac signals
    Frequency = Output.Frequencies(i);
    outputSignal1 = Amp * sin(Frequency * linspace(0, pi*2,session.Rate)'); % Sine wave with frequency based on the for loop
    outputSignal1 = cat(1, outputSignal1, [0; 0; 0]);
    
    queueOutputData(session,repmat(outputSignal1,Time,1));
    data = session.startForeground;
    
    Cutoff_fre = Frequency/(sample_rate/2)*10; % cutoff frequency for lowpass filter
    if Cutoff_fre>0.5 % just make sure that the frequncy is not passing 1
        Cutoff_fre = 0.5;
    end
    [k,l] = butter(3,Cutoff_fre);
    DATA = filtfilt(k,l,data(:,1));
    DATA = interp1(Calibration.Vdc_ch1, Calibration.Vdc, DATA,'linear', 'extrap');
    Output.DATA(i).data1 = DATA;
    Value1 = findpeaks(DATA); % finding the positive peaks in input signal
    Value1 = mean(Value1);
    Value2 = findpeaks((-1)*DATA); % finding the negative peaks in input signal
    Value2 = mean(Value2);
    Output.Vp_p1(i) = abs(Value1)+abs(Value2); % Volatge peak to peak calculation

    DATA = filtfilt(k,l,data(:,2));
    DATA = interp1(Calibration.Vdc_ch2, Calibration.Vdc, DATA,'linear', 'extrap');
    Output.DATA(i).data2 = DATA;
    Value1 = findpeaks(DATA); % finding the positive peaks in cuff signal
    Value1 = mean(Value1);
    Value2 = findpeaks((-1)*DATA); % finding the negative peaks in cuff signal
    Value2 = mean(Value2);
    Output.Vp_p2(i) = abs(Value1)+abs(Value2); % Volatge peak to peak calculation

    
    % finding Capacitor by estimation
    Gain = Output.Vp_p2(i)/Output.Vp_p1(i); % calculating gain for each frequency
    
    
    % finding Capacitor and RCT more accurate
    % find zero crossing
    [acor,lag] = xcorr(Output.DATA(i).data2, Output.DATA(i).data1);
    [~,I] = max(abs(acor));
    lagDiff = lag(I);
    
    Time_Diff = lagDiff/session.Rate;
    Period_Time = 1/Frequency;
    Phase_Diff = Time_Diff/Period_Time*2*pi;
    Phase_Diff = tan(Phase_Diff);
    
    
    % Calculation for measured magnitude and phase of impedance
    Betha = atan(-Output.Vp_p2(i)*sin(Phase_Diff)/(Output.Vp_p1(i)-Output.Vp_p2(i)*cos(Phase_Diff)));
    Re_Voltage = sqrt(Output.Vp_p2(i)^2+Output.Vp_p1(i)^2-2*Output.Vp_p1(i)*Output.Vp_p2(i)*cos(Phase_Diff));
    Output.Measured_Imp(i) = Output.Vp_p2(i)/(Re_Voltage/R);
    Output.Measured_Phase(i) = -wrapToPi(Phase_Diff-Betha);
    
    
    % Creating out put
    Output.GAIN(i) = Gain;
    Output.FREQUENCY(i) = Frequency;
end
