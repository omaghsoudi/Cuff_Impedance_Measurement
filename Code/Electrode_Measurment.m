function Output = Electrode_Measurment(max_sample_rate, MaxFre, MinFre, LBS, ai1, ai2, ao, Time, Amp, R, GainFactor, CPE, mode, Calibration)

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
Output.Input_Parameters.GainFactor = GainFactor;
Output.Input_Parameters.CPE = CPE;

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


% DC Test
outputSignal1 = Amp*ones(tot_samples,1); % Time secs DC Signal
outputSignal1 = cat(1,outputSignal1, [0; 0; 0]);
queueOutputData(session,repmat(outputSignal1,Time,1));
% session.startForeground;
% queueOutputData(session,repmat(outputSignal1,Time,1)); % Put the DC signal on Queue to generate the signal
data_dc = session.startForeground; % Generate and acquire signal for 5 seconds
Output.Vdc_in = interp1(Calibration.Vdc_ch1, Calibration.Vdc, mean(data_dc(:,1)),'linear', 'extrap'); % Input voltage
Output.Vdc_Cuff = interp1(Calibration.Vdc_ch2, Calibration.Vdc, mean(data_dc(:,2)),'linear', 'extrap'); % Output voltage

% DC zero level
outputSignal1 = zeros(tot_samples,1);
outputSignal1 = cat(1,outputSignal1, [0; 0; 0]);
queueOutputData(session,repmat(outputSignal1,Time,1));
session.startForeground; % Generate and acquire signal for 1 seconds


if mode == 0 || mode == 2
    % High Frequency test
    % Ideal case is having a device with a 1G samples per second and then the frequency of
    % this signal will be around 3MHz and then Rs and Rct will be calculated more accuarte 
    % but a total sample rate of 1,500,000 for all three channels is also acceptable
    % for just estimating capacitance, Rct, and Rs. It is enough to say that
    % cuff works properly or not.
    Frequency = max_sample_rate/50; % Time secs max AC Signal which looks like to a sine wave in output
    outputSignal1 = Amp*sin(Frequency*linspace(0,pi*2,session.Rate)'); % Sine wave
    outputSignal1 = cat(1,outputSignal1, [0; 0; 0]);
    queueOutputData(session,repmat(outputSignal1,Time,1)); 
    data = session.startForeground;
    [b,a] = butter(3,0.8); % lowpass filter with high cutoff frequency to just remove real fast signals
    DATA = filtfilt(b,a,data(:,1));
    Value1 = findpeaks(DATA); % finding the positive peaks in input signal
    Value1 = mean(Value1);
    Value2 = findpeaks((-1)*DATA); % finding the negative peaks in input signal
    Value2 = mean(Value2);
    Output.Vac_in = abs(Value1)+abs(Value2); % Volatge peak to peak calculation
    Output.Vac_in = interp1(Calibration.Vdc_ch1, Calibration.Vdc, Output.Vac_in,'linear', 'extrap');

    DATA = filtfilt(b,a,data(:,2));
    Value1 = findpeaks(DATA); % finding the positive peaks in cuff signal
    Value1 = mean(Value1);
    Value2 = findpeaks((-1)*DATA); % finding the negative peaks in cuff signal
    Value2 = mean(Value2);
    Output.Vac_Cuff = abs(Value1)+abs(Value2); % Volatge peak to peak calculation
    Output.Vac_Cuff = interp1(Calibration.Vdc_ch2, Calibration.Vdc, Output.Vac_Cuff,'linear', 'extrap');


    % Resistors Calculation
    %Estimation
    Output.Gain_Inf = Output.Vac_Cuff/Output.Vac_in;
    Rs = Output.Gain_Inf*GainFactor * R /(1-Output.Gain_Inf*GainFactor); % 100 here is just to have a better estimation for infinite ac voltage;
    Rs = Rs/2; %make it half as two are in series
    Output.Rs = Rs;
end


% however, the Gain_Inf can go down to even 1/1000 when the frequency goes to
% 3MHz or higher (Rs is estimated if the sample rate is not high enough)
Output.Gain_dc = Output.Vdc_Cuff/Output.Vdc_in;
RT = Output.Gain_dc * R / (1-Output.Gain_dc);
if mode == 0 % Estimated Cuff
    Rct = RT-Rs; % Rs is small and Rct is almost 1000 times larger; therefore, an estimation for Rs is not affecting Rct 
    Rct = Rct/2; % make it half as two are in series
    Output.Rct = Rct;
    
else % two equations
    %Accurate
    syms Rct positive
    Rs = (RT-2*Rct)/2;
end
    

% Swipping frequencies to find Capacitance
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


syms IMP positive % just positive response for IMP (capacitance) is acceptable
for i = 1: length(Output.Frequencies)
    % Data acquisition for ac signals
    Frequency = Output.Frequencies(i);
    outputSignal1 = Amp*sin(Frequency*linspace(0,pi*2,session.Rate)'); % Sine wave with frequency based on the for loop
    outputSignal1 = cat(1, outputSignal1, 0);
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
    f = Frequency;
    
    
    [acor,lag] = xcorr(Output.DATA(i).data2, Output.DATA(i).data1);
    [~,I] = max(abs(acor));
    lagDiff = lag(I);
    
    
    Time_Diff = lagDiff/session.Rate;
    Period_Time = 1/Frequency;
    Phase_Diff = Time_Diff/Period_Time*2*pi;
    Phase_Diff = tan(Phase_Diff);
    
    
    % a,b,c,d based on paper
    a = 2*Rct + 2*Rs + 2*Rct*Rs*IMP*(2*pi*f)^CPE*cos(CPE*pi/2);
    b = 2*Rct*Rs*IMP*(2*pi*f)^CPE*sin(CPE*pi/2);
    c = 2*Rct + (2*Rs + R)*(Rct*IMP*(2*pi*f)^CPE*cos(CPE*pi/2)+1);
    d = (2*Rs + R)*Rct*IMP*(2*pi*f)^CPE*sin(CPE*pi/2);
    
    
    
    Frequency_eqn = (Gain^2 * (c^2+d^2)^2 == (a*c+b*d)^2+(b*c-a*d)^2);% symbolic equation
    ResZ = solve(Frequency_eqn, IMP, 'Real', true); % Solver with just keeping the real answers
    
    if mode == 0
        %ResZ = solve(Frequency_eqn, IMP, 'Real', true); % Real
        ResZ = vpa(ResZ); % Numeric estimation for the symbolic solution
        if ~isempty(ResZ) % To check if the result is made or not
            Output.Capacitor(i) = ResZ(1);
        else
            Output.Capacitor(i) = nan;
        end
        Output.Rct(i) = Rct;
        Output.Rs(i) = Rs;
        
    else
        Phase_eqn = (Phase_Diff * a*c+b*d == b*c-a*d);
    
        ResZ = solve([Frequency_eqn, Phase_eqn], [IMP,Rct], 'Real', true);
        ResZIMP = vpa(ResZ.IMP);
        ResZRCT = vpa(ResZ.Rct);
        if ~isempty(ResZIMP) && ~isempty(ResZRCT) % To check if the result is made or not
            Output.Rct(i) = ResZRCT(end);
            Output.Capacitor(i) = ResZIMP(end);
            Output.Rs(i) = RT-Output.Rct(i);
        else
            Output.Rct(i) = nan;
            Output.Capacitor(i) = nan;
            Output.Rs(i) = nan;
        end
    end
    
    
    
    % Calculation for measured magnitude and phase of impedance
    Betha = atan(-Output.Vp_p2(i)/(Output.Vp_p1(i)-Output.Vp_p2(i)*cos(Phase_Diff)));
    Re_Voltage = sqrt(Output.Vp_p2(i)^2+Output.Vp_p1(i)^2-2*Output.Vp_p1(i)*Output.Vp_p2(i));
    Output.Measured_Imp(i) = Output.Vp_p2(i)/(Re_Voltage/R);
    Output.Measured_Phase(i) = Phase_Diff-Betha;
    
    
    % Calculation for modeled measured magnitude and phase of impedance
    e = Output.Rct(i) + Output.Rs(i) + Output.Rct(i)*Output.Rs(i)*Output.Capacitor(i)*(2*pi*f)^CPE*cos(CPE*pi/2);
    f = Output.Rct(i)*Output.Rs(i)*Output.Capacitor(i)*(2*pi*f)^CPE*sin(CPE*pi/2);
    g = Output.Rct(i)*Output.Capacitor(i)*(2*pi*f)^CPE*cos(CPE*pi/2)+1;
    h = Output.Rct(i)*Output.Capacitor(i)*(2*pi*f)^CPE*sin(CPE*pi/2);
    Output.Calculated_Impedance(i) = ((e*h+f*g)^2+(e*g-f*h))/(g^2+h^2);
    Output.Calculated_Phase(i) = atan((e*g-f*h)/(e*h+f*g));
    
    
    % Creating out put
    Output.GAIN(i) = Gain;
    Output.FREQUENCY(i) = Frequency;
    
    Output.CPE_Factor(i) = CPE;
end
