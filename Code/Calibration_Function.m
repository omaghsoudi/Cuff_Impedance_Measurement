function Output = Calibration_Function(max_sample_rate, ai1, ai2, ao)

Time = 1; % one second data acquiring
Amp = 5;

Devices = daq.getDevices;
Device_Name = Devices.ID;

%Create acquisition session
session = daq.createSession('ni');

sample_rate = max_sample_rate/5; % 3 channels are needed here
tot_samples = 1*sample_rate; % Total number of samples based on 1 second signal
%set acquisition parameters
session.NumberOfScans = tot_samples;
session.DurationInSeconds = 1;
session.Rate = sample_rate;
%set input channels
ch1 = addAnalogInputChannel(session,Device_Name, ai1, 'Voltage');
ch2 = addAnalogInputChannel(session,Device_Name, ai2, 'Voltage');
ch1.TerminalConfig = 'Differential';
ch2.TerminalConfig = 'Differential';
%set output channel
addAnalogOutputChannel(session,Device_Name,ao,'Voltage');


% DC Test
outputSignal1 = Amp*ones(tot_samples,1); % Time secs DC Signal
queueOutputData(session,repmat(outputSignal1,Time,1));
data_dc = session.startForeground; % Generate and acquire signal for 1 second
Output.Vdc_ch1(1) = mean(data_dc(:,1)); % Input voltage
Output.Vdc_ch2(1) = mean(data_dc(:,2)); % Output voltage
Output.Vdc(1) = Amp;




% DC zero level
outputSignal1 = zeros(tot_samples,1); 
queueOutputData(session,repmat(outputSignal1,Time,1));
data_dc = session.startForeground; % Generate and acquire signal for 1 second
Output.Vdc_ch1(2) = mean(data_dc(:,1)); % Input voltage
Output.Vdc_ch2(2) = mean(data_dc(:,2)); % Output voltage
Output.Vdc(2) = 0;
