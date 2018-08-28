function Output = PostProcessing(Output_Data, GainFactor, CPE, mode, R)
% this is not working for fast method and electrode mode. Just working for cuff estimated and cuff two
% equation modes.


% Get a copy of inputs for output saving
Output.Input_Parameters.max_sample_rate = Output_Data.Output.Input_Parameters.max_sample_rate;
Output.Input_Parameters.MaxFre = Output_Data.Output.Input_Parameters.MaxFre;
Output.Input_Parameters.MinFre = Output_Data.Output.Input_Parameters.MinFre;
Output.Input_Parameters.LBS = Output_Data.Output.Input_Parameters.LBS;
Output.Input_Parameters.ai1 = Output_Data.Output.Input_Parameters.ai1;
Output.Input_Parameters.ai2 = Output_Data.Output.Input_Parameters.ai2;
Output.Input_Parameters.ao = Output_Data.Output.Input_Parameters.ao;
Output.Input_Parameters.Time = Output_Data.Output.Input_Parameters.Time;
Output.Input_Parameters.Amp = Output_Data.Output.Input_Parameters.Amp;
Output.Input_Parameters.R = R;
Output.Input_Parameters.GainFactor = GainFactor;
Output.Input_Parameters.CPE = CPE;


sample_rate = Output.Input_Parameters.max_sample_rate/5; % 3 channels are needed here

MaxFre = Output_Data.Output.Input_Parameters.MaxFre;
MinFre = Output_Data.Output.Input_Parameters.MinFre;
LBS = Output_Data.Output.Input_Parameters.LBS;


Output.Vdc_in = Output_Data.Output.Vdc_in;
Output.Vdc_Cuff = Output_Data.Output.Vdc_Cuff;


if mode == 0
    % High Frequency test
    % Ideal case is having a device with a 1G samples per second and then the frequency of
    % this signal will be around 3MHz and then Rs and Rct will be calculated more accuarte 
    % but a total sample rate of 1,500,000 for all three channels is also acceptable
    % for just estimating capacitance, Rct, and Rs. It is enough to say that
    % cuff works properly or not.
    Output.Vac_in = Output_Data.Output.Vac_in;
    Output.Vac_Cuff = Output_Data.Output.Vac_Cuff;

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
if mode == 0 || mode == 2
    Rct = RT-Rs; % Rs is small and Rct is almost 1000 times larger; therefore, an estimation for Rs is not affecting Rct 
    Rct = Rct/2; % make it half as two are in series
    Output.Rct = Rct;
else
    %Accurate
    syms Rct positive
    Rs = (RT-2*Rct)/2;
end
    

% Swipping frequencies to find Capacitance
Output.Frequencies = logspace(log10(MinFre), log10(MaxFre), LBS); % Frequencies to test the cuff


LEN = length(Output.Frequencies);
% Initialize the output matrix
Output.GAIN = zeros(LEN,1);
Output.FREQUENCY = zeros(LEN,1);
Output.Rs = zeros(LEN,1);
Output.Rct = zeros(LEN,1);
Output.Capacitor = zeros(LEN,1);
Output.Vp_p1 = zeros(LEN,1);
Output.Vp_p2 = zeros(LEN,1);
Output.Measured_Imp = zeros(LEN,1);
Output.Measured_Phase = zeros(LEN,1);
Output.Calculated_Impedance = zeros(LEN,1);
Output.Calculated_Phase = zeros(LEN,1);
Output.CPE_Factor = zeros(LEN,1);


syms IMP positive % just positive response for IMP (capacitance) is acceptable
for i = 1: length(Output.Frequencies)
    % Data acquisition for ac signals
    Frequency = Output.Frequencies(i);
    
    Output.DATA(i).data1 = Output_Data.Output.DATA(i).data1;
    Value1 = findpeaks(Output.DATA(i).data1); % finding the positive peaks in input signal
    Value1 = mean(Value1);
    Value2 = findpeaks((-1)*Output.DATA(i).data1); % finding the negative peaks in input signal
    Value2 = mean(Value2);
    Output.Vp_p1(i) = abs(Value1)+abs(Value2); % Volatge peak to peak calculation
    Output.DATA(i).DATA1 = diff(sign(Output.DATA(i).data1)); % for zero crossing

    Output.DATA(i).data2 = Output_Data.Output.DATA(i).data2;
    Value1 = findpeaks(Output.DATA(i).data2); % finding the positive peaks in cuff signal
    Value1 = mean(Value1);
    Value2 = findpeaks((-1)*Output.DATA(i).data2); % finding the negative peaks in cuff signal
    Value2 = mean(Value2);
    Output.Vp_p2(i) = abs(Value1)+abs(Value2); % Volatge peak to peak calculation
    Output.DATA(i).DATA2 = diff(sign(Output.DATA(i).data2));
    
    [acor,lag] = xcorr(Output.DATA(i).data2, Output.DATA(i).data1);
    [~,I] = max(abs(acor));
    lagDiff = lag(I);
    
    % finding Capacitor by estimation
    Gain = Output.Vp_p2(i)/Output.Vp_p1(i); % calculating gain for each frequency
    f = Frequency;
    
    % it is not being used now, cross correlation has been using instead
%     % finding Capacitor and RCT more accurate
%     % find zero crossing
%     indx1 = find(Output.DATA(i).DATA1>0);
%     indx2 = find(Output.DATA(i).DATA2>0);
%     minlength = min(length(indx1), length(indx2));
%     Sample_Diff = mean(indx2(1:minlength)-indx1(1:minlength));
    
    Time_Diff = lagDiff/sample_rate;
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
    
    if mode == 0 || mode == 2
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
    Betha = atan(-Output.Vp_p2(i)*sin(Phase_Diff)/(Output.Vp_p1(i)-Output.Vp_p2(i)*cos(Phase_Diff)));
    Re_Voltage = sqrt(Output.Vp_p2(i)^2+Output.Vp_p1(i)^2-2*Output.Vp_p1(i)*Output.Vp_p2(i)*cos(Phase_Diff));
    Output.Measured_Imp(i) = Output.Vp_p2(i)/(Re_Voltage/R);
    Output.Measured_Phase(i) = -wrapToPi(Phase_Diff-Betha);
    
    
    % Calculation for modeled measured magnitude and phase of impedance
    E = 2*Output.Rct(i) + 2*Output.Rs(i) + 2*Output.Rct(i)*Output.Rs(i)*Output.Capacitor(i)*(2*pi*f)^CPE*cos(CPE*pi/2);
    F = 2*Output.Rct(i)*Output.Rs(i)*Output.Capacitor(i)*(2*pi*f)^CPE*sin(CPE*pi/2);
    G = Output.Rct(i)*Output.Capacitor(i)*(2*pi*f)^CPE*cos(CPE*pi/2)+1;
    H = Output.Rct(i)*Output.Capacitor(i)*(2*pi*f)^CPE*sin(CPE*pi/2);
    Output.Calculated_Impedance(i) = sqrt((E*G+F*H)^2+(F*G-E*H)^2)/(G^2+H^2);
    Output.Calculated_Phase(i) = wrapToPi(atan((F*G-E*H)/(E*G+F*H)));
    
    
    % Creating out put
    Output.GAIN(i) = Gain;
    Output.FREQUENCY(i) = Frequency;
    
    Output.CPE_Factor(i) = CPE;
end