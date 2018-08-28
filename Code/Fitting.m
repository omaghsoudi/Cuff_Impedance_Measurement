function Output = Fitting(Output_Data, Inputs, mode)
% this is not working for fast method and electrode mode. Just working for cuff estimated and cuff two
% equation modes.
Output.Output = Output_Data.Output;


%Inputs
%[Phase error degree badsed range, Phase factor,                 Magnitude factor,                frequency range]
Frequencies = Output_Data.FREQUENCY';
Magnitudes = Output_Data.Measured_Imp';
Phases = Output_Data.Measured_Phase';

if Output_Data.Rs(1)<10
    Output_Data.Rs(1) = 11;
end
params = [Output_Data.Rs(1) Output_Data.Rct(1) 10e-6 0.7 Output_Data.Output.Input_Parameters.R];

LEN = length(Frequencies);
    
Output.FREQUENCY = nan(LEN,1);
Output.Rs = nan(LEN,1);
Output.Rct = nan(LEN,1);
Output.Capacitor = nan(LEN,1);
Output.Measured_Imp = nan(LEN,1);
Output.Measured_Phase = nan(LEN,1);
Output.Calculated_Impedance = nan(LEN,1);
Output.Calculated_Phase = nan(LEN,1);
Output.CPE_Factor = nan(LEN,1);

Steps = 100;
Range = linspace(0.01,1,50);
Range = cat(2, Range, linspace(1.01,10,50));

options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000, 'Display', 'off');
if length(Inputs) ==3
    for i = 1:LEN
        w = Frequencies(i);
        imd = Magnitudes(i);
        ipd = Phases(i);
        
        optcn = fminsearch( @(x) impError(x, imd, ipd, params, w, Inputs, mode), [params(1); params(2); params(3); params(4)], options);
        
        [fim, fip] = evalCuff([optcn(1) optcn(2) optcn(3) optcn(4) params(5)], w, mode);
        
        Output.FREQUENCY(i) = w;
        Output.Rs(i) = optcn(1);
        Output.Rct(i) = optcn(2);
        Output.Capacitor(i) = optcn(3);
        Output.CPE_Factor(i) = optcn(4);
        Output.Measured_Imp(i) = imd;
        Output.Measured_Phase(i) = ipd;
        Output.Calculated_Impedance(i) = fim;
        Output.Calculated_Phase(i) = fip;        
    end

elseif length(Inputs) ==6
    Temp.FREQUENCY = nan(LEN,Steps);
    Temp.Rs = nan(LEN,Steps);
    Temp.Rct = nan(LEN,Steps);
    Temp.Capacitor = nan(LEN,Steps);
    Temp.Measured_Imp = nan(LEN,Steps);
    Temp.Measured_Phase = nan(LEN,Steps);
    Temp.Calculated_Impedance = nan(LEN,Steps);
    Temp.Calculated_Phase = nan(LEN,Steps);
    Temp.CPE_Factor = nan(LEN,Steps);
    Temp.G_Error = nan(1,Steps);
    Temp.L_Error = nan(1,Steps);
    
    Index = find(Frequencies <= Inputs(5) & Frequencies >= Inputs(4));
    
    w = Frequencies(Index);

    imd = Magnitudes(Index);
    ipd = Phases(Index);
        
    for j = 1:Steps
        
        [optcn,locerror] = fminsearch( @(x) impError(x, imd, ipd, params, w, [Inputs(1) Range(j) 1], mode), [params(1); params(2); params(3); params(4)], options);
        
        [fim, fip] = evalCuff([optcn(1) optcn(2) optcn(3) optcn(4) params(5)], w, mode);
        
        Temp.FREQUENCY(:,j) = Frequencies(:);
        Temp.Rs(Index,j) = optcn(1);
        Temp.Rct(Index,j) = optcn(2);
        Temp.Capacitor(Index,j) = optcn(3);
        Temp.CPE_Factor(Index,j) = optcn(4);
        Temp.Measured_Imp(:,j) = Magnitudes(:);
        Temp.Measured_Phase(:,j) = Phases(:);
        Temp.Calculated_Impedance(Index,j) = fim;
        Temp.Calculated_Phase(Index,j) = fip;
        
        Temp.L_Error(1,j) = locerror;
        Temp.G_Error(1,j) = impError([optcn(1) optcn(2) optcn(3) optcn(4)], imd, ipd, params, w, Inputs, mode);
    end
    
    Loc = find(Temp.G_Error==min(Temp.G_Error)); Loc = Loc(1);
    
    Output.FREQUENCY(:) = Temp.FREQUENCY(:, Loc);
    Output.Rs(Index) = Temp.Rs(Index, Loc);
    Output.Rct(Index) = Temp.Rct(Index, Loc);
    Output.Capacitor(Index) = Temp.Capacitor(Index, Loc);
    Output.CPE_Factor(Index) = Temp.CPE_Factor(Index, Loc);
    Output.Measured_Imp(:) = Temp.Measured_Imp(:, Loc);
    Output.Measured_Phase(:) = Temp.Measured_Phase(:, Loc);
    Output.Calculated_Impedance(Index) = Temp.Calculated_Impedance(Index, Loc);
    Output.Calculated_Phase(Index) = Temp.Calculated_Phase(Index, Loc);
    Output.Output.G_Error = Temp.G_Error;
    Output.Output.L_Error = Temp.L_Error;
    
else
    Index = find(Frequencies <= Inputs(5) & Frequencies >= Inputs(4));
    
    w = Frequencies(Index);
    
    imd = Magnitudes(Index);
    ipd = Phases(Index);
    
    optcn = fminsearch( @(x) impError(x, imd, ipd, params, w, Inputs, mode), [params(1); params(2); params(3); params(4)], options);
    
    [fim, fip] = evalCuff([optcn(1) optcn(2) optcn(3) optcn(4) params(5)], w, mode);

    Output.FREQUENCY(:) = Frequencies(:);
    Output.Rs(Index) = optcn(1);
    Output.Rct(Index) = optcn(2);
    Output.Capacitor(Index) = optcn(3);
    Output.CPE_Factor(Index) = optcn(4);
    Output.Measured_Imp(:) = Magnitudes(:);
    Output.Measured_Phase(:) = Phases(:);
    Output.Calculated_Impedance(Index) = fim;
    Output.Calculated_Phase(Index) = fip;
end

end







function err = impError( x, imdata, ipdata, params, w, Inputs, mode )
[im, ip]=evalCuff([[x(1) x(2)] x(3) x(4) params(5)], w, mode);

% If values getting negative makes the error 100 to show that it is going
% wrong. Also the Rs cannot be less than 5 and CPE value and exponent cannot be higher
% than 1
if x(1)<5 || x(3)>1 || x(4)>1|| any(x<0)
    err = 1000;
else
    Factors = abs(ipdata)*180/pi;
    Factors(Factors<Inputs(1)) = Inputs(1);
    Factors = Factors*pi/180;
    A = (Inputs(3)*sum( abs(imdata - im)./abs(imdata)));
    B = (Inputs(2)*sum( abs(ipdata-ip)./Factors ));
    err =  sqrt((A+B)^2+(A-B)^2);
end

end


function [impmag,impphase] = evalCuff( params, w, mode )  % cuff model
Rs = params(1);
Rct = params(2);
C = params(3);
n = params(4);
Re = params(5);
    
if mode == 0
    % try impedance version:
    e = 2*Rct + 2*Rs +2*Rct*Rs*C*(w.^n)*cos(n*pi/2);
    f = 2*Rct*Rs*C*(w.^n)*sin(n*pi/2);
    g = Rct*C*(w.^n)*cos(n*pi/2) + 1;
    h = Rct*C*(w.^n)*sin(n*pi/2);
else
    e = Rct + Rs + Rct*Rs*C*(w.^n)*cos(n*pi/2);
    f = Rct*Rs*C*(w.^n)*sin(n*pi/2);
    g = Rct*C*(w.^n)*cos(n*pi/2) + 1;
    h = Rct*C*(w.^n)*sin(n*pi/2);
end

ReImp = (e.*g + f.*h) ./ (g.^2 + h.^2);
ImImp = (f.*g - e.*h) ./ (g.^2 + h.^2);

impmag = sqrt( ReImp.^2 + ImImp.^2 );
impphase = wrapToPi(atan( ImImp ./ ReImp ));

end
