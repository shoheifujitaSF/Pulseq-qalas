%--------------------------------------------------------------------------
% Define high-level parameters
%--------------------------------------------------------------------------

% Set system parameters
sys = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m', ...
    'MaxSlew', 160, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% Load ky-kz trajectory data
traj = readmatrix('csv_imports/QALAS_1x1x3mm_R1_TR58_ETL127_revised.txt');
traj = traj(:,[3,2,1]);
nTR = 58;
nETL = 127;

% Create sequence object and define FOV, resolution, and other high-level parameters
seq = mr.Sequence(sys);         % Create a new sequence object
fov=[192 168 168]*1e-3;         % Define FOV and resolution
N = [192 168 56];

Nx = N(1);
Ny = N(2);
Nz = N(3);

bandwidth = 340;                % Hz/pixel
flip_angle = 4;                 % degrees
rfSpoilingInc=117;              % RF spoiling increment
Tread = 1 / bandwidth;          % readout duration in sec
os_factor = 1;                  % readout oversampling amount

%--------------------------------------------------------------------------
% Prepare sequence blocks
%--------------------------------------------------------------------------

deltak = 1 ./ fov;
gx = mr.makeTrapezoid('x',sys,'FlatArea',-Nx*deltak(1),'FlatTime',ceil((Tread+sys.adcDeadTime)/sys.gradRasterTime)*sys.gradRasterTime,'system',sys);    % readout gradient
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2);                          % Gx prewinder
gxSpoil = mr.makeTrapezoid('x',sys,'Area',gx.area,'Duration',1.8e-3);         % Gx spoiler 1.8ms to match IDEA

% Make trapezoids for inner loop to save computation

gyPre = mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxPre), 'system',sys);
gzPre = mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxPre), 'system',sys);

gyReph = mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxSpoil), 'system',sys);
gzReph = mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxSpoil), 'system',sys);

stepsY=((traj(:,2)-1)-Ny/2)/Ny*2;
stepsZ=((traj(:,1)-1)-Nz/2)/Nz*2;

% Create non-selective pulse
rf = mr.makeBlockPulse(flip_angle*pi/180, sys, 'Duration', 0.1e-3);

% Spoilers after t2prep and IR prep
gslSp_t2prep = mr.makeTrapezoid('z','Amplitude',-42.58*8*1e3,'Risetime',0.85e-3,'Duration',0.85e-3+8e-3+0.85e-3,'system',sys); %Amplitude in Hz/m
gslSp_IRprep = mr.makeTrapezoid('z','Amplitude',-42.58*8*1e3,'Risetime',1e-3,'Duration',1e-3+8e-3+1e-3,'system',sys); %Amplitude in Hz/m

% Analog to digital convertersys
adc = mr.makeAdc(Nx * os_factor,'Dwell',floor((gx.flatTime)/N(1)/sys.adcRasterTime)*sys.adcRasterTime,'Delay',gx.riseTime);

% Prep pulses
% T2 prep and IR prep pulse are imported from external txt files
% txt file is in mag and phase, while mr.makeArbitraryRf assumes real and imag
t2prep = readmatrix('csv_imports/HyperbolicSecantRF_T2prep.txt');
Re = t2prep(:,1) .* cos(t2prep(:,2));
Im = t2prep(:,1) .* sin(t2prep(:,2));
t2prep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 380.4*pi/180, 'system',sys, 'dwell', 1e-6);

text = readmatrix('csv_imports/rf90.txt');
Re = text(:,1) .* cos(text(:,2));
Im = text(:,1) .* sin(text(:,2));
rf90 = mr.makeArbitraryRf((Re+Im*1i).', pi/2, 'system',sys, 'dwell', 1e-6);
rf90_180PhaseOffset = mr.makeArbitraryRf((Re+Im*1i).', pi/2, 'system',sys, 'dwell', 1e-6, 'PhaseOffset',-180*pi/180);

IRprep = readmatrix('csv_imports/mp2rage_invpulse_nograd.txt');
Re = IRprep(:,1) .* cos(IRprep(:,2));
Im = IRprep(:,1) .* sin(IRprep(:,2));
IRprep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 850*pi/180, 'system',sys, 'dwell', 1e-5);
%IRprep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 1500*pi/180, 'system',sys, 'dwell', 1e-5);

%-------------------------------------------------------------------------
% Hard-coded sequence timings to match IDEA implementation
%--------------------------------------------------------------------------

% SF: these  timings are hard coded to match IDEA implementation

esp = 5.8e-3;
te  = 2.29e-3;
gap_between_readouts = 900e-3;

delay_1_t2prep  =   11e-3 + 80e-6;      % 11.08 ms
delay_2_t2prep  =   25e-3;              % 25 ms
delay_3_t2prep  =   14e-3 - 80e-6;      % 13.92 ms
delay_IRprep    =   100e-3 - mr.calcDuration(IRprep_pulse)/2;         % gap between end of inversion and start of readout#2 

delay_TE        =   te - (mr.calcDuration(rf)/2 + mr.calcDuration(gxPre) + mr.calcDuration(gx)/2);     % Set this delay to 0 for minimum echo time achievement 
delay_TRinner   =   esp - (mr.calcDuration(rf) + delay_TE + mr.calcDuration(gxPre)+mr.calcDuration(gx)+mr.calcDuration(gxSpoil));         
delay_TRouter   =   gap_between_readouts - esp*nETL;
delT_M3_M4      =   gap_between_readouts - esp*nETL - mr.calcDuration(IRprep_pulse) - delay_IRprep;     % between end of readout#1 and start of inversion
delT_M3_M4      =   delT_M3_M4 - 0.22e-3;   % this 0.22 ms is needed to match IDEA
delT_M13_2end   =   53.5e-3;

%--------------------------------------------------------------------------
% Build sequence
%--------------------------------------------------------------------------

for iZ = 1:nTR
    rf_phase=0;
    rf_inc=0;

    % T2 prep pulse
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    % Contrast 1
    [rf_phase, rf_inc] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, te, gyPre,gyReph,gzPre,gzReph);
    seq.addBlock(mr.makeDelay(delT_M3_M4));

    % IR prep
    seq.addBlock(IRprep_pulse);
    seq.addBlock(gslSp_IRprep,mr.makeDelay(delay_IRprep));

    % Contrast 2
    [rf_phase, rf_inc] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc,stepsZ,stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, te, gyPre,gyReph,gzPre,gzReph);
    seq.addBlock(mr.makeDelay(delay_TRouter))

    % Contrast 3
    [rf_phase, rf_inc] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, te, gyPre,gyReph,gzPre,gzReph);
    seq.addBlock(mr.makeDelay(delay_TRouter));

    % Contrast 4
    [rf_phase, rf_inc] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, te, gyPre,gyReph,gzPre,gzReph);
    seq.addBlock(mr.makeDelay(delay_TRouter));

    % Contrast 5
    [rf_phase, rf_inc] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, te, gyPre,gyReph,gzPre,gzReph);
    seq.addBlock(mr.makeDelay(delT_M13_2end));

end

%--------------------------------------------------------------------------
% Check timing and write sequence
%--------------------------------------------------------------------------

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('/autofs/cluster/berkin/fujita/matlab_codes/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma

if (pns_ok)
     fprintf('PNS check passed successfully\n');
else
     fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end

gradSpectrum

% Set definitions
seq.setDefinition('FOV', fov);

% plot
seq.plot('TimeRange',[0 1.5],'timeDisp','ms');

% Write to pulseq file
seq.write(['qalas_1x1x3mm.seq']);

%--------------------------------------------------------------------------
% Functions
%--------------------------------------------------------------------------

function [rf_phase, rf_inc] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner,te, gyPre,gyReph,gzPre,gzReph)
for iY = 1:nETL

    % Calculate index for ky-kz look up table
    index = iY + (iZ-1)*nETL;

    % RF spoiling
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);       %increment RF phase

    % Excitation
    seq.addBlock(rf);

    % Delay included to attain the desired echo time
    seq.addBlock(mr.makeDelay(delay_TE));

    % Encoding
    seq.addBlock(gxPre, ...
        mr.scaleGrad(gyPre,-stepsY(index)), ...
        mr.scaleGrad(gzPre,-stepsZ(index)));    % Gz, Gy blips, Gx pre-winder

    seq.addBlock(gx, adc);                      % Gx readout

    seq.addBlock(gxSpoil, ...
        mr.scaleGrad(gyReph,stepsY(index)), ...
        mr.scaleGrad(gzReph,stepsZ(index)));    % -Gz, -Gy blips, Gx spoiler

    seq.addBlock(mr.makeDelay(delay_TRinner));  % wait until desired echo spacing
end
end
