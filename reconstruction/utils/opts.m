function out=opts(varargin)
%OPTS Set gradient limits of the MR system.
%   g=OPTS() Return the default amplitude and slew limits.
%
%   g=OPTS('maxGrad',30,'gradUnit','mT/m') Set the maximum gradient to
%   30mT/m.

persistent parser
validGradUnits={'Hz/m','mT/m','rad/ms/mm'};
validSlewUnits={'Hz/m/s','mT/m/ms','T/m/s','rad/ms/mm/ms'};
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'opts';
    parser.addParamValue('gradUnit',validGradUnits{1},...
        @(x) any(validatestring(x,validGradUnits)));
    parser.addParamValue('slewUnit',validSlewUnits{1},...
        @(x) any(validatestring(x,validSlewUnits)));
    parser.addParamValue('maxGrad',[],@isnumeric);
    parser.addParamValue('maxSlew',[],@isnumeric);
    parser.addParamValue('riseTime',[],@isnumeric);
    parser.addParamValue('rfDeadTime',0,@isnumeric);
    parser.addParamValue('rfRingdownTime',0,@isnumeric);
    parser.addParamValue('adcDeadTime',0,@isnumeric);
    parser.addParamValue('adcRasterTime',100e-9,@isnumeric);
    parser.addParamValue('rfRasterTime',1e-6,@isnumeric);
    parser.addParamValue('gradRasterTime',10e-6,@isnumeric);
    parser.addParamValue('blockDurationRaster',10e-6,@isnumeric);
    parser.addParamValue('gamma',42.576e6,@isnumeric); % Hz/T
    parser.addParamValue('B0',1.5,@isnumeric); % T
end
parse(parser,varargin{:});
opt = parser.Results;

if isempty(opt.maxGrad)
    maxGrad = mr.convert(40,'mT/m','gamma',opt.gamma);    % Default: 40 mT/m
else
    maxGrad = mr.convert(opt.maxGrad,opt.gradUnit,'Hz/m','gamma',opt.gamma);
end
if isempty(opt.maxSlew)
    maxSlew=mr.convert(170,'T/m/s','gamma',opt.gamma);	% Default: 170 mT/m/ms
else
    maxSlew = mr.convert(opt.maxSlew,opt.slewUnit,'Hz/m','gamma',opt.gamma);
end
if ~isempty(opt.riseTime)
    %maxSlew=[];
    maxSlew=maxGrad/opt.riseTime;
end

out.maxGrad = maxGrad;
out.maxSlew = maxSlew;
out.riseTime = opt.riseTime;
out.rfDeadTime = opt.rfDeadTime;
out.rfRingdownTime = opt.rfRingdownTime;
out.adcDeadTime = opt.adcDeadTime;
out.adcRasterTime = opt.adcRasterTime;
out.rfRasterTime = opt.rfRasterTime;
out.gradRasterTime = opt.gradRasterTime;
out.blockDurationRaster = opt.blockDurationRaster;
out.gamma=opt.gamma;
out.B0=opt.B0;

end