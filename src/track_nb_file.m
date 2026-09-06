function output = track_nb_file(inputMAT, varargin)
%TRACK_NB_FILE Track a bounded window of existing localization lists.
% Does not run SVD or localization. By default processes at most 192 frames.
% Required geometry: PixelSizeUm=[axial,lateral] in the localization grid.
% For this L22 2 kHz, 128-to-1024 endpoint-interpolated grid:
% result=track_nb_file('localizations.mat','PixelSizeUm',[12.32 12.414467],...
%     'FrameRateHz',2000,'MaxFrames',192,'OutputMAT','new_nb_tracks.mat');
% Expected SR_Localizations: cell vector of N-by-2 [z,x] pixel coordinates.
% Uses meta.dsRange (original one-based frame indices) or frame_times when
% present. Otherwise assumes contiguous frames at FrameRateHz. Explicit
% FrameTimesSeconds overrides metadata and must cover the full input list.
% TrackerOptions is a cell of additional nbtracker name/value pairs.
% Output global indices refer to output.SR_Localizations (selected window).
% Output files must be new; source and previous results are not overwritten.

p=inputParser;
p.addParameter('PixelSizeUm',[],@(x)isnumeric(x)&&isreal(x)&&numel(x)==2&&all(isfinite(x(:)))&&all(x(:)>0));
p.addParameter('FrameRateHz',2000,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&isreal(x)&&x>0);
p.addParameter('StartFrame',1,@positiveInteger);
p.addParameter('MaxFrames',192,@positiveInteger);
p.addParameter('FrameTimesSeconds',[],@(x)isnumeric(x)&&isreal(x)&&(isvector(x)||isempty(x)));
p.addParameter('TrackerOptions',{},@(x)iscell(x)&&mod(numel(x),2)==0);
p.addParameter('OutputMAT','',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
p.parse(varargin{:});opt=p.Results;
assert(~isempty(opt.PixelSizeUm),'track_nb_file:GeometryRequired',...
    'Supply PixelSizeUm for the localization grid; native and interpolated lateral spacing differ.');
if strlength(string(opt.OutputMAT))>0
    assert(~isfile(opt.OutputMAT),'track_nb_file:ExistingOutput','OutputMAT already exists; choose a new filename.');
end
% Load only saved localization variables, never RF or IQ arrays.
info=whos('-file',inputMAT);names={info.name};
assert(any(strcmp(names,'SR_Localizations')),'track_nb_file:MissingLocalizations',...
    'Input must contain the saved SR_Localizations cell vector.');
wanted=intersect({'SR_Localizations','meta','frame_times'},names,'stable');
S=load(inputMAT,wanted{:});
assert(iscell(S.SR_Localizations)&&(isvector(S.SR_Localizations)||isempty(S.SR_Localizations)),'track_nb_file:MissingLocalizations',...
    'Input must contain the saved SR_Localizations cell vector.');
P=S.SR_Localizations(:);T=numel(P);
assert(opt.StartFrame<=T,'track_nb_file:FrameRange','StartFrame exceeds available frames.');
selected=(opt.StartFrame:min(T,opt.StartFrame+opt.MaxFrames-1))';
times=opt.FrameTimesSeconds(:);timeSource='explicit FrameTimesSeconds';
if isempty(times)
    if isfield(S,'frame_times')
        times=S.frame_times(:);timeSource='saved frame_times';
    elseif isfield(S,'meta')&&isfield(S.meta,'dsRange')
        times=(double(S.meta.dsRange(:))-1)/opt.FrameRateHz;timeSource='meta.dsRange / acquisition FrameRateHz';
    else
        times=(0:T-1)'/opt.FrameRateHz;timeSource='assumed contiguous frames at FrameRateHz';
    end
end
assert(numel(times)==T&&isreal(times)&&all(isfinite(times))&&all(diff(times)>0),...
    'track_nb_file:InvalidTimes','Timestamps must be increasing and cover all input frames.');
reserved={'PixelSizeUm','FrameRateHz','FrameTimesSeconds'};
for k=1:2:numel(opt.TrackerOptions)
    key=char(opt.TrackerOptions{k});
    assert(~any(strcmpi(key,reserved)),'track_nb_file:ReservedOption',...
        'Pass geometry and timestamps as wrapper options, not inside TrackerOptions.');
end
code=fileparts(mfilename('fullpath'));addpath(code);
tic;[tracks,adjacency_tracks,details]=nbtracker(P(selected),...
    'PixelSizeUm',opt.PixelSizeUm,'FrameRateHz',opt.FrameRateHz,...
    'FrameTimesSeconds',times(selected),opt.TrackerOptions{:});elapsed=toc;
sourceMeta=struct();if isfield(S,'meta'),sourceMeta=S.meta;end
output=struct('SR_Localizations',{P(selected)},'tracks',{tracks},...
    'adjacency_tracks',{adjacency_tracks},'details',details,...
    'sourceFile',char(inputMAT),'sourceMeta',sourceMeta,...
    'sourceFrameIndices',selected,'timeSource',timeSource,'seconds',elapsed);
if strlength(string(opt.OutputMAT))>0,save(opt.OutputMAT,'-struct','output','-v7.3');end
fprintf('NB tracker: %d saved frames, %d retained tracks, %.2f s (%s).\n',...
    numel(selected),numel(tracks),elapsed,timeSource);
end

function b=positiveInteger(x)
b=isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>0&&x==floor(x);
end
