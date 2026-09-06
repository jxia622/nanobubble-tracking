function [tracks, adjacency_tracks, details] = nbtracker(points, varargin)
%NBTRACKER Motion-predicted assignment for intermittent ultrasound detections.
% Input: cell vector of N-by-2 [z,x] observed localization coordinates, pixels.
% Output: same track/global-index convention as simpletracker. Only observed
% detections are returned; predictions through gaps are never counted as data.
%
% [tracks, adjacency_tracks, details] = nbtracker(points, ...
%   'FrameRateHz',2000, 'PixelSizeUm',[12.32 12.4145], ...
%   'FrameTimesSeconds',times, 'MinDetections',16);
%
% State prediction is part of assignment, unlike post-hoc Kalman smoothing.
% Missed observations preserve the state and grow its uncertainty. Exclusive
% assignments minimize Gaussian innovation cost including log-determinant,
% with explicit unmatched alternatives. Connected components bound Hungarian
% work. A small confirmation window suppresses transient false tracks.
% No localization coordinates, SVD cutoffs, amplitudes or PSFs are changed.
% Requires munkres.m in this source directory (no Tracking Toolbox needed).
% LocalizationSigmaUm=[] estimates uncertainty from nearest displacements;
% supply a measured standard deviation when available (see estimateSigma).
% Evidence is a spatial likelihood score, not a calibrated probability.

p=inputParser;
p.addParameter('FrameRateHz',2000,@positiveScalar);
p.addParameter('PixelSizeUm',[12.32 12.4145],@(x)isnumeric(x)&&isreal(x)&&numel(x)==2&&all(isfinite(x(:)))&&all(x(:)>0));
p.addParameter('FrameTimesSeconds',[],@(x)isnumeric(x)&&isreal(x)&&isvectorOrEmpty(x));
p.addParameter('LocalizationSigmaUm',[],@(x)isempty(x)||positiveScalar(x));
p.addParameter('AccelerationStdUmPerSec2',5e5,@positiveScalar);
p.addParameter('InitialVelocityStdUmPerSec',2e4,@positiveScalar);
p.addParameter('MaxLinkingDistance',8,@positiveScalar);
p.addParameter('MaxMissedFrames',5,@nonnegativeInteger);
p.addParameter('MinDetections',16,@positiveInteger);
p.addParameter('ConfirmationHits',3,@positiveInteger);
p.addParameter('ConfirmationWindow',5,@positiveInteger);
p.addParameter('GateChiSquare',13.8155,@positiveScalar); % 99.9%, 2 dimensions
p.addParameter('UnmatchedCost',6,@positiveScalar);
p.addParameter('DetectionProbability',.8,@(x)positiveScalar(x)&&x<1);
p.addParameter('MinLogEvidence',log(100),@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x));
p.addParameter('UseEvidenceFilter',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('UseMotionPrediction',true,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opt=p.Results;
assert(exist('munkres','file')==2,'nbtracker:MissingDependency',...
    'Add this repository''s src folder to the MATLAB path (munkres.m).');
assert(iscell(points)&&isvectorOrEmpty(points),'nbtracker:InvalidPoints','points must be a cell vector.');
assert(opt.ConfirmationHits<=opt.ConfirmationWindow,'nbtracker:InvalidConfirmation','ConfirmationHits exceeds window.');
points=points(:); T=numel(points);
for f=1:T
    if isempty(points{f}), points{f}=zeros(0,2); end
    assert(isnumeric(points{f})&&isreal(points{f})&&size(points{f},2)==2&&all(isfinite(points{f}(:))),...
        'nbtracker:InvalidPoints','Each frame must contain finite real N-by-2 coordinates.');
    points{f}=double(points{f});
end
times=opt.FrameTimesSeconds(:);
if isempty(times), times=(0:T-1)'/opt.FrameRateHz; end
assert(numel(times)==T&&all(isfinite(times))&&all(diff(times)>0),...
    'nbtracker:InvalidTimes','FrameTimesSeconds must be finite, strictly increasing and match frame count.');
pixel=double(opt.PixelSizeUm(:)'); unit=mean(pixel); dt0=1/opt.FrameRateHz;
if isempty(opt.LocalizationSigmaUm)
    opt.LocalizationSigmaUm=estimateSigma(points,pixel,unit);
end
sigma=opt.LocalizationSigmaUm/unit; R=eye(2)*sigma^2;
accel=opt.AccelerationStdUmPerSec2*dt0^2/unit;
initialVelocity=opt.InitialVelocityStdUmPerSec*dt0/unit;
H=[eye(2) zeros(2)]; I=eye(4);
maxAge=(opt.MaxMissedFrames+1)*dt0;
count=cellfun(@(x)size(x,1),points); offset=[0;cumsum(count)];
allpoints=vertcat(points{:});
if isempty(allpoints), area=1; else
    extent=(max(allpoints,[],1)-min(allpoints,[],1)).*pixel/unit;
    area=prod(max(extent,8));
end
template=struct('state',zeros(4,1),'P',zeros(4),'lastTime',0,'lastSeen',0,...
    'indices',zeros(0,1),'frames',zeros(0,1),'history',zeros(0,1),'confirmed',false,'logEvidence',0);
active=repmat(template,0,1); finished=repmat(template,0,1);
stats=struct('assignments',0,'births',0,'expired',0,'maxActiveTracks',0,'componentMaxSize',0);

for f=1:T
    % Expire by physical elapsed time before association, including irregular
    % acquisition timestamps. The boundary admits exactly MaxMissedFrames.
    expired=arrayfun(@(s)times(f)-s.lastSeen>maxAge+dt0*1e-8,active);
    finished=[finished;active(expired)]; %#ok<AGROW>
    stats.expired=stats.expired+sum(expired); active=active(~expired);
    n=numel(active); obs=points{f}.*pixel/unit; m=size(obs,1);
    costs=Inf(n,m);
    evidence=-Inf(n,m);
    for k=1:n
        dt=(times(f)-active(k).lastTime)/dt0;
        [F,Q]=motionModel(dt,accel);
        if ~opt.UseMotionPrediction, active(k).state(3:4)=0; end
        active(k).state=F*active(k).state;
        active(k).P=F*active(k).P*F'+Q;
        active(k).lastTime=times(f);
        residual=obs-active(k).state(1:2)';
        S=active(k).P(1:2,1:2)+R;
        % Covariance has independent axes but use a full solve for stability
        % if the state model is later extended to include correlated errors.
        d2=sum((residual/S).*residual,2);
        pixelDistance=sqrt(sum((residual*unit./pixel).^2,2));
        valid=d2<=opt.GateChiSquare & pixelDistance<=opt.MaxLinkingDistance;
        logdet=2*sum(log(diag(chol(S))))-2*sum(log(diag(chol(R))));
        costs(k,valid)=.5*(d2(valid)+logdet);
        % Conservative spatial null: all current detections contribute to
        % clutter density, including real bubbles. This likelihood ratio
        % rewards coherent hits and charges missed observations; a long
        % chance chain is not accepted merely for reaching MinDetections.
        density=max(m,1)/area;
        evidence(k,valid)=log(opt.DetectionProbability/(2*pi*sqrt(det(S))*density))-.5*d2(valid);
    end
    [pairs,componentSize]=assignComponents(costs,opt.UnmatchedCost);
    stats.componentMaxSize=max(stats.componentMaxSize,componentSize);
    assignedTrack=false(n,1); assignedDetection=false(m,1);
    for a=1:size(pairs,1)
        k=pairs(a,1); d=pairs(a,2);
        K=active(k).P*H'/(H*active(k).P*H'+R);
        active(k).state=active(k).state+K*(obs(d,:)'-H*active(k).state);
        % Joseph covariance update keeps uncertainty positive semidefinite.
        active(k).P=(I-K*H)*active(k).P*(I-K*H)'+K*R*K';
        active(k).P=(active(k).P+active(k).P')/2;
        active(k).lastSeen=times(f);
        active(k).logEvidence=active(k).logEvidence+evidence(k,d);
        active(k).indices(end+1,1)=offset(f)+d;
        active(k).frames(end+1,1)=f;
        assignedTrack(k)=true; assignedDetection(d)=true;
    end
    for k=1:n
        if ~assignedTrack(k), active(k).logEvidence=active(k).logEvidence+log1p(-opt.DetectionProbability); end
        active(k).history=[active(k).history;double(assignedTrack(k))];
        active(k).history=active(k).history(max(1,end-opt.ConfirmationWindow+1):end);
        if sum(active(k).history)>=opt.ConfirmationHits, active(k).confirmed=true; end
    end
    % Tentative tracks get the same physical gap allowance, but must gather
    % enough evidence within their first confirmation window to stay alive.
    rejected=arrayfun(@(s)~s.confirmed&&numel(s.history)>=opt.ConfirmationWindow,active);
    finished=[finished;active(rejected)]; %#ok<AGROW>
    active=active(~rejected);
    for d=find(~assignedDetection)'
        s=template; s.state=[obs(d,:)';0;0];
        s.P=diag([sigma^2 sigma^2 initialVelocity^2 initialVelocity^2]);
        s.lastTime=times(f); s.lastSeen=times(f);
        s.indices=offset(f)+d; s.frames=f; s.history=1;
        s.confirmed=opt.ConfirmationHits==1;
        active(end+1,1)=s; %#ok<AGROW>
    end
    stats.assignments=stats.assignments+size(pairs,1);
    stats.births=stats.births+sum(~assignedDetection);
    stats.maxActiveTracks=max(stats.maxActiveTracks,numel(active));
end
finished=[finished;active];
keep=arrayfun(@(s)s.confirmed&&numel(s.indices)>=opt.MinDetections&&...
    (~opt.UseEvidenceFilter||s.logEvidence>=opt.MinLogEvidence),finished);
finished=finished(keep);
if ~isempty(finished)
    [~,order]=sort(arrayfun(@(s)s.indices(1),finished)); finished=finished(order);
end
tracks=cell(numel(finished),1); adjacency_tracks=cell(numel(finished),1);
recordTemplate=struct('id',0,'globalIndices',[],'frameIndices',[],'timesSeconds',[],...
    'positionsPixels',[],'smoothedPositionsPixels',[],'velocityUmPerSec',[],...
    'positionStdUm',[],'maxMissingFrames',0,'logEvidence',0);
records=repmat(recordTemplate,numel(finished),1);
for k=1:numel(finished)
    s=finished(k); adjacency_tracks{k}=s.indices;
    tr=NaN(T,1); tr(s.frames)=s.indices-offset(s.frames); tracks{k}=tr;
    pos=allpoints(s.indices,:); tt=times(s.frames);
    [smooth,vel,stdev]=smoothTrack(pos,tt,pixel,unit,dt0,sigma,accel,initialVelocity);
    records(k)=struct('id',k,'globalIndices',s.indices,'frameIndices',s.frames,...
        'timesSeconds',tt,'positionsPixels',pos,'smoothedPositionsPixels',smooth,...
        'velocityUmPerSec',vel,'positionStdUm',stdev,'maxMissingFrames',max([0;diff(s.frames)-1]),...
        'logEvidence',s.logEvidence);
end
details=struct('algorithm','NB motion-predicted assignment v1','options',opt,...
    'frameTimesSeconds',times,'records',records,'statistics',stats);
end

function [F,Q]=motionModel(dt,accel)
F=[eye(2) dt*eye(2);zeros(2) eye(2)];
% Continuous white acceleration, dt in nominal frame intervals. Unlike a
% per-observation noise term, this composes consistently through empty frames.
Q=accel^2*[dt^3/3*eye(2) dt^2/2*eye(2);dt^2/2*eye(2) dt*eye(2)];
end

function sigma=estimateSigma(points,pixel,unit)
% Robust initial uncertainty estimate from observed displacement magnitudes.
% Independent isotropic localization errors have nearest-pair radial median
% 2*sqrt(log(2))*sigma in a sparse, slowly moving scene. Missing detections and
% motion bias upward; competing dense detections bias downward. An explicit
% measured LocalizationSigmaUm overrides this fallback. No ground truth used.
distances=zeros(0,1);
for f=1:min(numel(points)-1,32)
    a=points{f}.*pixel; b=points{f+1}.*pixel;
    if isempty(a)||isempty(b), continue; end
    % Bound calibration work on dense frames, using deterministic subsampling.
    take=unique(round(linspace(1,size(a,1),min(128,size(a,1)))));
    for k=take
        dd=sum((b-a(k,:)).^2,2);
        distances(end+1,1)=sqrt(min(dd)); %#ok<AGROW>
    end
end
if isempty(distances), sigma=unit; else
    sigma=max(.05*unit,median(distances)/(2*sqrt(log(2))));
end
end

function [pairs,maxSize]=assignComponents(costs,unmatched)
% Solve each connected gated bipartite component independently. An explicit
% unmatched alternative prevents forced associations inside a loose gate.
[n,m]=size(costs); valid=isfinite(costs)&costs<2*unmatched;
remaining=any(valid,2); pairs=zeros(0,2); maxSize=0;
while any(remaining)
    rr=find(remaining,1); cc=find(any(valid(rr,:),1));
    while true
        nextRows=find(any(valid(:,cc),2));
        nextCols=find(any(valid(nextRows,:),1));
        if numel(nextRows)==numel(rr)&&numel(nextCols)==numel(cc), break; end
        rr=nextRows; cc=nextCols;
    end
    remaining(rr)=false; a=numel(rr); b=numel(cc); maxSize=max(maxSize,a+b);
    block=costs(rr,cc); block(~valid(rr,cc))=Inf;
    augmented=Inf(a+b,a+b);
    augmented(1:a,1:b)=block;
    augmented(sub2ind(size(augmented),(1:a)',b+(1:a)'))=unmatched;
    augmented(sub2ind(size(augmented),a+(1:b)',(1:b)'))=unmatched;
    augmented(a+1:end,b+1:end)=1e-12;
    assignment=munkres(augmented);
    for k=1:a
        target=assignment(k);
        if target>=1&&target<=b&&valid(rr(k),cc(target))
            pairs(end+1,:)=[rr(k),cc(target)]; %#ok<AGROW>
        end
    end
end
end

function [position,velocity,stdev]=smoothTrack(pos,t,pixel,unit,dt0,sigma,accel,initialVelocity)
% Rauch-Tung-Striebel smoothing at observed timestamps only. Prediction spans
% the true elapsed time across gaps; the first observation is assimilated once.
L=size(pos,1); z=pos.*pixel/unit; R=eye(2)*sigma^2; H=[eye(2) zeros(2)];
xf=zeros(4,L); Pf=zeros(4,4,L); xp=xf; Pp=Pf;
xf(:,1)=[z(1,:)';0;0]; Pf(:,:,1)=diag([sigma^2 sigma^2 initialVelocity^2 initialVelocity^2]);
for k=2:L
    [F,Q]=motionModel((t(k)-t(k-1))/dt0,accel);
    xp(:,k)=F*xf(:,k-1); Pp(:,:,k)=F*Pf(:,:,k-1)*F'+Q;
    K=Pp(:,:,k)*H'/(H*Pp(:,:,k)*H'+R);
    xf(:,k)=xp(:,k)+K*(z(k,:)'-H*xp(:,k));
    Pf(:,:,k)=(eye(4)-K*H)*Pp(:,:,k)*(eye(4)-K*H)'+K*R*K';
end
xs=xf; Ps=Pf;
for k=L-1:-1:1
    [F,~]=motionModel((t(k+1)-t(k))/dt0,accel);
    C=Pf(:,:,k)*F'/Pp(:,:,k+1);
    xs(:,k)=xf(:,k)+C*(xs(:,k+1)-xp(:,k+1));
    Ps(:,:,k)=Pf(:,:,k)+C*(Ps(:,:,k+1)-Pp(:,:,k+1))*C';
end
position=xs(1:2,:)'*unit./pixel; velocity=xs(3:4,:)'*unit/dt0;
stdev=zeros(L,2);
for k=1:L, stdev(k,:)=sqrt(max(0,diag(Ps(1:2,1:2,k))))'*unit; end
end

function b=positiveScalar(x), b=isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>0; end
function b=positiveInteger(x), b=positiveScalar(x)&&x==floor(x); end
function b=nonnegativeInteger(x), b=isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x)&&x>=0&&x==floor(x); end
function b=isvectorOrEmpty(x), b=isvector(x)||isempty(x); end
