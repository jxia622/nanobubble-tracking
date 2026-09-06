function test_nbtracker()
% Behavioral tests, independent of the randomized benchmark and its scores.
root=fileparts(mfilename('fullpath'));code=fullfile(fileparts(root),'src');
addpath(code);
base={'UseEvidenceFilter',false,'LocalizationSigmaUm',1,'MinDetections',1,...
    'ConfirmationHits',1};
checks={};
[tr,adj,d]=nbtracker(cell(0,1),base{:});
assert(isempty(tr)&&isempty(adj)&&isempty(d.records));
[~,adj]=nbtracker(repmat({zeros(0,2)},12,1),base{:});assert(isempty(adj));
[tr,adj,d]=nbtracker({[12 15]},base{:});
assert(isequal(adj{1},1)&&isequal(tr{1},1)&&isequal(d.records.positionsPixels,[12 15]));
checks{end+1}='empty frames and singleton';

P=repmat({[40 60]},80,1);original=P;P(21:25)={zeros(0,2)};
[tr,adj,d]=nbtracker(P,base{:},'MaxMissedFrames',5);
assert(numel(adj)==1&&numel(adj{1})==75&&all(isnan(tr{1}(21:25))));
assert(isequal(d.records.frameIndices,find(~cellfun(@isempty,P))));
checkIndices(P,tr,adj,d);checks{end+1}='five-frame gap contains no synthetic observations';
P(26)={zeros(0,2)};[~,adj]=nbtracker(P,base{:},'MaxMissedFrames',5);
assert(numel(adj)==2);checks{end+1}='six missing frames expire identity';
P=original;P(21)={zeros(0,2)};[~,adj]=nbtracker(P,base{:},'MaxMissedFrames',0);
assert(numel(adj)==2);checks{end+1}='zero-gap boundary';
P={[40 60];[40 60]};[~,adj]=nbtracker(P,base{:},'FrameTimesSeconds',[0 7/2000]);
assert(numel(adj)==2);[~,adj]=nbtracker(P,base{:},'FrameTimesSeconds',[0 6/2000]);
assert(numel(adj)==1);checks{end+1}='elapsed-time expiry for irregular acquisition';

rng(3401);T=144;P=cell(T,1);ids=cell(T,1);
for f=1:T
    q=[50,20+.45*(f-1);50.2,80-.45*(f-1)]+.08*randn(2);
    order=randperm(2);P{f}=q(order,:);ids{f}=order(:);
end
before=P;[tr,adj,d]=nbtracker(P,'LocalizationSigmaUm',1,'UseEvidenceFilter',false);
labels=vertcat(ids{:});assert(numel(adj)==2);
for k=1:2, assert(max(histcounts(labels(adj{k}),[.5 1.5 2.5]))/numel(adj{k})>.98);end
checkIndices(P,tr,adj,d);assert(isequal(P,before));
shifted=cellfun(@(x)x+[1000 -300],P,'UniformOutput',false);
[~,shiftAdj]=nbtracker(shifted,'LocalizationSigmaUm',1,'UseEvidenceFilter',false);
assert(isequal(adj,shiftAdj));checks{end+1}='crossing identity, exclusive assignments, immutable input and translation invariance';

pixel=[10 40];t=(0:95)'/2000;t=t(mod((1:96)',7)~=0);v=[3000 -5000];
exact=[30 70]+t*v./pixel;observed=exact+.025*randn(size(exact));
P=mat2cell(observed,ones(numel(t),1),2);
[tr,adj,d]=nbtracker(P,'PixelSizeUm',pixel,'FrameTimesSeconds',t,...
    'LocalizationSigmaUm',1,'UseEvidenceFilter',false);
assert(numel(adj)==1&&numel(adj{1})==numel(t));
assert(max(abs(mean(d.records.velocityUmPerSec,1)-v)./abs(v))<.02);
assert(mean((d.records.smoothedPositionsPixels-exact).^2,'all')<mean((observed-exact).^2,'all'));
assert(isequal(d.records.timesSeconds,t));checkIndices(P,tr,adj,d);
checks{end+1}='physical velocity, anisotropic pixels, irregular timestamps and smoother accuracy';

mustFail(@()nbtracker({[NaN 1]}));mustFail(@()nbtracker({[1 2 3]}));
mustFail(@()nbtracker({[1 2];[1 2]},'FrameTimesSeconds',[1 1]));
mustFail(@()nbtracker({[1 2]},'FrameTimesSeconds',1i));
mustFail(@()nbtracker({[1 2]},'PixelSizeUm',[1 1i]));
checks{end+1}='invalid coordinates, geometry and timestamps rejected';
report=struct('passed',true,'checks',{checks},'matlab',version);
fid=fopen(fullfile(root,'unit_test_results.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
fprintf('PASS: %d behavioral test groups\n',numel(checks));
end

function checkIndices(P,tr,adj,d)
counts=cellfun(@(x)size(x,1),P);frames=repelem((1:numel(P))',counts);
offset=[0;cumsum(counts)];allPoints=vertcat(P{:});allIdx=vertcat(adj{:});
assert(numel(unique(allIdx))==numel(allIdx));
for k=1:numel(adj)
    idx=adj{k};f=frames(idx);assert(all(diff(f)>0));
    assert(isequal(d.records(k).positionsPixels,allPoints(idx,:)));
    assert(isequal(tr{k}(f),idx-offset(f)));
end
end

function mustFail(fn)
failed=false;try,fn();catch,failed=true;end
assert(failed,'Expected invalid input to be rejected');
end
