function run_holdout()
root=fileparts(mfilename('fullpath')); repo=fileparts(root);
addpath(fullfile(repo,'src'));
dest=fullfile(root,'holdout_results');if ~isfolder(dest),mkdir(dest);end
jobs=jsondecode(fileread(fullfile(root,'holdout_manifest.json')));
for c=1:numel(jobs)
    job=jobs(c);out=fullfile(dest,[job.key '.mat']);if isfile(out),continue;end
    S=load(fullfile(root,'holdout_inputs',[job.key '.mat']));P=S.SR_Localizations;
    result=struct();
    for g=[1 6]
        tic;[~,adj]=simpletracker(P,'MaxLinkingDistance',8,'MaxGapClosing',g,'Debug',false);
        adj=adj(cellfun(@numel,adj)>15);
        name='legacy_default';if g==6,name='legacy_gap6';end
        result.(name)=struct('tracks',{adj},'seconds',toc);
    end
    for mode=1:3
        name='nb';extra={};
        if mode==1,name='nb_position_only';extra={'UseMotionPrediction',false};end
        if mode==2,name='nb_no_evidence';extra={'UseEvidenceFilter',false};end
        tic;[~,adj,details]=nbtracker(P,'FrameTimesSeconds',S.frame_times,...
            'PixelSizeUm',S.scale_um,extra{:});
        result.(name)=struct('tracks',{adj},'seconds',toc,'details',details);
    end
    save(out,'result','-v7');
    fprintf('%d/%d %s | old %.2fs NB %.2fs\n',c,numel(jobs),job.key,...
        result.legacy_default.seconds,result.nb.seconds);
end
end
