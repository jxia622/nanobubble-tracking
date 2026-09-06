function test_file_workflow()
root=fileparts(mfilename('fullpath'));code=fullfile(fileparts(root),'src');addpath(code);
work=tempname;mkdir(work);cleanup=onCleanup(@()rmdir(work,'s')); %#ok<NASGU>
SR_Localizations=cell(96,1);for f=1:96,SR_Localizations{f}=[40 40+.1*f];end
meta=struct('dsRange',101:196);input=fullfile(work,'input.mat');dest=fullfile(work,'output.mat');
save(input,'SR_Localizations','meta','-v7.3');
out=track_nb_file(input,'PixelSizeUm',[12.32 12.414467],...
    'StartFrame',5,'MaxFrames',64,'OutputMAT',dest);
assert(isequal(out.SR_Localizations,SR_Localizations(5:68)));
assert(isequal(out.sourceFrameIndices,(5:68)'));
assert(isequal(out.details.frameTimesSeconds,(104:167)'/2000));
assert(numel(out.adjacency_tracks)==1);
saved=load(dest);assert(isequal(saved.adjacency_tracks,out.adjacency_tracks));
original=load(input);assert(isequal(original.SR_Localizations,SR_Localizations));
failed=false;try,track_nb_file(input,'PixelSizeUm',[12.32 12.414467],'OutputMAT',dest);catch,failed=true;end
assert(failed,'Existing output must not be overwritten.');
failed=false;try,track_nb_file(input);catch,failed=true;end
assert(failed,'Pixel geometry must be supplied.');
report=struct('passed',true,'checks',{{'bounded frame selection','original timestamps retained',...
    'v7.3 saved output reloads','localizations unchanged','existing output protected','geometry required'}});
fid=fopen(fullfile(root,'workflow_test_results.json'),'w');fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
fprintf('PASS: saved-localization workflow\n');
end
