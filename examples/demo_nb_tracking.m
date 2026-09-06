function output=demo_nb_tracking(inputMAT)
%DEMO_NB_TRACKING Run the tracker on a saved localization MAT file.
% The input contains saved detections; no SVD/localization or IQ load occurs.
if nargin<1
    error('demo_nb_tracking:InputRequired', ...
        'Pass a MAT file containing SR_Localizations as demo_nb_tracking(inputMAT).');
end
code=fileparts(mfilename('fullpath'));src=fullfile(fileparts(code),'src');addpath(src);
output=track_nb_file(inputMAT,'PixelSizeUm',[12.32 12.41446725317693],...
    'FrameRateHz',2000,'MaxFrames',192);
figure('Name','NB tracker: 128 saved L22 frames');ax=axes;hold(ax,'on');
for k=1:numel(output.details.records)
    p=output.details.records(k).positionsPixels;
    % Dots show measured nodes only. Lines connect observations, including gaps.
    plot(ax,p(:,2),p(:,1),'.-','MarkerSize',3,'LineWidth',.5);
end
axis(ax,'equal');set(ax,'YDir','reverse');xlabel(ax,'Lateral pixel');ylabel(ax,'Axial pixel');
title(ax,{'New tracker on a small saved L22 subset',...
    'Retained observed detections; real-data identity accuracy is unknown'});
end
