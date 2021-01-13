function subroutine_fixReference(unfixed_reference_filename)
% subroutine_fixReference This function rotates any image and translates it
% so that it is centered. The user has to manually draw a line following
% the medial the line and click it twice to signal it finished. After, the
% user has to encircle the whole window. 
%
% Developed by SAM on 11/15/2020. Based on code developed by LFM.
% Last updated 12/27/2020.

unfixed_reference = imread(unfixed_reference_filename);
mouse = strsplit(unfixed_reference_filename, {'.','_'});
mouse = mouse{1};

%% Figure options 
fig_opt = {'Color', [0 0 0], 'Position', [200 200 725 725]};

ax_opt = {'xticklabel',[],'yticklabel',[], 'XColor',[0 0 0], 'YColor',[0 0 0],...
    'Colormap',colormap('gray'), 'DataAspectRatio', [1 1 1]};
close(gcf)
t_opt = {'FontName', 'CMU Bright', ...
    'Color', [1 1 1], 'FontSize', 30};

%%
isSuccess = false;
while ~isSuccess
    
    figure(fig_opt{:}, 'Position', [200 200 1200 600]),
    subplot(1, 2, 1)
    imagesc(unfixed_reference),
    set(gca, ax_opt{:}),
    title('Original Reference', t_opt{:});
    
    try
        h = drawline('Position',h.Position,'Color','w');
    catch
        h = drawline('Color','w');
    end
    
    customWait(h);
    
    % Rotation
    line_param = polyfit(h.Position(:,1), h.Position(:,2), 1);
    angle = 90 - atand(line_param(1));
    fixed_reference = imrotate(unfixed_reference, -angle + 180, 'crop');
     
    % Centering
    try
        c = drawcircle('Center',c.Center,'Radius',c.Radius, ...
            'Color', 'w', 'FaceAlpha', 0.4);
    catch
        c = drawcircle('Color', 'w', 'FaceAlpha', 0.4);
    end
    
    customWait(c);
    
    translation = round(size(fixed_reference) / 2) - c.Center;
    fixed_reference = imtranslate(fixed_reference, translation);
    
    subplot(1, 2, 2)
    imagesc(fixed_reference),
    set(gca, ax_opt{:}),
    title('Processed', t_opt{:});
     
    pause(1.5)
    
    repetition_flag = questdlg('Is everything as expected?',...
        'Checkpoint','Yes', 'No', 'Yes');
    
    if strcmp(repetition_flag,'Yes')
        isSuccess = true;
    end
    
end
close(gcf)
imwrite(fixed_reference, strcat(mouse,'_ref_fixed.tif'),'tif');

    function pos = customWait(hROI)
        
        % Listen for mouse clicks on the ROI
        l = addlistener(hROI,'ROIClicked',@clickCallback);
        
        % Block program execution
        uiwait;
        
        % Remove listener
        delete(l);
        
        % Return the current position
        pos = hROI.Position;
        
    end

    function clickCallback(~,evt)
        
        if strcmp(evt.SelectionType,'double')
            uiresume;
        end
        
    end

end