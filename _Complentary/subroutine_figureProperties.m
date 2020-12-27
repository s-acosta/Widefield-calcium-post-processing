function [fig_opt, ax_opt, t_opt] = subroutine_figureProperties()

fig_opt = {'Color', [0 0 0], 'Position', [200 200 725 725]};

ax_opt = {'xticklabel',[],'yticklabel',[], 'XColor',[0 0 0], 'YColor',[0 0 0],...
    'Colormap',colormap('gray'), 'DataAspectRatio', [1 1 1]};
close(gcf)
t_opt = {'FontName', 'CMU Bright', ...
    'Color', [1 1 1], 'FontSize', 30};

end