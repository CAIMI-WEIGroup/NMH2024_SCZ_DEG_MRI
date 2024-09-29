function y_plot_permut(Xnull, Xraw, xname, yname, fig_title, xsize, ysize, outname, xlim)

if nargin <=8
    isXlim = 0;
else
    isXlim = 1;
end

clear g
g = gramm('x', Xnull);
g.stat_bin('geom', 'stairs', 'fill', 'edge'); % 'edges', (min(Xnull)-0.02):0.01:(max(Xnull)+0.02)Default fill is edge
g.set_color_options('chroma', 0, 'lightness', 40); 
g.axe_property('FontSize', 6); % , 'XGrid', 'on', 'YGrid', 'on'
g.set_names('x', xname, 'y', yname);   
g.geom_vline('xintercept', Xraw, 'style', 'r--');
g.set_text_options('font', 'sans-serif');

if isXlim == 1
    g.axe_property('XLim', xlim);
end

if ~isempty(fig_title)
    g.set_title(fig_title, "FontSize", 9, "FontName", 'Arial');
end

figure('Unit', 'centimeters', 'Position', [0 0 xsize ysize]);
g.draw();
    
saveas(gcf, outname);
end