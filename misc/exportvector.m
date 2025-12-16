% - EXPORTVECTOR {NAME}
function exportvector(name, fig)
    if (exist('fig', 'var'))
        n = str2double(fig);
        if ~isempty(n)
            figure(n);
        end
    end
    exportgraphics(gcf, ['misc/IGN_graphics/', name, '.pdf'], 'ContentType', 'vector');
end