% - EXPORTVECTOR {NAME}
function exportvector(name)
    exportgraphics(gcf, ['misc/IGN_graphics/', name, '.pdf'], 'ContentType', 'vector');
end