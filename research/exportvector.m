% - EXPORTVECTOR {NAME}
function exportvector(name)
    exportgraphics(gca, ['research/IGN_graphics/', name, '.pdf'], 'ContentType', 'vector');
end