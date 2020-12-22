function [] = graph_out(F1,graph_name,graph_size)

% Export figure

set(F1,'Units','Inches','position',[0 0 graph_size]);
set(F1,'PaperPositionMode','Auto');
set(F1,'Clipping','off')

print(F1,graph_name,'-dpng');

close(F1);