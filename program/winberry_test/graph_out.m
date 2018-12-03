function [] = graph_out(F1,graph_tag,label,SFont,graph_size,lgd)
figure(F1)

set(gca,'FontSize',SFont,'FontWeight','bold')

if isempty(lgd) == 0
    h1 = legend(lgd,...
        'orientation','horizontal','color','none','FontSize',SFont);
    legend('boxoff')
    
    p = get(h1,'Position');
    p(1) = 0.2;
    p(2) = 0;
    set(h1,'Position',p)
end

set(F1,'Units','Inches','position',[0 0 graph_size]);
set(F1,'PaperPositionMode','Auto');
set(F1,'Clipping','off')

print(F1,[graph_tag label],'-dpng');