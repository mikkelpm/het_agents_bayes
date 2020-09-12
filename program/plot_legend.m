% plot legend for layout 2, no need to be in final package

co = [zeros(1,3); get(0, 'DefaultAxesColorOrder')];
SFont = 12;
graph_size = [6 .5];

F1 = figure;
hold on
for i_model = 1:5
    plot(0,0,'LineWidth',2,'Color',co(i_model,:))
end
hold off

axis([1 2 1 2]) 
axis off 

h1 = legend({'FI','macro','3 mom','2 mom','1 mom'},... % change if needed
            'orientation','horizontal','color','none','FontSize',SFont);
legend('boxoff')
h1.Position(1) = 0.5 - h1.Position(3)/2; 
h1.Position(2) = 0.5 - h1.Position(4)/2;

graph_out(F1,'results/hh_lik_legend_layout2',graph_size)