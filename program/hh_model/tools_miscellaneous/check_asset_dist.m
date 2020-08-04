figure;
subplot(1,2,1)
ksdensity(asset_aux(1,:))
title('unemp')
subplot(1,2,2)
ksdensity(asset_aux(2,:))
title('emp')

figure;
count = 1;
for i = 1:2
    for j = 1:3
        subplot(2,3,count)
        eval(['plot(sim_struct.measureCoefficient_' num2str(i) '_' num2str(j) ')'])
        title(['coef' num2str(i) num2str(j)])
        count = count+1;
    end
end
    