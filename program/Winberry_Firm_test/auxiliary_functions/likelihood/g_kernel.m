function gg = g_kernel(prod,logk,moment,measureCoefficient)
    
    moment_aux = moment;
    moment_aux(1:2) = 0;
    prod_demean = prod-moment(1);
    logk_demean = logk-moment(2);
    
    nMeasure_all = length(moments);
    aux = nan(1,nMeasure_all);
    counter = 0;
    for i_Measure = 1:nMeasure_all
        aux(counter+(1:(i_Measure+1))) = (prod_demean.^((nMeasure:-1:0)')) ...
            .*(logk_demean.^((0:nMeasure)'));
        counter = counter+i_Measure+1;
    end
    gg = measureCoefficient*(aux-moment_aux');
    
end