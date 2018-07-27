function gg = g_kernel(prod,logk,moment,measureCoefficient,nMeasure)
    
    moment_aux = moment;
    moment_aux(1:2) = 0;
    prod_demean = prod-moment(1);
    logk_demean = logk-moment(2);
    
    gg = zeros(size(prod));
    counter = 0;
    for i_Measure = 1:nMeasure
        for j_Measure = 1:(i_Measure+1)
            gg = gg+measureCoefficient(counter+j_Measure)...
                *((prod_demean.^(i_Measure-j_Measure+1)) ...
                .*(logk_demean.^(j_Measure-1))-moment_aux(counter+j_Measure));
        end
        counter = counter+i_Measure+1;
    end
    
end