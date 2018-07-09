function F = exp_poly_dist(param)

global nMeasure g2 g3;

mm = param(1:3);
gg = [param(4) g2 g3];
mm_aux = mm;
mm_aux(1) = 0;
fl = @(l) exp(gg*((l-mm(1)).^((1:nMeasure)')-mm_aux'));
fl0 = integral(fl, -Inf, Inf);

F = nan(nMeasure+1,1);
F(1) = integral(@(l) fl(l)*l, -Inf, Inf)/fl0-mm(1);
for i_Measure = 1:nMeasure
    F(i_Measure) = integral(@(l) fl(l)*(l-mm(1))^i_Measure, -Inf, Inf)/fl0-mm(i_Measure);
end
F(end) = integral(@(x) fl(log(x)), 0, Inf)/fl0-1;

