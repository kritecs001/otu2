// Run all eight frequency tests using the editable Xcos model.
loadScicos();
loadXcosLibs();

function trace = capture_frequency_trial(trial)
    scicos_simulate(trial,list(),struct(),"nw");
    trace = lab2_freq_trace;
endfunction

folder = get_absolute_file_path("lab2_frequency_sweep.sce");
model = xcosDiagramToScilab(folder+"lab2_frequency.xcos");
[parameters, context_error] = script2var(model.props.context,struct());
if context_error <> 0 then error("Check the Xcos diagram context."); end
period_row = grep(model.props.context,"/^P = /","r");
if size(period_row,"*") <> 1 then error("The diagram context needs one line starting with P =."); end
K = parameters.K;
T = parameters.T;
ds = [parameters.d1,parameters.d2];
if T <= 0 | min(ds) <= 0 then error("Use positive T, d1 and d2 for this frequency experiment."); end
periods = T*[20;10;5;2;1;0.5;0.2;0.1];
frequency = 2*%pi ./ periods;
amplitude = zeros(8,2);
phase = zeros(8,2);
for i = 1:size(periods,1)
    P = periods(i);
    trial = model;
    trial.props.context(period_row) = msprintf("P = %.16g;",P);
    trial.props.tf = (ceil(12*T/(min(ds)*P))+3)*P;
    trace = capture_frequency_trial(trial);
    tail = find(trace.time > trace.time($)-3*P);
    t = trace.time(tail);
    basis = [sin(frequency(i)*t),cos(frequency(i)*t)];
    c = basis \ trace.values(tail,:);
    amplitude(i,:) = sqrt(sum(c(:,2:3).^2,"r")) / norm(c(:,1));
    angle_out = atan(c(2,2:3),c(1,2:3))-atan(c(2,1),c(1,1));
    phase(i,:) = atan(sin(angle_out),cos(angle_out))*180/%pi;
end

window = scf();
window.figure_size = [1000,720];
subplot(2,1,1);
magnitude_db = 20*log10(amplitude);
plot2d("ln",frequency,magnitude_db,style=[-9,-5]);
gca().data_bounds = [min(frequency)/1.4,min(magnitude_db)-4;max(frequency)*1.4,max(magnitude_db)+4];
xtitle("Measured magnitude from Xcos","Angular frequency, rad/s","Magnitude, dB");
legend(msprintf("d1 = %g",ds(1)),msprintf("d2 = %g",ds(2)));
xgrid();
subplot(2,1,2);
plot2d("ln",frequency,phase,style=[-9,-5]);
gca().data_bounds = [min(frequency)/1.4,min(phase)-8;max(frequency)*1.4,max(phase)+8];
xtitle("Measured phase from Xcos","Angular frequency, rad/s","Phase, deg");
legend(msprintf("d1 = %g",ds(1)),msprintf("d2 = %g",ds(2)));
xgrid();

results = [periods,frequency,amplitude(:,1),phase(:,1),amplitude(:,2),phase(:,2)];
mprintf("P (s)       omega (rad/s)    A(d1)       phi(d1)     A(d2)       phi(d2)\n");
for i=1:size(results,1)
    mprintf("%-11.4g %-16.6g %-11.6g %-11.5f %-11.6g %.5f\n",results(i,:));
end
