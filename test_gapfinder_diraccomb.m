clear
close all

N = 2000;
m = 5;
height = m^2;

A = generate_dirac_comb(N, m, height);
ev = sort(eig(full(A)));
n = size(A, 1);

% normalize matrix so that spectrum is in [0, 10]:
lmax = max(ev);
lmin = min(ev);
I = speye(n);
A = 10*(A - lmin*I) / (lmax - lmin);

tic;
ee = eig(full(A), 'vector');
t_eig = toc;
ee = sort(ee);
fprintf("Eigenvalue computation time: %.2f\n", t_eig)

mus = linspace(min(ee), max(ee), 1000); 
delta = 1e-2;
its_lanc = 150;
bound_type = "diffsafe";
c = 2;
d = 3;

rng(0);
tic;
[gaps, trest, trest_upper, trest_lower, ~] = gapfinder_main(A, mus, delta, its_lanc, bound_type, d, c);
t_diff = toc;

bound_type = "residue";
rng(0);
tic;
[gaps_res, trest_res, trest_upper_res, trest_lower_res, ~] = gapfinder_main(A, mus, delta, its_lanc, bound_type, d, c);
t_res = toc;

% Find truegaps corresponding to found gaps:
for j = 1:length(gaps)
	mu = 0.5*(gaps{j}(2) + gaps{j}(1));
	idx = find(ee > mu, 1);
	truegaps{j} = [ee(idx-1), ee(idx)];

end

% Plots:
matrixname = "diraccomb";

f1 = figure;
p1 = plot(ee, 'o'); 
p1.MarkerSize = 3;
xlim([-400, length(ee) + 400]);
ylim([min(ee)-.5, max(ee)+.5]);
% yticks([-5, 0, 5, 10, 15]);

ax = gca;
set(ax, 'Units', 'normalized');
set(ax, 'Position', [0.2 0.2 0.7 0.7]);
set(ax, 'LooseInset', [0,0,0,0]);

fname = strcat("plots/",matrixname,"_eig");
print(fname, '-depsc2');

pause(0.1);

f2 = figure;
plot(mus, trest, '-.r'); hold on;
plot(mus, trest_upper, '-.b'); hold on;
plot(mus, trest_lower, '-.b'); hold on;
xlim([min(mus) - .5, max(mus) + .5]);
ylim([min(trest) - 400, max(trest) + 400]);
xlabel("\mu");

ax = gca;
set(ax, 'Units', 'normalized');
set(ax, 'Position', [0.2 0.2 0.7 0.7]);
set(ax, 'LooseInset', [0,0,0,0]);

fname = strcat("plots/gapfinder_",matrixname,"_diff");
print(fname, '-depsc2');

pause(0.1);

f3 = figure;
plot(mus, trest_res, '-.r'); hold on;
plot(mus, trest_upper_res, '-.b'); hold on;
plot(mus, trest_lower_res, '-.b'); hold on;
xlim([min(mus) - .5, max(mus) + .5]);
ylim([min(trest) - 400, max(trest) + 400]);
xlabel("\mu");

ax = gca;
set(ax, 'Units', 'normalized');
set(ax, 'Position', [0.2 0.2 0.7 0.7]);
set(ax, 'LooseInset', [0,0,0,0]);

fname = strcat("plots/gapfinder_",matrixname,"_res");
print(fname, '-depsc2');




% Construct table:
% first gap | second gap | third gap | fourth gap | time
tablerow1 = sprintf("\\texttt{eig} & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", truegaps{[2, 4, 5, 6]}, t_eig);
tablerow2 = sprintf("\\cref{algorithm:eigenvalue-gap-finder--final} w/~consec.~diff. & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", gaps{[2, 4, 5, 6]}, t_diff);
tablerow3 = sprintf("\\cref{algorithm:eigenvalue-gap-finder--final} w/~\\cref{prop:lanczos-a-posteriori-error-bound} & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", gaps_res{[2, 4, 5, 6]}, t_res);
table = sprintf("%s\n%s\n%s", tablerow1, tablerow2, tablerow3);

fname = strcat("tables/table_gapfinder_",matrixname,".txt");
fid = fopen(fname, 'w');
fprintf(fid, "%s\n\n", table);
fclose(fid);
