clear
close all

load DFT-examples/h2_chain_1250mol_6-31G_lda_matrices.mat
Z=chol(S);
% 

% 
A = (Z')\H_DFT/Z;  % Hamiltonian w.r.t. an orthogonal basis
% 
Z = sparse(Z);
H_DFT = sparse(H_DFT);
% 
tic;
ee = eig(full(A), 'vector');
ee = sort(ee);		% should coincide with eigs_Ha
t_eig = toc;
fprintf("Eigenvalue computation time: %.2f\n", t_eig)

% Define function handle:
Afun = @(x) Z' \ (H_DFT * (Z \ x));
A = matrix_handle(Afun, size(Z,1));
n = size(A, 1);

matrixname = strcat("H2chain",int2str(n));

rng(0);
fprintf("n = %d\n", n);

mus = linspace(min(ee), max(ee), 1000); 
delta = 1e-2;
its_lanc = 100;
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
close all
figure;
p1 = plot(ee, 'o'); 
p1.MarkerSize = 3;
xlim([-100, length(ee) + 100]);
ylim([min(ee)-.1, max(ee)+.1]);
% yticks([-5, 0, 5, 10, 15]);

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = strcat("plots/",matrixname,"_eig");
print(fname, '-depsc2');

pause(0.2)

figure;
plot(mus, trest, '-.r'); hold on;
plot(mus, trest_upper, '-.b'); hold on;
plot(mus, trest_lower, '-.b'); hold on;
xlim([min(mus) - .1, max(mus) + .1]);
ylim([min(trest) - 100, max(trest) + 100]);
xlabel("\mu");

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = strcat("plots/gapfinder_",matrixname,"_diff");
print(fname, '-depsc2');

pause(0.2)

figure;
plot(mus, trest_res, '-.r'); hold on;
plot(mus, trest_upper_res, '-.b'); hold on;
plot(mus, trest_lower_res, '-.b'); hold on;
xlim([min(mus) - .1, max(mus) + .1]);
ylim([min(trest) - 100, max(trest) + 100]);
xlabel("\mu");

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize', [6 4]);
set(gcf, "PaperPosition", [0 0 6 4]);

fname = strcat("plots/gapfinder_",matrixname,"_res");
print(fname, '-depsc2');

% Construct table:
% first gap | second gap | third gap | time
tablerow1 = sprintf("\\texttt{eig} & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", truegaps{1:3}, t_eig);
tablerow2 = sprintf("\\cref{algorithm:eigenvalue-gap-finder--final} w/~consec.~diff. & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", gaps{1:3}, t_diff);
tablerow3 = sprintf("\\cref{algorithm:eigenvalue-gap-finder--final} w/~\\cref{prop:lanczos-a-posteriori-error-bound} & [%.3f, %.3f] & [%.3f, %.3f] & [%.3f, %.3f] & %.3f \\\\", gaps_res{1:3}, t_res);
table = sprintf("%s\n%s\n%s", tablerow1, tablerow2, tablerow3);

fname = strcat("tables/table_gapfinder_",matrixname,".txt");
fid = fopen(fname, 'w');
fprintf(fid, "%s\n\n", table);
fclose(fid);
