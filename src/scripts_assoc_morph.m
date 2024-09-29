% Compute associations between gene expression and brain morphology
clc, clear, close
addpath(genpath('./myPLS/'));
addpath(genpath('./gramm/'));
addpath(genpath('./simpleBrainPlot'));

selectIMG = 1; % 1: VOL, 2: SA, 3: CT
isGrouping = 1;
isRunPLS = 0;

% 1. Load data
load('../data/data_all_lausanne120.mat');
dataGE = ge.GE_deg;
dataMRI = squeeze(mri.dataMRI(:, selectIMG, :))';
dataCOV = clin.cov;
group = mri.groupMRI;
subjects = clin.id;

% 2. PLS analysis
if isRunPLS == 0
    load(['../results/results_pls_gene_morph', num2str(selectIMG), '.mat']);

    r = corr(res.X(res.grouping==1, :), res.Y(res.grouping==1, :));
    writematrix(r, '../data/simulation/r_scz.txt', 'Delimiter', ',');
    writematrix(res.X(res.grouping==1, :), '../data/simulation/X_scz.txt', 'Delimiter', ',');
    writematrix(res.Y(res.grouping==1, :), '../data/simulation/Y_scz.txt', 'Delimiter', ',');

    r = corr(res.X(res.grouping==0, :), res.Y(res.grouping==0, :));
    writematrix(r, '../data/simulation/r_hc.txt', 'Delimiter', ',');

else
    X = dataGE;
    Y = dataMRI;
    C = dataCOV;
    G = group;
    
    % exclude genes with more than 3 nans
    II_gene_incl = sum(double(X==0), 1) <= 3;
    X = X(:, II_gene_incl);
    
    % exclude subjects if outliers of expression data
    II_subj_incl = ~isoutlier(sum(X, 2), 'quartiles');
    X = X(II_subj_incl, :);
    Y = Y(II_subj_incl, :);
    C = C(II_subj_incl, :);
    G = G(II_subj_incl, :);
    subjects = subjects(II_subj_incl);
    
    % replace 0 by group median
    for ii = 1:size(X, 2)
        tmp1 = X(G==0, ii);
        tmp1(tmp1==0) = median(nonzeros(tmp1));
        X(G==0, ii) = tmp1;
    
        tmp1 = X(G==1, ii);
        tmp1(tmp1==0) = median(nonzeros(tmp1));
        X(G==1, ii) = tmp1;
    end
        
    % regress out cov from X and Y
    for ii = 1:size(X, 2)
        tmp = regstats(X(:, ii), C, 'linear', {'r', 'tstat'});
        X1(:, ii) = tmp.tstat.beta(1) + tmp.r;
    end
    for ii = 1:size(Y, 2)
        tmp = regstats(Y(:, ii), C, 'linear', {'r', 'tstat'});
        Y1(:, ii) = tmp.tstat.beta(1) + tmp.r;
    end
    
    X1 = X1(randperm(size(X1, 1)), :);

    % PLS correlation analysis
    if isGrouping == 0
        [input, pls_opts, save_opts] = y_pls_input(X1, Y1, G, ...
            '../results/PLS/', 5000, 1000);
    else
        [input, pls_opts, save_opts] = y_pls_input_withGrouping(X1, Y1, G, ...
            '../results/PLS/', 100, 100);
    end
    
    [input, pls_opts, save_opts] = myPLS_initialize(input, pls_opts, save_opts);
    res = myPLS_analysis(input, pls_opts);
    % 
    save(['../results/results_pls_gene_morph', num2str(selectIMG), '.mat'], ...
        'res', 'subjects', 'X', 'Y', 'X1', 'Y1', 'C', 'G', 'II_gene_incl');
end


%% 3. Results: Is there any significant component?
for ii = 1:numel(res.explCovLC)
    sumExpl(ii, 1) = sum(res.explCovLC(1:ii));
end
[~, nComp] = min(pdist2(sumExpl, 0.80)); % 80% explained variance
idx_sigLC = find(mafdr(res.LC_pvals(1:nComp), 'BHFDR', true) < 0.05); % fdr correction

disp('# index of significant LC: ');
disp(idx_sigLC);

disp('How many variance explained by LC1?')
disp(sumExpl(1))

% plot correlation composite scores
clear g;
II = 1; % LC1
grouplabel = cell(size(res.grouping));
grouplabel(res.grouping == 0) = repmat({'HC'}, nnz(res.grouping == 0), 1);
grouplabel(res.grouping == 1) = repmat({'SCZ'}, nnz(res.grouping == 1), 1);
g = gramm('x', res.Lx(:, II), 'y', res.Ly(:, II), 'color', grouplabel);
g.geom_point();
g.stat_glm();
g.set_color_options('map','d3_20');
g.set_point_options('base_size', 3);
g.set_names('x', 'Transcription composite score', ...
    'y', 'Brain composite score', 'color', 'Group');
g.axe_property('FontSize', 8, 'FontName', 'sans-serif');
g.set_text_options('font', 'sans-serif', 'legend_title_scaling', 1, ...
    'legend_scaling', 1);
figure('Unit', 'centimeters', 'Position', [0 0 6.5 5]);
g.draw();
saveas(gcf, ['../figures/PLS_LX_LY_LC',num2str(II),'.svg']);

% HC
disp('within HC: Lx vs Ly');
[r,p] = corr(res.Lx(res.grouping == 0, 1), res.Ly(res.grouping == 0, 1))

% SCZ
disp('within SCZ: Lx vs Ly');
[r,p] = corr(res.Lx(res.grouping == 1, 1), res.Ly(res.grouping == 1, 1))

% plot
y_plot_permut(res.Sp_vect(II, :), res.S(II,II), 'Singluar value', 'Count', ...
    '', 3, 3, ['../figures/PLS_permut_LC',num2str(II),'.svg'], [0, 3000])


%% 3.1 Contributions from X: Z-score 
for ii = 1:numel(idx_sigLC)
    II = idx_sigLC(ii);
    X_salience(:, II) = res.V(:, II);
    X_loadings(:, II) = res.LC_img_loadings(:, II);
    X_zval(:, II) = X_salience(:, II) ./ res.boot_results.Vb_std(:, II);
    X_pval(:, II) = 2 * normcdf(-abs(X_zval(:, II)));
    X_pvaladj(:, II) = mafdr(X_pval(:, II), 'BHFDR', true);
end

% how many significant genes 
nnz(X_pvaladj(:, II) < 0.05) 
nnz(X_pvaladj(:, II) < 0.05) ./ numel(X_pvaladj(:, II) < 0.05) 

nnz((X_pvaladj(:, II) < 0.05) & (X_salience(:, II) > 0))
nnz((X_pvaladj(:, II) < 0.05) & (X_salience(:, II) < 0))


% plot X salience
clear g;    
II = 1;
g = gramm('x', 1:numel(X_salience(:,II)), 'y', X_zval(:, II), ...
    'color', X_pvaladj(:, II) < 0.05);
g.geom_point();
g.set_names('x', 'Genes', 'y', 'z-score', 'color', 'Significance');
g.set_color_options('map',  [0.7,0.7,0.7;    0.1804    0.4824    0.7216],...
            'n_color', 2, 'n_lightness', 1);
g.no_legend();
figure('Unit', 'centimeters', 'Position', [0 0 11 4.5]);
g.draw();
saveas(gcf, ['../figures/PLS_volume_gene_zval_', num2str(II),'.svg']);


%% 3.2 Contributions from Y: Z-score 
for ii = 1:numel(idx_sigLC)
    II = idx_sigLC(ii);

    Y_salience(:, II) = res.U(:, II);
    Y_loadings(:, II) = res.LC_behav_loadings(:, II);    
    Y_zval(:, II) = Y_salience(:, II) ./ res.boot_results.Ub_std(:, II);

    Y_pval(:, II) = 2 * normcdf(-abs(Y_zval(:, II)));
    Y_pvaladj(:, II) = mafdr(Y_pval(:, II), 'BHFDR', true);

    if isGrouping == 0
        clear g;    
        g=gramm('x', [1:numel(Y_zval(:,II))], 'y', Y_zval(:, II), ...
            'color', Y_pvaladj(:, II)<0.05);
        g.geom_point();
        g.set_names('x', 'Regions', 'y', 'z-score', 'color', 'Significance');
        g.set_title('color');
        figure('Position', [100 100 400 300]);
        g.draw();
    else
        clear g;
        tmpgroup = [repmat({'HC'}, size(res.Y, 2), 1); repmat({'SCZ'}, size(res.Y, 2), 1)];
        g=gramm('x', [1:numel(Y_zval(:,II))], 'y', Y_zval(:, II), ...
            'color', tmpgroup);
        g.geom_point();
        g.set_names('x', 'Regions', 'y', 'z-score', 'color', 'Group');
        g.set_title('color');
        figure('Position', [100 100 400 300]);
        g.draw();
    end
end

% plot brain plot
II = 1;
Nregion = size(Y,2);
regionDescriptions = mri.regionDescriptions;
val = Y_zval(1:Nregion, II);
val = val .* (Y_pvaladj(1:Nregion, II) < 0.05);
plotBrain(regionDescriptions, val, flipud(cbrewer('div', 'RdBu', 1000)), ...
    'limits', [-10,10], 'savePath',  ['../figures/PLS_HC_IMG_LC', num2str(II)]);
regionDescriptions(find(val~=0))

val = Y_zval((1+Nregion):end, II);
val = val .* (Y_pvaladj((1+Nregion):end, II) < 0.05);
plotBrain(regionDescriptions, val, flipud(cbrewer('div', 'RdBu', 1000)), ...;
    'limits', [-10,10], 'savePath',  ['../figures/PLS_SCZ_IMG_LC', num2str(II)]);
regionDescriptions(find(val~=0))

val = Y_zval((1+Nregion):end, II);
plotBrain(regionDescriptions, val, flipud(cbrewer('div', 'RdBu', 1000)), ...;
    'limits', [-10,10], 'savePath',  ['../figures/PLS_SCZ_IMG_ALL_LC', num2str(II)]);

save('../results/PLS_deg_morph_LCs.mat', 'Y_salience', 'X_salience', ...
    "Y_loadings", "X_loadings", "Y_zval", "X_zval", "Y_pval", "Y_pvaladj", ...
    "X_pvaladj", "X_pval", "res", "subjects", "II_gene_incl");

table(regionDescriptions, Y_loadings((1+Nregion):end), Y_pvaladj((1+Nregion):end))



%% plot scatter
II = 1;
% region = 'ctx-lh-middletemporal_2';
region = 'ctx-rh-rostralmiddlefrontal_2';
% region = 'ctx-lh-cuneus_1';

[~, III] = ismember(region, regionDescriptions);
grouplabel = cell(size(res.grouping));
grouplabel(res.grouping == 0) = repmat({'HC'}, nnz(res.grouping == 0), 1);
grouplabel(res.grouping == 1) = repmat({'SCZ'}, nnz(res.grouping == 1), 1);
g = gramm('x', res.Lx(:, II), 'y', res.Y(:, III), 'color', grouplabel);
g.geom_point();
g.stat_glm();
g.set_color_options('map','d3_20');
g.set_point_options('base_size', 3);
g.set_names('x', 'Transcription composite score', ...
    'y', 'Adjusted brain volume', 'color', 'Group');
g.axe_property('FontSize', 9, 'Font', 'sans-serif');
figure('Unit', 'centimeters', 'Position', [0 0 6.5 4.5]);
g.draw();
saveas(gcf, ['../figures/PLS_morph_Y_LY_LC',num2str(II),'_',region,'.svg']);

X = res.Lx(:, II);
Y = res.Y(:, III);
[r,p] = corr(X(res.grouping == 1),Y(res.grouping == 1))
return




%% NOT REPORTED
% nsdata = load('~/projects/project_gamba/data/neurosynth/results/result_make_neurosynth_QC_05092019.mat');
nsdata = load('~/projects/project_gamba/data/neurosynth/results/result_brainmapfMRI_QC_05092019.mat');

XX = Y_loadings(115:end, 1);

for ii = 1:numel(nsdata.features)
    [rr(ii,1), pp(ii,1)] = corr(squeeze(nsdata.regionFeatures(:,2,ii)), XX,...
        'rows', 'complete');
end
[rrs, rid] = sort(rr,'descend');
tt = table(rrs, pp(rid), mafdr(pp(rid), 'BHFDR', true), nsdata.features(rid));
disp(head(tt));

[~, J] = ismember('hallucination', tt.Var3)


%% permut analysis
isRunPermt = true;
if isRunPermt
    load(['../results/PLS/results_pls_morphology_null_', num2str(selectIMG), '.mat']);
    LC_img_null_HC = [];    
    LC_img_null_SZ = [];

    for i = 1:numel(resnull)
        tmpnull = resnull{i};
        LC_img_null_HC(:, i) = tmpnull.LC_behav_loadings(1:114, 1);
        LC_img_null_SZ(:, i) = tmpnull.LC_behav_loadings(115:end, 1);
    end
    
    for i = 1:size(LC_img_null_HC, 1)
        PP_HC(i, 1) = nnz(LC_img_null_HC(i, :) > res.LC_behav_loadings(i, 1)) ./ size(LC_img_null_HC, 2);
        PP_SZ(i, 1) = nnz(LC_img_null_SZ(i, :) > res.LC_behav_loadings(i+114, 1)) ./ size(LC_img_null_HC, 2);
    end
    
    plotBrain(regionDescriptions, PP_HC, cbrewer('seq', 'Reds', 1000));
    plotBrain(regionDescriptions, PP_SZ, cbrewer('seq', 'Reds', 1000));
end    


%% PLS null-random gene
isRunPermt = false;
if isRunPermt
    for i = 1:1000
        IIgene_rand = IIgene(randperm(numel(IIgene)));
        GE_rand = GE_all(:, IIgene_rand);
        dataGErand = GE_rand(ia, :);

        X = dataGErand;
        Y = dataMRI;
        C = cov;
        G = group;

        II_gene_incl = sum(double(X==0), 1) <= 3;
        X = X(:, II_gene_incl);

        II_subj_incl = ~isoutlier(sum(X, 2), 'quartiles');
        X = X(II_subj_incl, :);
        Y = Y(II_subj_incl, :);
        C = C(II_subj_incl, :);
        G = G(II_subj_incl, :);

        % replace 0 by group median
        for ii = 1:size(X, 2)
            tmp1 = X(G==0, ii);
            tmp1(tmp1==0) = median(nonzeros(tmp1));
            X(G==0, ii) = tmp1;
            tmp1 = X(G==1, ii);
            tmp1(tmp1==0) = median(nonzeros(tmp1));
            X(G==1, ii) = tmp1;
        end

        % Regress out cov from X and Y
        X1 = [], Y1 = [];
        for ii = 1:size(X, 2)
            tmp = regstats(X(:, ii), C, 'linear', {'r', 'tstat'});
            X1(:, ii) = tmp.tstat.beta(1) + tmp.r;
        end
        for ii = 1:size(Y, 2)
            tmp = regstats(Y(:, ii), C, 'linear', {'r', 'tstat'});
            Y1(:, ii) = tmp.tstat.beta(1) + tmp.r;
        end

        % PLS correlation analysis
        if isGrouping == 0
            [input, pls_opts, save_opts] = y_pls_input(X1, Y1, G, '../results/PLS/', 2, 2);
        else
            [input, pls_opts, save_opts] = y_pls_input_withGrouping(X1, Y1, G, '../results/PLS/', 2, 2);
        end

        [input, pls_opts, save_opts] = myPLS_initialize(input, pls_opts, save_opts);
        resnull{i, 1} = myPLS_analysis(input, pls_opts);
        resnull{i, 1}.boot_results = [];        
    end
    
    save(['../results/PLS/results_pls_morphology_null_', num2str(selectIMG), '.mat'], 'resnull', '-v7.3');
end

return
