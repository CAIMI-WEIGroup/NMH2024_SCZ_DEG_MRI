clc, clear, close all

filename = '../results/DESeq2_deg_results_all.csv';
tbl = readtable(filename);

% down-regulated genes
gs_targ{1} = tbl.Row_names((tbl.log2FoldChange < 0) & ...
    (tbl.padj < 0.05));

% up-regulated genes
gs_targ{2} = tbl.Row_names((tbl.log2FoldChange > 0) & ...
    (tbl.padj < 0.05));

% both genes
gs_targ{3} = tbl.Row_names((abs(tbl.log2FoldChange) > 0) & ...
    (tbl.padj < 0.05));

labels = {'down', 'up', 'all'};

% total number
M = size(tbl, 1);

%% load hess and bergon
tt = readtable('../data/DEGs_significant.xlsx');
II = ismember(tt.EnsemblGeneID, gs_targ{3});
tt = tt(II, :);

X1 = nnz(tt.GenesOverlappedWithHessEtAl_)
P = 1 - hygecdf(X1, M, 1807, 2145)

X2 = nnz(tt.GenesOverlappedWithBergonEtAl_)
P = 1 - hygecdf(X2, M, 1807, 247)

%% Load gwas gene sets
% scz gwas
tblgwas = readtable('../data/SCZ_GWAS/genes.txt');
II = ismember(tblgwas.ensg, tbl.Row_names);
tblgwas(~II, :) = [];
gs_ref{1} = tblgwas.ensg; % all significant genes
gs_ref{2} = tblgwas.ensg(tblgwas.posMapSNPs > 0); % positional mapping
gs_ref{3} = tblgwas.ensg(tblgwas.eqtlMapSNPs > 0); % eqtl mapping
gs_ref{4} = tblgwas.ensg(ismember(tblgwas.ciMap, 'Yes')); % ci mapping
gs_label = {'scz_all','scz_pos','scz_eqtl','scz_ci'};

% TWAS scz
twas = readtable('../data/TWAS_scz_gusev.xlsx', 'ReadVariableNames', false);
twas_genes = unique(twas.Var2(ismember(twas.Var11, 1)));
[~, J] = ismember(twas_genes, tmp.Symbol);
J(J==0) = [];
gs_ref{5} = tmp.id(J);
gs_label = [gs_label, 'TWAS'];

% hypergeometric testing
for ii = 1:numel(gs_targ)
    K = numel(gs_targ{ii});
    for jj = 1:numel(gs_ref)
        x(ii,jj) = numel(gs_ref{jj});
        X(ii,jj) = numel(intersect(gs_targ{ii}, gs_ref{jj}));
        N = numel(gs_ref{jj});
        P(ii, jj) = 1 - hygecdf(X(ii,jj), M, K, N);
    end
end

