% This script prepare the input file for magma gene property analysis
% i.e., gene-set analysis for continous variables
% the required file is a table with the first col indicating gene id
% and a header of property names (e.g., tissue name / region name)

% --------------------------- settings ---------------------------
% data matrix of the size of #gene by #properties (tissues/regions) 

filename = '../data/DESeq2_deg_results_all.csv';

data = readtable(filename, 'Delimiter', ',');
data.log2FoldChange(isnan(data.padj)) = 0;
datamatrix = [abs(data.log2FoldChange), -log10(data.pvalue)];

tmp = readtable('../data/genes.expression.xls.txt'); % This file is available upon request
[~, J] = ismember(data.Row_names, tmp.id);
gene_symbols = tmp.Symbol(J);

% gene_symbols = data.gene_symbols;
tissue_names = {'log2fc', 'pval'};

% Reference file
reffile = '../data/magma/NCBI37.3.gene.loc';
outputfile = '../data/genecovar.txt'; % This file is available now


% -----------------------------------------------------------------
% Load ref data
tblref = readtable(reffile, 'FileType', 'text');

% based on ensg id
tbl = array2table(datamatrix);
tbl.Properties.RowNames = data.Row_names;
tbl.Properties.VariableNames = tissue_names;

% write table
writetable(tbl, 'tmp.txt', 'WriteRowNames', true, 'Delimiter', '\t');

% replace NaN with NA
fid = fopen('tmp.txt', 'r');
fid1 = fopen(outputfile, 'w');

while ~feof(fid)
   tline = fgetl(fid);
   tline = strrep(tline, 'NaN', 'NA');
   fprintf(fid1, '%s\n', tline);
end
fclose(fid);
fclose(fid1);
delete('tmp.txt');