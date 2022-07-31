tic
clear all
close all
addpath('/data/matlab_functions/')

savefig_flag = 0;
savefig_pdf = 0;

sam_cell = {'08_1','08_2','10_1','10_2','18_1','18_2','19_1','23_1','23_2','23_3','41_1','41_2','42_1',...
    '42_2','45_1','45_2','45_3','46_1','46_2','46_3','74_1','74_2','75_1','75_2','75_3','76_1',...
    '76_2'};
sam_folder = {'out10X08_1_210219','out10X08_2_190709','out10X10_1_210219','out10X10_2_210219'....
    ,'out10X18_1_210219','out10X18_2_210219','out10X19_1_210219','out10X23_1_210219','out10X23_2_210219'...
    ,'out10X23_3_210219','out10X41_1_210220','out10X41_2_200318','out10X42_1_200319','out10X42_2_210220',...
    'out10X45_1_200319','out10X45_2_210220','out10X45_3_200319','out10X46_1_210220','out10X46_2_200319',....
    'out10X46_3_210220','out10X74_1_210221','out10X74_2_210221','out10X75_1_210221','out10X75_2_210221',...
    'out10X75_3_210104','out10X76_1_210221','out10X76_2_210221'};
for i=1:length(sam_cell)
    fprintf(['loading sample 10X',sam_cell{i},'\n']);
    mtx_file = ['/data/runs/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/matrix.mtx'];
    bc_file = ['/data/runs/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/barcodes.tsv'];
    gene_file = ['/data/runs/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/features.tsv'];
    command1 = ['cp ',mtx_file,' /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/'];
    command2 = ['mv /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/matrix.mtx /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/',sam_cell{i},'_matrix.mtx'];
    eval(['system(''',command1,''')']);
    eval(['system(''',command2,''')']);
    command1 = ['cp ',bc_file,' /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/'];
    command2 = ['mv /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/barcodes.tsv /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/',sam_cell{i},'_barcodes.tsv'];
    eval(['system(''',command1,''')']);
    eval(['system(''',command2,''')']);
    command1 = ['cp ',gene_file,' /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/'];
    command2 = ['mv /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/features.tsv /data/Technion_analysis/Amygdala/raw_cellranger_files_Amy_220306/',sam_cell{i},'_features.tsv'];
    eval(['system(''',command1,''')']);
    eval(['system(''',command2,''')']);
end













