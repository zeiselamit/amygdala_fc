tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;

sam_cell = {'08_1','08_2','10_1','10_2','18_1','18_2','19_1','23_1','23_2','23_3','41_1','41_2','42_1',...
    '42_2','45_1','45_2','45_3','46_1','46_2','46_3','74_1','74_2','75_1','75_2','75_3','76_1',...
    '76_2'};
sam_folder = {'out10X08_1_220331','out10X08_2_190709','out10X10_1_210219','out10X10_2_210219'....
    ,'out10X18_1_210219','out10X18_2_210219','out10X19_1_210219','out10X23_1_210219','out10X23_2_210219'...
    ,'out10X23_3_210219','out10X41_1_210220','out10X41_2_220330','out10X42_1_220331','out10X42_2_210220',...
    'out10X45_1_220331','out10X45_2_210220','out10X45_3_220331','out10X46_1_210220','out10X46_2_220331',....
    'out10X46_3_210220','out10X74_1_210221','out10X74_2_210221','out10X75_1_210221','out10X75_2_210221',...
    'out10X75_3_210104','out10X76_1_210221','out10X76_2_210221'};
for i=1:length(sam_cell)
    fprintf(['loading sample 10X',sam_cell{i},'\n']);
    mtx_file = ['/data/runs/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/matrix.mtx'];
    bc_file = ['/data/runs/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/barcodes.tsv'];
    gene_file = ['/data/runs/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/features.tsv'];
    eval(['[data_',sam_cell{i},', geneid_',sam_cell{i},', barcodes_',sam_cell{i},'] = load10xMtxFile(mtx_file,bc_file,gene_file,1000,5e4);']);
end
str = [];
for i=1:length(sam_cell)
    str = [str,'repmat({''',sam_cell{i},'''},length(barcodes_',sam_cell{i},'),1);'];
end
eval(['sample = [',str,'];'])
sample = regexprep(sample,'_','-');
sample_uni = unique(sample);
str = [];
for i=1:length(sam_cell)
    str = [str,'data_',sam_cell{i},','];
end
eval(['data = [',str,'];'])
data(data<0) = 0;
str = [];
for i=1:length(sam_cell)
    str = [str,'barcodes_',sam_cell{i},';'];
end
eval(['cellid = [',str,'];'])
geneid = geneid_76_1(:,1);

cellid = cellfun(@(x,y) [x,'_',y], cellid, sample,'UniformOutput',0);

fc_time = zeros(size(cellid));
fc_time(strcmpi(sample,'10-1') | strcmpi(sample,'10-2')....
    | strcmpi(sample,'18-1') | strcmpi(sample,'18-2')  | ....
    strcmpi(sample,'19-1') | strcmpi(sample,'76-1') | strcmpi(sample,'76-2')) = 2;
fc_time(strcmpi(sample,'42-1') | strcmpi(sample,'42-2')....
    | strcmpi(sample,'46-1') | strcmpi(sample,'46-2')) = 8;
fc_time(strcmpi(sample,'41-1') | strcmpi(sample,'41-2')....
    | strcmpi(sample,'45-1') | strcmpi(sample,'45-2')) = 24;
fc_time(strcmpi(sample,'75-1') | strcmpi(sample,'75-2')....
    | strcmpi(sample,'75-3')) = 28;

amy_flag = strcmpi(sample,'08-1') | strcmpi(sample,'10-1') | strcmpi(sample,'23-1') | strcmpi(sample,'23-3') |....
    strcmpi(sample,'18-1') | strcmpi(sample,'19-1') | strcmpi(sample,'42-1') | strcmpi(sample,'42-2')....
    | strcmpi(sample,'46-1') | strcmpi(sample,'46-2') | strcmpi(sample,'46-3') | strcmpi(sample,'41-1') | strcmpi(sample,'41-2')....
    | strcmpi(sample,'45-1') | strcmpi(sample,'45-2') | strcmpi(sample,'45-3')....
    | strcmpi(sample,'74-1') | strcmpi(sample,'74-2') ....
    | strcmpi(sample,'75-1') | strcmpi(sample,'75-2') | strcmpi(sample,'75-3')....
    | strcmpi(sample,'76-1') | strcmpi(sample,'76-2') ;
piri_flag = strcmpi(sample,'08-2') | strcmpi(sample,'10-2') | strcmpi(sample,'23-2') | strcmpi(sample,'18-2');
batch_flag = zeros(size(fc_time));
batch_flag(strcmpi(sample,'08-1') | strcmpi(sample,'08-2') | strcmpi(sample,'10-1')....
    | strcmpi(sample,'10-2') | strcmpi(sample,'18-1') | strcmpi(sample,'18-2') |....
    strcmpi(sample,'19-1') | strcmpi(sample,'23-1') | strcmpi(sample,'23-2') | strcmpi(sample,'23-3')) = 1;
batch_flag(strcmpi(sample,'74-1') | strcmpi(sample,'74-2') | strcmpi(sample,'75-1')....
    | strcmpi(sample,'75-2') | strcmpi(sample,'75-3') | strcmpi(sample,'76-1') | strcmpi(sample,'76-2') |....
    strcmpi(sample,'76-3')) = 3;
batch_flag(batch_flag==0) = 2;

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);

% validcells = (tot_mol>4000 & tot_mol<5e4 & tot_genes>2000 & sum(data(loc,:)>0)<=2 & (stmn2>0 | snap25>0) & amy_flag');
validcells = (tot_mol>3000 & tot_mol<5e4 & tot_genes>2500 & amy_flag'); % all cells
for i=1:length(sample_uni)
   fprintf(['valid cells in ',sample_uni{i},' = ', num2str(sum(validcells(strcmpi(sample,sample_uni{i})))),'\n']);
end
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
sample = sample(validcells);
fc_time = fc_time(validcells);
amy_flag = amy_flag(validcells);
piri_flag = piri_flag(validcells);
batch_flag = batch_flag(validcells);

stmn2 = data(strcmpi(geneid,'Stmn2'),:);
snap25 = data(strcmpi(geneid,'Snap25'),:);
for i=1:length(sample_uni)
   fprintf(['snap25 in before ',sample_uni{i},' = ', num2str(sum(snap25(strcmpi(sample,sample_uni{i}))>0)),'\n']);
end

save(['afterloading_QC3000_FC_',date],'data','geneid','cellid','sample','fc_time','batch_flag','amy_flag','piri_flag','-v7.3');


