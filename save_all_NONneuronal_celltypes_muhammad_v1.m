tic
clear all
close all
addpath('/data/matlab_functions/')

savefig_flag = 0;
savefig_pdf = 0;

load afterloading_noneuronal_QC2000_FC_16-Nov-2021

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);


geneselection_oligos = loadCellFile('gene_correlation_oligos_16-Nov-2021.txt');
geneselection_astro = loadCellFile('gene_correlation_astro_16-Nov-2021.txt');
geneselection_vascular = loadCellFile('gene_correlation_vascular_16-Nov-2021.txt');
geneselection_immune = loadCellFile('gene_correlation_immune_16-Nov-2021.txt');


cellidclsuter_oligos = loadCellFile('cellidCluster_step2_oligos_FC_16-Nov-2021.txt');
oligos_clustername = loadCellFile('clustername_oligos_step2.txt');
cellidclsuter_oligos = [cellidclsuter_oligos,cell(length(cellidclsuter_oligos),1)];
c = cell2mat(cellidclsuter_oligos(:,2));
for i=1:length(oligos_clustername)
    ind = find(c==oligos_clustername{i});
    cellidclsuter_oligos(ind,4) = repmat(oligos_clustername(i,2),length(ind),1);
end
cellidclsuter_astro = loadCellFile('cellidCluster_step2_astro_FC_16-Nov-2021.txt');
astro_clustername = loadCellFile('clustername_astro_step2.txt');
cellidclsuter_astro = [cellidclsuter_astro,cell(length(cellidclsuter_astro),1)];
c = cell2mat(cellidclsuter_astro(:,2));
for i=1:length(astro_clustername)
    ind = find(c==astro_clustername{i});
    cellidclsuter_astro(ind,4) = repmat(astro_clustername(i,2),length(ind),1);
end
cellidclsuter_vascular = loadCellFile('cellidCluster_step2_vascular_FC_16-Nov-2021.txt');
vascular_clustername = loadCellFile('clustername_vascular_step2.txt');
cellidclsuter_vascular = [cellidclsuter_vascular,cell(length(cellidclsuter_vascular),1)];
c = cell2mat(cellidclsuter_vascular(:,2));
for i=1:length(vascular_clustername)
    ind = find(c==vascular_clustername{i});
    cellidclsuter_vascular(ind,4) = repmat(vascular_clustername(i,2),length(ind),1);
end
cellidclsuter_immune = loadCellFile('cellidCluster_step2_immune_FC_16-Nov-2021.txt');
immune_clustername = loadCellFile('clustername_immune_step2.txt');
cellidclsuter_immune = [cellidclsuter_immune,cell(length(cellidclsuter_immune),1)];
c = cell2mat(cellidclsuter_immune(:,2));
for i=1:length(immune_clustername)
    ind = find(c==immune_clustername{i});
    cellidclsuter_immune(ind,4) = repmat(immune_clustername(i,2),length(ind),1);
end


cellid_clusters = cell(length(cellid),4);
[~,loc]= ismember(cellidclsuter_oligos(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_oligos;
[~,loc]= ismember(cellidclsuter_astro(:,1),cellid);
cellid_clusters(loc,:) = cellidclsuter_astro;
[~,loc]= ismember(cellidclsuter_vascular(:,1),cellid);
cellid_clusters(loc,:) = cellidclsuter_vascular;
[~,loc]= ismember(cellidclsuter_immune(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_immune;
rmv = find(strcmpi(cellid_clusters(:,4),'exc'));
cellid_clusters(rmv,:) = cell(length(rmv),4);

validcells = cell2mat(cellfun(@(x) ~isempty(x), cellid_clusters(:,1),'UniformOutput',0));
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
cellid_clusters = cellid_clusters(validcells,:);
sample = sample(validcells);
fc_time = fc_time(validcells);
amy_flag = amy_flag(validcells);
batch_flag = batch_flag(validcells);

[~,xi] = sort(cellid_clusters(:,4));
data = data(:,xi);
cellid = cellid(xi);
cellid_clusters = cellid_clusters(xi,:);
sample = sample(xi);
fc_time = fc_time(xi);
amy_flag = amy_flag(xi);
batch_flag = batch_flag(xi);



sample_uni = {'08-1','23-1', '23-3',  '45-3','46-3', '74-1', '74-2',.....
   '10-1', '18-1', '19-1', '76-1', '76-2',.....
   '42-1', '42-2','46-1', '46-2',.....
   '41-1', '41-2',  '45-1', '45-2',....
   '75-1', '75-2', '75-3'};


c = cellid_clusters(:,4);

save(['amyFC_NONneurons_sc_celltypes_muhammad_',date],'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time','-v7.3')

    