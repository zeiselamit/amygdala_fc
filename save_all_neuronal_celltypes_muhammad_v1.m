tic
clear all
close all
addpath('/data/matlab_functions/')

load afterloading_QC3000_FC_06-Oct-2021

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);

cellidclsuter_gaba = loadCellFile('cellidCluster_step2_GABA_FC_07-Oct-2021.txt');
cellidclsuter_gaba = [cellidclsuter_gaba, cellfun(@(x,y) ['gaba_',num2str(x),'_',num2str(y)],cellidclsuter_gaba(:,2),cellidclsuter_gaba(:,3),'UniformOutput',0)]; 
cellidclsuter_glut1 = loadCellFile('cellidCluster_step2_glut1_FC_07-Oct-2021.txt');
cellidclsuter_glut1 = [cellidclsuter_glut1, cellfun(@(x,y) ['glut1_',num2str(x),'_',num2str(y)],cellidclsuter_glut1(:,2),cellidclsuter_glut1(:,3),'UniformOutput',0)]; 
cellidclsuter_glut2 = loadCellFile('cellidCluster_step2_Glut2_FC_07-Oct-2021.txt');
cellidclsuter_glut2 = [cellidclsuter_glut2, cellfun(@(x,y) ['glut2_',num2str(x),'_',num2str(y)],cellidclsuter_glut2(:,2),cellidclsuter_glut2(:,3),'UniformOutput',0)]; 
geneselection_gaba = loadCellFile('gene_correlation_GABA_07-Oct-2021.txt');
geneselection_glut1 = loadCellFile('gene_correlation_Glut1_07-Oct-2021.txt');
geneselection_glut2 = loadCellFile('gene_correlation_Glut2_07-Oct-2021.txt');

gaba_clustername = loadCellFile('FC_GABA_cluster_list_211007.txt');
gaba_clustername(:,2) = regexprep(gaba_clustername(:,2),'_gaba','');
gaba_clustername(:,2) = cellfun(@(x) ['GABA_',x], gaba_clustername(:,2),'UniformOutput',0);
gaba_clustername(:,2) = regexprep(gaba_clustername(:,2),'GABA_exc','exc');
glut1_clustername = loadCellFile('FC_Glut1_cluster_list_211007');
glut1_clustername(:,2) = regexprep(glut1_clustername(:,2),'_glut1','');
glut1_clustername(:,2) = cellfun(@(x) ['Glut1_',x], glut1_clustername(:,2),'UniformOutput',0);
glut1_clustername(:,2) = regexprep(glut1_clustername(:,2),'Glut1_exc','exc');
glut2_clustername = loadCellFile('FC_Glut2_cluster_list_211007');
glut2_clustername(:,2) = regexprep(glut2_clustername(:,2),'_glut2','');
glut2_clustername(:,2) = cellfun(@(x) ['Glut2_',x], glut2_clustername(:,2),'UniformOutput',0);
glut2_clustername(:,2) = regexprep(glut2_clustername(:,2),'Glut2_exc','exc');

cellidclsuter_gaba = [cellidclsuter_gaba,cell(length(cellidclsuter_gaba),1)];
c = cell2mat(cellidclsuter_gaba(:,2));
for i=1:length(gaba_clustername)
    ind = find(c==gaba_clustername{i});
    cellidclsuter_gaba(ind,5) = repmat(gaba_clustername(i,2),length(ind),1);
end
cellidclsuter_glut1 = [cellidclsuter_glut1,cell(length(cellidclsuter_glut1),1)];
c = cell2mat(cellidclsuter_glut1(:,2));
for i=1:length(glut1_clustername)
    ind = find(c==glut1_clustername{i});
    cellidclsuter_glut1(ind,5) = repmat(glut1_clustername(i,2),length(ind),1);
end
cellidclsuter_glut2 = [cellidclsuter_glut2,cell(length(cellidclsuter_glut2),1)];
c = cell2mat(cellidclsuter_glut2(:,2));
for i=1:length(glut2_clustername)
    ind = find(c==glut2_clustername{i});
    cellidclsuter_glut2(ind,5) = repmat(glut2_clustername(i,2),length(ind),1);
end

cellid_clusters = cell(length(cellid),5);
[~,loc]= ismember(cellidclsuter_gaba(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_gaba;
[~,loc]= ismember(cellidclsuter_glut1(:,1),cellid);
cellid_clusters(loc,:) = cellidclsuter_glut1;
[~,loc]= ismember(cellidclsuter_glut2(:,1),cellid);
cellid_clusters(loc,:) = cellidclsuter_glut2;
rmv = find(strcmpi(cellid_clusters(:,5),'exc'));
cellid_clusters(rmv,:) = cell(length(rmv),5);

validcells = cell2mat(cellfun(@(x) ~isempty(x), cellid_clusters(:,1),'UniformOutput',0));
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
cellid_clusters = cellid_clusters(validcells,:);
sample = sample(validcells);
fc_time = fc_time(validcells);
amy_flag = amy_flag(validcells);
batch_flag = batch_flag(validcells);

[~,xi] = sort(cellid_clusters(:,5));
data = data(:,xi);
cellid = cellid(xi);
cellid_clusters = cellid_clusters(xi,:);
sample = sample(xi);
fc_time = fc_time(xi);
amy_flag = amy_flag(xi);
batch_flag = batch_flag(xi);
gaba_glut_flag = zeros(size(cellid));
gaba_glut_flag(cell2mat(cellfun(@(x) contains(x,'GABA_'),cellid_clusters(:,5),'UniformOutput',0))) = 1;
gaba_glut_flag(cell2mat(cellfun(@(x) contains(x,'Glut1_'),cellid_clusters(:,5),'UniformOutput',0))) = 2;
gaba_glut_flag(cell2mat(cellfun(@(x) contains(x,'Glut2_'),cellid_clusters(:,5),'UniformOutput',0))) = 3;

cellid_clusters(gaba_glut_flag==1,3) = cellfun(@(x) ['GABA_',num2str(x)], cellid_clusters(gaba_glut_flag==1,3),'UniformOutput',0);
cellid_clusters(gaba_glut_flag==2,3) = cellfun(@(x) ['Glut1_',num2str(x)], cellid_clusters(gaba_glut_flag==2,3),'UniformOutput',0);
cellid_clusters(gaba_glut_flag==3,3) = cellfun(@(x) ['Glut2_',num2str(x)], cellid_clusters(gaba_glut_flag==3,3),'UniformOutput',0);

sample_uni = {'08-1','23-1', '23-3',  '45-3','46-3', '74-1', '74-2',.....
   '10-1', '18-1', '19-1', '76-1', '76-2',.....
   '42-1', '42-2','46-1', '46-2',.....
   '41-1', '41-2',  '45-1', '45-2',....
   '75-1', '75-2', '75-3'};
c = cellid_clusters(:,5);

save(['amyFC_neurons_sc_celltypes_muhammad_',date],'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time','-v7.3')


