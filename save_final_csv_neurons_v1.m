tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 1;
savefig_pdf = 1;


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
    
data = normalize(data,'norm');
data = round(data./repmat(sum(data),length(data(:,1)),1)*20e3);
[median(sum(data(:,batch_flag==1))),median(sum(data(:,batch_flag==2))),median(sum(data(:,batch_flag==3)))]
[median(sum(data(:,batch_flag==1)>0)),median(sum(data(:,batch_flag==2)>0)),median(sum(data(:,batch_flag==3)>0))]

IEG_markers = loadCellFile('IEG_cluster.txt');
not_relevant = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
    ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};%,'Snca','Gm9843','Gm12027','Gm12031'
in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.8 & ~ismember(geneid(:,1),sex_genes)  ....
    & ~ismember(geneid(:,1),[IEG_markers]) & ~ismember(geneid(:,1),[not_relevant]));

corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);
g = unique([geneid(in(corr_filt));geneselection_gaba(:,1);geneselection_glut1(:,1);geneselection_glut2(:,1)]);
[~,corr_filt] = ismember(g,geneid);

% % % % % % % % % % % 
z = linkage(data((corr_filt),:),'ward','correlation');
tic
idx = cluster(z,'maxclust',50);
toc
d = corr_mat(data((corr_filt),:)');
dsum = sum(d>0.2,2);
leaforder = optimalleaforder(z,squareform(1-d,'tovector'));

% [~,xi] = sort(idx);
figure;
set(gcf,'color','w','position',[20,20,900,800])
h1 = axes('position',[0.2,0.1,0.7,0.85]);
imagesc(d(leaforder,leaforder),[0,0.3]);
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',geneid((corr_filt(leaforder))))
freezeColors(h1)
% axis tight
% axis equal
% text([1:length(idx)],[1:length(idx)],num2str(idx(leaforder)))
h2 = axes('position',[0.91,0.1,0.04,0.85]);
imagesc(idx(leaforder));
set(h2,'ydir','reverse','ytick',[]);
colormap('prism');
h3 = axes('position',[0.95,0.1,0.05,0.85]);
plot(sum(d(leaforder,leaforder)>0.2,2),[1:length(leaforder)],'.');
set(h3,'ydir','reverse','ytick',[]);
linkaxes([h1,h2,h3],'y')
% % % % % % % % % % %
% % % % % % % % % % % 

data_orig_all = data;
geneid_all = geneid;
data = data((corr_filt),:);
geneid = geneid((corr_filt));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% 
% % % % % % % % % % % % % % % % % % %
gabaorder = loadCellFile('ClusterOrder_step3_GABA_FC_04-Nov-2021.txt');
glut1order = loadCellFile('ClusterOrder_step3_Glut1_FC_12-Dec-2021.txt');
glut2order = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');
clusteruni = [gabaorder(:,2);glut1order(:,2);glut2order(:,2)];

% clusteruni = unique(cellid_clusters(:,5));
T_cells_tmp = zeros(length(cellid),2);
clustersize = zeros(length(clusteruni),1);
l2cluster = [];
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,5),clusteruni{i}));
    T_cells_tmp(ind,1) = i;
    if ~strcmpi(l2cluster,unique(cellid_clusters(ind,3)))
        l2cluster = [l2cluster;unique(cellid_clusters(ind,3))];
        T_cells_tmp(ind,2) = max(T_cells_tmp(:,2))+1;
    else 
        T_cells_tmp(ind,2) = max(T_cells_tmp(:,2));
    end
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp(:,1));

[~,xi] = sort(T_cells_tmp(:,1));


data_orig_all = data_orig_all(:,xi);
cellid = cellid(xi);
cellid_clusters = cellid_clusters(xi,:);
sample = sample(xi);
fc_time = fc_time(xi);
amy_flag = amy_flag(xi);
batch_flag = batch_flag(xi);
T_cells_tmp = T_cells_tmp(xi);

table1 = [ [cellid';sample';m2c(fc_time)';m2c(batch_flag)';clusteruni(T_cells_tmp)'];m2c(data_orig_all) ];
table1 = [[cell(5,1);geneid_all],table1];

fcell2csv(['Amy_all_final_neurons_UMI_',date,'.csv'],table1);

