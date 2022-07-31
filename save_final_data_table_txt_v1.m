tic
clear all
close all
addpath('/data/matlab_functions/')


load geneid
load(['amyFC_NONneurons_sc_celltypes_muhammad_25-Nov-2021.mat'])%%%,'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time')
data_nn = data;
cellid_nn = cellid;
c_nn = c;
batch_flag_nn = batch_flag;
fc_time_nn = fc_time;
sample_nn = sample;

load(['amyFC_neurons_sc_celltypes_muhammad_25-Nov-2021.mat'])%%%,'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time')
data = [data,data_nn];
data = normalize(data,'norm');
data = round(data./repmat(sum(data),length(data(:,1)),1)*20e3);
cellid = [cellid;cellid_nn];
c = [c;c_nn];
batch_flag = [batch_flag;batch_flag_nn];
fc_time = [fc_time;fc_time_nn];
sample = [sample;sample_nn];


gabahex = loadCellFile_turbo('GABA-hexa.txt',1);
gabahex = gabahex(2:end,:);
glut1hex = loadCellFile_turbo('Glut1-hexa.txt',1);
glut1hex = glut1hex(2:end,:);
glut2hex = loadCellFile_turbo('Glut2-hexa.txt',1);
glut2hex = glut2hex(2:end,:);
nnhex = loadCellFile_turbo('non-neuronal-hexa.txt',1);
nnhex = nnhex(2:end,:);

gabaorder = loadCellFile('ClusterOrder_step3_GABA_FC_04-Nov-2021.txt');
glut1order = loadCellFile('ClusterOrder_step3_Glut1_FC_12-Dec-2021.txt');
glut2order = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');

clusteruni = [gabaorder(:,2);glut1order(:,2);glut2order(:,2)];

cluster_info = loadCellFile('Clusters_info_hannah.txt');
cluster_info(1,:) = [];
clusteruni = [clusteruni,cellfun(@(x,y) [x,'-',y] , cluster_info(:,4),cluster_info(:,5),'UniformOutput',0)];
for i=1:length(clusteruni)
    in = find(strcmpi(clusteruni{i,1},c));
    c(in) = repmat(clusteruni(i,2),length(in),1);
end


table1 = [geneid,m2c(data)];

metadata = [ [{'cellID';'celltype';'sample';'FC time';'batch'}], [cellid';c';sample';m2c([fc_time';batch_flag'])] ];

saveCellFile([metadata;table1],['Amy_FC_allcells_with_metadata_',date,'.txt'])      
