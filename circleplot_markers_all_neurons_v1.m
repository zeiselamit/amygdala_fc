tic
clear all
close all
addpath('/data/matlab_functions/')

%%
load geneid
load(['amyFC_NONneurons_sc_celltypes_muhammad_25-Nov-2021.mat'])%%%,'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time')
data_nn = data;
cellid_nn = cellid;
c_nn = c;

load(['amyFC_neurons_sc_celltypes_muhammad_25-Nov-2021.mat'])%%%,'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time')
data = [data,data_nn];
cellid = [cellid;cellid_nn];
c = [c;c_nn];


gabahex = loadCellFile_turbo('GABA-hexa.txt',1);
gabahex = gabahex(2:end,:);
glut1hex = loadCellFile_turbo('Glut1-hexa.txt',1);
glut1hex = glut1hex(2:end,:);
glut2hex = loadCellFile_turbo('Glut2-hexa.txt',1);
glut2hex = glut2hex(2:end,:);
nnhex = loadCellFile_turbo('non-neuronal-hexa.txt',1);
nnhex = nnhex(2:end,:);

gabaorder = loadCellFile('ClusterOrder_step3_GABA_FC_04-Nov-2021.txt');
glut1order = loadCellFile('ClusterOrder_step3_Glut1_FC_24-Oct-2021.txt');
glut2order = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');

cuni = [gabaorder(:,2);glut1order(:,2);glut2order(:,2);....
    {'Astro';'Astro_agt';'Astro_SC';'Epend';'OPC_cycling';'OPC';'COP';'OL';'EC';'Peri';'VSM';'VLMC';'microglia';'pvm'}];
r1 = gabaorder(:,2);
r2 = glut1order(:,2);
r3 = glut2order(:,2);
r4 = {'Astro';'Astro_agt';'Astro_SC';'Epend';'OPC_cycling';'OPC';'COP';'OL';'EC';'Peri';'VSM';'VLMC';'microglia';'pvm'};
nn_colors = distinguishable_colors(length(r4));
in = find(sum(data>0,2)>50 );
logflag = 0;
%%
data_med = zeros(length(geneid),length(cuni));
data_mean = zeros(length(geneid),length(cuni));
data_frac = zeros(length(geneid),length(cuni));
data_75 = zeros(length(geneid),length(cuni));
data_90 = zeros(length(geneid),length(cuni));
for i=1:length(cuni)
    i
    gr2 = find(strcmpi(cuni{i},c));
    data_med(:,i) = median(data(:,gr2),2);
    data_mean(:,i) = mean(data(:,gr2),2);
    data_frac(:,i) = mean(data(:,gr2)>0,2);
    data_75(:,i) = prctile(data(:,gr2),75,2);
    data_90(:,i) = prctile(data(:,gr2),90,2);
end
%%
% list = loadCellFile('circleplot-Amy-markers.txt');
% list = list(2:end);
list = loadCellFile('Circleplot-modular.txt');
list = loadCellFile('circleplot-genes_0309.txt');
list = loadCellFile_csv('top_genes_by_cluster_neurons_07-Mar-2022.csv');
list = list(2:end,1);
% list = {'Nos1','Adcyap1','Cck','Grp','Tac1','Tac2','Cartpt','Sst','Pvalb','Vip','Crh','Npy',...
%     'Penk','Pdyn','Nts','Gal','Trh','Avp','Oxt'...
%     ,'Chat','Slc5a7','Slc18a2','Slc18a3','Ddc','Th','Slc6a3','Slc6a9','Slc6a5'...
%     ,'Slc6a1','Pnoc','Gad1','Gad2','Slc32a1','Slc1a1','Slc1a2','Slc1a6','Slc17a7','Slc17a6'};
[~,loc] = ismember(list,geneid);
colorvec = [gabahex(:,3);glut1hex(:,3);glut2hex(:,3)];
colorvec = cellfun(@(x) ['#',x], colorvec,'UniformOutput',0);
colorvec = hex2rgb(colorvec);
hf1 = figure;
set(gcf,'color','w','position',[20,20,1000,900],'Visible','on');
for jjj=1:length(list)
    jjj
    
    g = loc(jjj);%find(strcmpi(gn,geneid));
    gn = list{jjj}
    sz = log2(data_75(g,1:end-14)+1);%.*data_90(g,1:end-14);
    sz = sz/max(sz);
    scatter([1:length(cuni)-14],jjj*ones(1,length(cuni)-14),5*sz+0.001,colorvec,'filled'); hold on;%,'markeredgecolor','k'
%     plot([1,length(cuni)-14],jjj*[1,1],'-k','linewidth',0.2);
%     text(length(cuni)-14,jjj,num2str(round(max(data_75(g,:)))),'HorizontalAlignment','left','fontsize',8)
end
text(-1*ones(length(list),1),[1:length(list)],list,'HorizontalAlignment','right','FontSize',7)
axis off 
set(gca,'ydir','reverse')
% axis equal
eval(['export_fig circleplot_top_genes_by_cluster_',date,'.pdf']);
% eval(['export_fig circleplot_p75_FC_Glut2_',date,'.pdf']);
