tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;


% load afterloading_QC3000_FC_06-Oct-2021
load afterloading_QC3000_FC_31-Mar-2022

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);
% %%
% cellidclsuter_gaba = loadCellFile('cellidCluster_step3_GABA_FC_24-Oct-2021.txt');
% gaba_clustername = loadCellFile('ClusterOrder_step3_GABA_FC_24-Oct-2021.txt');
% cellidclsuter_gaba = [cellidclsuter_gaba,cell(length(cellidclsuter_gaba),1)];
% c = cell2mat(cellidclsuter_gaba(:,2));
% for i=1:length(gaba_clustername)
%     ind = find(c==gaba_clustername{i});
%     cellidclsuter_gaba(ind,end) = repmat(gaba_clustername(i,2),length(ind),1);
% end
% 

%%
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
moldata = data;

datalog_tmp = cent_norm([(log2(moldata+1))]);
data_tsne = datalog_tmp;
% data_tsne = cent_norm(log2(datamarkers+1));
initial_dims = length(corr_filt);
[prj,m,D,V,Q] = pca_wis(data_tsne',initial_dims);
D = diag(D);
initial_dims = findknee(D);
figure;
subplot(1,2,1)
plot(cumsum(D)); hold on;
plot(initial_dims,sum(D(1:initial_dims)),'sk');
subplot(1,2,2)
plot(D); hold on;
plot(initial_dims,D(initial_dims),'sk');
title(['opt PC = ',num2str(initial_dims)]);
% initial_dims = 30;
prj = prj(:,1:initial_dims);
init = prj(:,1:2)/std(prj(:,1))*1e-4;
init = init-repmat(mean(init),length(init),1);
% % % % % % % % % % % % % % % % % % % % % test harmony py
s = prj;
batchid = batch_flag;
usepylib = 1;
[sout]=harmonypy(s,batchid,usepylib);
prj = sout;
% % % % % % % % % % % % % % % % % % % 
% D = squareform(pdist(prj,'correlation'),'tomatrix');
% [Dsort,XI] = sort(D,'ascend');
% per_range = 500;
% x = 1:per_range;%length(D);
% optk = zeros(length(D),1);
% for i=1:length(D)
%     y = Dsort(1:per_range,i);
%     x = x(:);
%     y = y(:);
%     a = atan((y(end)-y(1))/(x(end)-x(1)));
%     xn = x*cos(a) + y*sin(a);
%     yn = -x*sin(a) + y*cos(a);
%     [~,imax] = max(yn);
%     optk(i) = round((x(imax)));
% end
% 
% perplexity = median(optk);

% options = statset('MaxIter',1000);
% mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
%     'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
% toc
% 
% % % % % % % % % % % % % % % % % % %
gabaorder = loadCellFile('ClusterOrder_step3_GABA_FC_04-Nov-2021.txt');
glut1order = loadCellFile('ClusterOrder_step3_Glut1_FC_12-Dec-2021.txt');
glut2order = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');
clusteruni = [gabaorder(:,2);glut1order(:,2);glut2order(:,2)];

cluster_info = loadCellFile('Clusters_info_hannah.txt');
cluster_info(1,:) = [];
clusteruni = [clusteruni,cellfun(@(x,y) [x,'-',y] , cluster_info(:,4),cluster_info(:,5),'UniformOutput',0)];

tmp = loadCellFile_turbo('GABA-hexa.txt',1);
tmp = tmp(2:end,2:3);
tmp1 = loadCellFile_turbo('Glut1-hexa.txt',1);
tmp1 = tmp1(2:end,2:3);
clustername_color_hannah = loadCellFile_turbo('Glut2-hexa.txt',1);
clustername_color_hannah = clustername_color_hannah(2:end,2:3);
clustername_color_hannah = [tmp;tmp1;clustername_color_hannah];
[~,loc] = ismember(clusteruni(:,1),clustername_color_hannah(:,1));
clusteruni = [clusteruni,clustername_color_hannah(loc,2)];

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
T_cells_tmp_uni_l2 = unique(T_cells_tmp(:,2));
saveCellFile([clusteruni,m2c(clustersize)],['clustername_all_neurons_',date,'.txt'])
% % % % % % % % % % % % % % % % % % % 
%%
fc_time_uni = unique(fc_time);
sum_tp_cluster = zeros(length(clusteruni),length(fc_time_uni));
sum_sample_cluster = zeros(length(clusteruni),length(sample_uni));
for i=1:length(T_cells_tmp_uni)
    ind = T_cells_tmp(:,1)==T_cells_tmp_uni(i);
    for jj=1:length(fc_time_uni)
        sum_tp_cluster(i,jj) = sum(fc_time(ind)==fc_time_uni(jj));
    end
    for jj=1:length(sample_uni)
        sum_sample_cluster(i,jj) = sum(strcmpi(sample(ind),sample_uni{jj}));
    end
end
sum_tp_cluster_frac = sum_tp_cluster./repmat(sum(sum_tp_cluster,2),1,length(fc_time_uni));
figure;
set(gcf,'color','w','position',[20,20,900,800])
subplot(2,1,1);
bar(sum_tp_cluster_frac,'stacked');
legend({'0','2','8','24','recall'})

cm = [[linspace(0,1,length(sample_uni))',linspace(0,1,length(sample_uni))',linspace(1,1,length(sample_uni))']];
cm = colormap('parula');
cm = cm(round(linspace(1,length(cm),length(sample_uni))),:);
sum_sample_cluster_frac = sum_sample_cluster./repmat(sum(sum_sample_cluster,2),1,length(sample_uni));
subplot(2,1,2)
b = bar(sum_sample_cluster_frac,'stacked','FaceColor','flat');
for i=1:length(sample_uni)
    b(i).CData = cm(i,:);
end
legend(sample_uni)
if savefig_pdf==1
    eval(['export_fig bar_frac_sample_tp_per_cluster_',date,'.pdf']);
end
%%
% colors = distinguishable_colors(length(unique(T_cells_tmp))+1);
% figure;
% set(gcf,'color','w','position',[20,20,900,800])
% for i=unique(T_cells_tmp)'
%     if i==0
%         ii=find(T_cells_tmp==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'o','color',colors(i+1,:),'markersize',3); hold on;
%     elseif i>0
%         ii=find(T_cells_tmp==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i+1,:),'markersize',4); hold on;
%     end
% end
% axis tight;
% xl = get(gca,'xlim');
% yl = get(gca,'ylim');
% for i=T_cells_tmp'
%     if i>=0
%         in = T_cells_tmp==i;
%         ht = text(median(mapped_xy(in,1))-0.02*diff(xl),median(mapped_xy(in,2))-0.02*diff(yl),regexprep(clusteruni{i},'_','-'));
%         set(ht,'fontsize',8,'FontName','arial');%,'BackgroundColor',0.8*[1,1,1]
%     end
% end
% axis tight;
% axis equal
% axis off
% 
% title(['#C=',num2str(max(T_cells_tmp))],'fontsize',8);
% % 

% % % % % % % % % % % % %
% list = {'Slc17a6','Slc17a7','Slc32a1','Gad2','Sim1','Trhr','Trh','Penk','Tshz2','Ptk2b','Oxtr','Cxcl14','Sst','C1ql1','Rorb','Gpr88','Neurod6'};
% figure;
% set(gcf,'color','w','position',[20,20,1100,960])
% [ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% for i=1:length(list)
%     genePlot = list{i};
%     markergene = (data_orig_all(strcmpi(geneid_all,genePlot),:));
%     inpos = markergene>0;
%     tmpthlow = prctile(markergene(markergene>0),1);
%     tmpthhigh = prctile(markergene(markergene>0),99);
%     if tmpthlow==tmpthhigh
%         tmpthlow = 0;
%     end
%     markergene(markergene>tmpthhigh) = tmpthhigh;
%     markergene(markergene<tmpthlow) = tmpthlow;
%     c_rgb = [1,0,0];rand([1,3]);
%         markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%             ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
%     markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
%         interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
%         ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
%     axes(ha(i));
%     scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
%     scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
%     set(gca,'xlim',[-150,150],'ylim',[-150,150])
%     title(genePlot);
%     axis tight
%     axis equal
%     axis off
% end
% if savefig_flag==1
%     eval(['export_fig tsne_markers2_AmyPiri_FC_',date,'.pdf']);
%     savefig(gcf,['tsne_FC_Glut1_markers1_',date,'.fig'])
% end
% if savefig_pdf==1
%     eval(['export_fig tsne_FC_Glut1_markers1_',date,'.pdf']);
% end
% % % % % % % % % % % % %

meangr_mat1 = zeros(length(data(:,1)),length(T_cells_tmp_uni));
meangr_mat_all = zeros(length(data_orig_all(:,1)),length(T_cells_tmp_uni));
clust_cent = zeros(length(T_cells_tmp_uni),2);
prj_mean = zeros(initial_dims,length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    prj_mean(:,jjj) = mean(prj(T_cells_tmp(:,1)==T_cells_tmp_uni(jjj),:))';
    meangr_mat1(:,jjj) = mean(log2(data(:,T_cells_tmp(:,1)==T_cells_tmp_uni(jjj))+1),2);
    meangr_mat_all(:,jjj) = mean(log2(data_orig_all(:,T_cells_tmp(:,1)==T_cells_tmp_uni(jjj))+1),2);
%     clust_cent(jjj,:) = [median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),1)),median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),2))];
end
% save agg cluster
save(['all_neurons_cluster_ave'],'meangr_mat_all','clusteruni','geneid_all','geneid');

[prj,m,D,V,Q] = pca_wis(meangr_mat1',15);
% prj = prj_mean(1:10,:)';
% prj = meangr_mat1';

Zpca = linkage(prj,'ward','correlation');
Dpca = pdist(prj,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
figure;
set(gcf,'position',[100,100,1000,1000],'color','w')
h1 = axes('position',[0.03,0.03,0.2,0.93]);
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
h2 = axes('position',[0.37,0.03,0.63,0.93]);
x = squareform(Dpca); imagesc(x(leaforder_pca,leaforder_pca));
colormap('summer')
set(gca,'ytick',[1:length(clusteruni)],'YTickLabel',regexprep(clusteruni(leaforder_pca),'_','-'),'xtick',[],'fontsize',8,'ydir','normal')
linkaxes([h1,h2],'y');
if savefig_pdf==1
    eval(['export_fig tree_corrmat_neurons_per_cluster_',date,'.pdf']);
end

%%
init = prj(:,[1,2])/std(prj(:,1))*1e-4;
init = init-repmat(mean(init),length(init),1);
options = statset('MaxIter',1000,'Tolfun',1e-20);
perplexity = 30;
exagg = 1;
mapped_xy_clusters = tsne(prj,'Algorithm','exact','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',5,'Verbose',1,'Options',options,'Exaggeration',exagg);
% mapped_xy_clusters = prj(:,[2,3]);
% mapped_xy_clusters = tsne_vdM(prj, [], 2, length(prj(1,:)), perplexity, [500],2,1000);

% mapped_xy_clusters = mapped_xy_clusters-repmat(mean(mapped_xy_clusters),length(mapped_xy_clusters),1);
% mapped_xy_clusters = prj(:,[2,3]);

% colors = distinguishable_colors(length(T_cells_tmp_uni_l2));
colors = distinguishable_colors(3);
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(T_cells_tmp_uni)
    scatter(mapped_xy_clusters(i,1),mapped_xy_clusters(i,2),clustersize(i),hex2rgb(clusteruni{i,3}),'filled','o','MarkerEdgeColor','k'); hold on;
end
% for i=1:length(T_cells_tmp_uni)
% %     colorin = unique(T_cells_tmp(T_cells_tmp_uni(i)==T_cells_tmp(:,1),2));
% %     ind = unique(T_cells_tmp(find(T_cells_tmp_uni_l2(i)==T_cells_tmp(:,2))));
%     if contains(clusteruni{i},'GABA')
% %         plot(mapped_xy_clusters(i,1),mapped_xy_clusters(i,2),'.','color',colors(1,:),'markersize',10); hold on;
%         scatter(mapped_xy_clusters(i,1),mapped_xy_clusters(i,2),clustersize(i),colors(1,:),'filled','o','MarkerEdgeColor','k'); hold on;
%     elseif contains(clusteruni{i},'Glut1')
% %         plot(mapped_xy_clusters(i,1),mapped_xy_clusters(i,2),'.','color',colors(2,:),'markersize',10); hold on;
%         scatter(mapped_xy_clusters(i,1),mapped_xy_clusters(i,2),clustersize(i),colors(2,:),'filled','o','MarkerEdgeColor','k'); hold on;
%     elseif contains(clusteruni{i},'Glut2')
% %         plot(mapped_xy_clusters(i,1),mapped_xy_clusters(i,2),'.','color',colors(3,:),'markersize',10); hold on;
%         scatter(mapped_xy_clusters(i,1),mapped_xy_clusters(i,2),clustersize(i),colors(3,:),'filled','o','MarkerEdgeColor','k'); hold on;
%     end
% end
axis tight;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
rr = 0.05*diff(xl);
for i=1:length(T_cells_tmp_uni)
    if i>=0
        in = i;
        th = atan(mapped_xy_clusters(in,2)/mapped_xy_clusters(in,1)); %(rand*360)*pi/180;
        xo = rr*cos(th); yo = rr*sin(th);
        plot([mapped_xy_clusters(in,1),mapped_xy_clusters(in,1)+xo] , ...
            [mapped_xy_clusters(in,2),mapped_xy_clusters(in,2)+yo],'-','color',0.7*[1,1,1]);
        ht = text(mapped_xy_clusters(in,1)+xo,mapped_xy_clusters(in,2)+yo,regexprep(clusteruni{i,2},'_','-'));
        set(ht,'fontsize',8,'FontName','arial');%,'BackgroundColor',0.8*[1,1,1]
    end
end
% for i=1:length(T_cells_tmp_uni)
% %     ind = unique(T_cells_tmp(find(T_cells_tmp_uni_l2(i)==T_cells_tmp(:,2))));
%     ht = text(mapped_xy_clusters(i,1)-0.02*diff(xl),mapped_xy_clusters(i,2)-0.02*diff(yl),regexprep(clusteruni(ind),'_','-'));
%     set(ht,'fontsize',6);%,'BackgroundColor',0.8*[1,1,1]
% end
axis tight;
axis equal
axis off
title(['perplexity=',num2str(perplexity),',Exagg=',num2str(exagg),',#C=',num2str(length(T_cells_tmp_uni)),',#PC=',num2str(size(prj,2))])
if savefig_pdf==1
    eval(['export_fig tsne_neurons_per_cluster_',date,'.pdf']);
end
table1 = [clusteruni,m2c(clustersize),m2c(mapped_xy_clusters)];
tabel1 = [{'oldname','newname','hexcol','size','x','y'};table1];
saveCellFile(table1,['celltype_tsne_all_xy_size_colors_',date,'.txt'])
%%
% % % % % % % % % % % % % % % % % % % % % % 
% in = find(sum(data_orig_all>0,2)>5 & sum(data_orig_all>0,2)<length(data_orig_all(1,:))*0.5 & ~ismember(geneid_all(:,1),sex_genes)  ....
%     & ~ismember(geneid_all(:,1),[IEG_markers]) & ~ismember(geneid_all(:,1),[not_relevant]));
% top_g = 150;
% gr1 = find(gaba_glut_flag==2);
% gr2 = find(gaba_glut_flag==3);
% % ptt = ones(size(in));
% % for s=1:length(in)
% %  [~,ptt(s)] = kstest2(data_orig_all_sorted(in(s),gr1),data_orig_all_sorted(in(s),gr2),'tail','larger');
% % end
% x1 = mean(log2(data_orig_all(:,gr1)+1),2);
% x2 = mean(log2(data_orig_all(:,gr2)+1),2);
% d = x1-x2 ;
% [~,xi] = sort(d);
% figure('position',[200,200,1000,580],'color','w');
% [ha, pos] = tight_subplot(1, 2, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
% axes(ha(1))
% plot(x1, x2, '.');hold on;
% xmax = max(x1);
% plot([0,xmax],[0,xmax],'-k'); grid on
% plot([0,xmax],[0,xmax]+1,'--k'); grid on
% plot([1,xmax],[0,xmax-1],'--k'); grid on
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% xi = flipud(xi);
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% % title(['Sex comparison (average)'])
% xlabel(['mean (glut1)'])
% ylabel(['mean (glut2)'])
% axis tight
% 
% top_g = 250;
% x1 = mean(data_orig_all(:,gr1)>0,2);
% x2 = mean(data_orig_all(:,gr2)>0,2);
% d = x1-x2;
% d(x1>0.5 & x2>0.5) = 0;
% [~,xi] = sort(d);
% axes(ha(2))
% plot(x1, x2, '.');hold on;
% xmax = max(x1);
% plot([0,xmax],[0,xmax],'-k'); grid on
% plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
% plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% xi = flipud(xi);
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% % title(['Sex comparison (precentage)'])
% xlabel(['mean>0 (glut1)'])
% ylabel(['mean>0 (glut2)'])
% axis tight
% if savefig_pdf==1
%     eval(['export_fig scatter_glut1_vs_glut2_',date,'.pdf']);
% end
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% top_g = 150;
% gr1 = find(gaba_glut_flag==1);
% gr2 = find(gaba_glut_flag==2);
% % ptt = ones(size(in));
% % for s=1:length(in)
% %  [~,ptt(s)] = kstest2(data_orig_all_sorted(in(s),gr1),data_orig_all_sorted(in(s),gr2),'tail','larger');
% % end
% x1 = mean(log2(data_orig_all(:,gr1)+1),2);
% x2 = mean(log2(data_orig_all(:,gr2)+1),2);
% d = x1-x2 ;
% [~,xi] = sort(d);
% figure('position',[200,200,1000,580],'color','w');
% [ha, pos] = tight_subplot(1, 2, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
% axes(ha(1))
% plot(x1, x2, '.');hold on;
% xmax = max(x1);
% plot([0,xmax],[0,xmax],'-k'); grid on
% plot([0,xmax],[0,xmax]+1,'--k'); grid on
% plot([1,xmax],[0,xmax-1],'--k'); grid on
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% xi = flipud(xi);
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% % title(['Sex comparison (average)'])
% xlabel(['mean (gaba)'])
% ylabel(['mean (glut1)'])
% axis tight
% 
% top_g = 250;
% x1 = mean(data_orig_all(:,gr1)>0,2);
% x2 = mean(data_orig_all(:,gr2)>0,2);
% d = x1-x2;
% d(x1>0.5 & x2>0.5) = 0;
% [~,xi] = sort(d);
% axes(ha(2))
% plot(x1, x2, '.');hold on;
% xmax = max(x1);
% plot([0,xmax],[0,xmax],'-k'); grid on
% plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
% plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% xi = flipud(xi);
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% % title(['Sex comparison (precentage)'])
% xlabel(['mean>0 (gaba)'])
% ylabel(['mean>0 (glut1)'])
% axis tight
% if savefig_pdf==1
%     eval(['export_fig scatter_gaba_vs_glut1_',date,'.pdf']);
% end
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% top_g = 150;
% gr1 = find(gaba_glut_flag==1);
% gr2 = find(gaba_glut_flag==3);
% % ptt = ones(size(in));
% % for s=1:length(in)
% %  [~,ptt(s)] = kstest2(data_orig_all_sorted(in(s),gr1),data_orig_all_sorted(in(s),gr2),'tail','larger');
% % end
% x1 = mean(log2(data_orig_all(:,gr1)+1),2);
% x2 = mean(log2(data_orig_all(:,gr2)+1),2);
% d = x1-x2 ;
% [~,xi] = sort(d);
% figure('position',[200,200,1000,580],'color','w');
% [ha, pos] = tight_subplot(1, 2, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
% axes(ha(1))
% plot(x1, x2, '.');hold on;
% xmax = max(x1);
% plot([0,xmax],[0,xmax],'-k'); grid on
% plot([0,xmax],[0,xmax]+1,'--k'); grid on
% plot([1,xmax],[0,xmax-1],'--k'); grid on
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% xi = flipud(xi);
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% % title(['Sex comparison (average)'])
% xlabel(['mean (gaba)'])
% ylabel(['mean (glut2)'])
% axis tight
% 
% top_g = 250;
% x1 = mean(data_orig_all(:,gr1)>0,2);
% x2 = mean(data_orig_all(:,gr2)>0,2);
% d = x1-x2;
% d(x1>0.5 & x2>0.5) = 0;
% [~,xi] = sort(d);
% axes(ha(2))
% plot(x1, x2, '.');hold on;
% xmax = max(x1);
% plot([0,xmax],[0,xmax],'-k'); grid on
% plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
% plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% xi = flipud(xi);
% plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
% text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% % title(['Sex comparison (precentage)'])
% xlabel(['mean>0 (gaba)'])
% ylabel(['mean>0 (glut2)'])
% axis tight
% if savefig_pdf==1
%     eval(['export_fig scatter_gaba_vs_glut2_',date,'.pdf']);
% end
% 

toc

