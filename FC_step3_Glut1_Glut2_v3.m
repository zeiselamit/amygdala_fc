tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;


load afterloading_QC3000_FC_06-Oct-2021

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);

cellidclsuter_gaba = loadCellFile('cellidCluster_step2_glut1_FC_07-Oct-2021.txt');

cellidclsuter_glut1 = loadCellFile('cellidCluster_step2_glut1_FC_07-Oct-2021.txt');
cellidclsuter_glut1 = [cellidclsuter_glut1, cellfun(@(x,y) ['Glut1_',num2str(x),'_',num2str(y)],cellidclsuter_glut1(:,2),cellidclsuter_glut1(:,3),'UniformOutput',0)]; 
cellidclsuter_glut2 = loadCellFile('cellidCluster_step2_Glut2_FC_07-Oct-2021.txt');
cellidclsuter_glut2 = [cellidclsuter_glut2, cellfun(@(x,y) ['Glut2_',num2str(x),'_',num2str(y)],cellidclsuter_glut2(:,2),cellidclsuter_glut2(:,3),'UniformOutput',0)]; 

glut1_clustername = loadCellFile('FC_Glut1_cluster_list_211007');
glut1_clustername(:,2) = regexprep(glut1_clustername(:,2),'_glut1','');
glut1_clustername(:,2) = cellfun(@(x) ['Glut1_',x], glut1_clustername(:,2),'UniformOutput',0);
glut1_clustername(:,2) = regexprep(glut1_clustername(:,2),'Glut1_exc','exc');

glut2_clustername = loadCellFile('FC_Glut2_cluster_list_211007');
glut2_clustername(:,2) = regexprep(glut2_clustername(:,2),'_glut2','');
glut2_clustername(:,2) = cellfun(@(x) ['Glut2_',x], glut2_clustername(:,2),'UniformOutput',0);
glut2_clustername(:,2) = regexprep(glut2_clustername(:,2),'Glut2_exc','exc');



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
cellidclsuter_glut1 = [cellidclsuter_glut1;cellidclsuter_glut2];

cellid_clusters = cell(length(cellid),5);
[~,loc]= ismember(cellidclsuter_glut1(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_glut1;

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
not_relevant = {'Slc32a1','Gad1','Gad2','C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
    ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};%,'Snca','Gm9843','Gm12027','Gm12031'
in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),sex_genes)  ....
    & ~ismember(geneid(:,1),[IEG_markers]) & ~ismember(geneid(:,1),[not_relevant]));

corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[1000],1,0);
% g = unique([geneid(in(corr_filt));geneselection_gaba(:,1);]);
% [~,corr_filt] = ismember(g,geneid);

% % % % % % % % % % % 
z = linkage(data(in(corr_filt),:),'ward','correlation');
tic
idx = cluster(z,'maxclust',50);
toc
d = corr_mat(data(in(corr_filt),:)');
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
data = data(in(corr_filt),:);
geneid = geneid(in(corr_filt));

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
perplexity = 100;
options = statset('MaxIter',1000);
mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',5);
toc



% % % % % % % % % % % % % % % % % % % 
clusteruni = unique(cellid_clusters(:,5));
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,5),clusteruni{i}));
    T_cells_tmp(ind) = i;
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp);

% % % % % % % % % % % % % % % % % % % 
% colors = distinguishable_colors(length(T_cells_tmp_uni)+1);
tmp = loadCellFile_turbo('Glut1-hexa.txt',1);
tmp = tmp(2:end,2:3);
clustername_color_hannah = loadCellFile_turbo('Glut2-hexa.txt',1);
clustername_color_hannah = clustername_color_hannah(2:end,2:3);
clustername_color_hannah = [tmp;clustername_color_hannah];
[~,loc] = ismember(clusteruni,clustername_color_hannah(:,1));
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(T_cells_tmp_uni)
    ii = find(T_cells_tmp(:,1)==i); 
    h = plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',['#',clustername_color_hannah{loc(i),2}],'markersize',3); hold on;
end
axis tight;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
rr = 0.05*diff(xl);
for i=1:length(T_cells_tmp_uni)
    if i>=0
        in = T_cells_tmp==i;
        th = (rand*360)*pi/180;
        xo = rr*cos(th); yo = rr*sin(th);
        plot([median(mapped_xy(in,1)),median(mapped_xy(in,1))+xo] , ...
            [median(mapped_xy(in,2)),median(mapped_xy(in,2))+yo],'-','color',0.7*[1,1,1]);
        ht = text(median(mapped_xy(in,1))+xo,median(mapped_xy(in,2))+yo,regexprep(clusteruni{i},'_','-'));
        set(ht,'fontsize',8,'FontName','arial');%,'BackgroundColor',0.8*[1,1,1]
    end
end
axis tight;
axis equal
axis off
title(['#C=',num2str(max(T_cells_tmp)),', perplexity=',num2str(perplexity)],'fontsize',8);
if savefig_flag==1
    savefig(gcf,['tsne_final_step3_FC_Glut2_by_cluster_',date,'.fig'])    
end
if savefig_pdf==1
    eval(['export_fig tsne_final_step3_FC_Glut2_by_cluster_',date,'.pdf']);
end
% 
% sort the data by the clusters and remove outliers
[T_cells_tmp,xi] = sort(T_cells_tmp);
xi(T_cells_tmp==0) = [];
T_cells_tmp(T_cells_tmp==0) = [];
data_sorted_all = data(:,xi);
data_orig_all_sorted = data_orig_all(:,xi);
cellid_sorted = cellid((xi));
sample_sorted = sample((xi));
mapped_xy = mapped_xy(xi,:);
fc_time_sorted = fc_time(xi);
batch_sorted = batch_flag(xi);
amy_flag_sorted = amy_flag(xi);
piri_flag_sorted = piri_flag(xi);
prj_sorted = prj(xi,:);

no_dims = 1;
initial_dims = 10;
perplexity = 5;
epsilon = 100;
dist_flag = 2;
theta = 0.5;

rand_seed = 13;
data_tsne = cent_norm(log2(data_sorted_all+1));
xi = [1:length(T_cells_tmp)];
for i=1:length(T_cells_tmp_uni)
    i
    ind = find(T_cells_tmp==i);
    if length(ind)>20
        tmp1d = fast_tsne((data_tsne(:,ind))', no_dims, initial_dims, perplexity,theta, rand_seed);
        [~,xitmp] = sort(tmp1d);
        xi(ind) = xi(ind((xitmp)));
    end
end

data_sorted_all = data_sorted_all(:,xi);
data_orig_all_sorted = data_orig_all_sorted(:,xi);
cellid_sorted = cellid_sorted((xi));
sample_sorted = sample_sorted((xi));
mapped_xy = mapped_xy(xi,:);
fc_time_sorted = fc_time_sorted(xi);
amy_flag_sorted = amy_flag_sorted(xi);
batch_sorted = batch_sorted(xi);
piri_flag_sorted = piri_flag_sorted(xi);
prj_sorted = prj_sorted(xi,:);

meangr_mat = zeros(length(moldata(:,1)),length(T_cells_tmp_uni));
clust_cent = zeros(length(T_cells_tmp_uni),2);
for jjj=1:length(T_cells_tmp_uni)
    jjj
    meangr_mat(:,jjj) = mean(log2(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))+1),2);
    clust_cent(jjj,:) = [median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),1)),median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),2))];
end

meangr_mat1 = meangr_mat;
% % forcing the tree order from glu1 and glut2 separate analysis
glut1clusters = loadCellFile('ClusterOrder_step3_Glut1_FC_24-Oct-2021.txt');
glut2clusters = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');
[~,loc] = ismember([glut1clusters(:,2);glut2clusters(:,2)],clusteruni);
leaforder_pca = loc;%force the order of the glut1 glut2 separate
% % % 
[prj,m,D,V,Q] = pca_wis(meangr_mat1',initial_dims);
Zpca = linkage(prj,'ward','correlation');
Dpca = pdist(prj,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
figure;
set(gcf,'position',[100,100,1200,1000],'color','w')
h1 = axes('position',[0.03,0.03,0.21,0.93]);
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
h2 = axes('position',[0.35,0.03,0.6,0.93]);
x = corr_mat(prj'); imagesc(x(leaforder_pca,leaforder_pca));
cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);
colormap((cm))
set(gca,'ytick',[1:length(clusteruni)],'YTickLabel',regexprep(clusteruni(leaforder_pca),'_','-'),'xtick',[],'fontsize',8,'ydir','normal')
cb = colorbar('Position',[0.96,0.03,0.02,0.93]);
linkaxes([h1,h2],'y');
if savefig_flag==1
    savefig(gcf,['tree_final_step3_FC_Glut2_',date,'.fig'])    
end
if savefig_pdf==1
    eval(['export_fig tree_final_step3_FC_Glut2_',date,'.pdf']);
end
leaforder = leaforder_pca;
Zpca_post = linkage(prj(leaforder_pca,:),'ward','correlation');

bigclust = cluster(Zpca_post,'maxclust',10);
T_cells_tmp_new = zeros(length(T_cells_tmp),2);
for i=1:length(leaforder)
    T_cells_tmp_new(T_cells_tmp==T_cells_tmp_uni(leaforder(i)),1) = i;
    T_cells_tmp_new(T_cells_tmp==T_cells_tmp_uni(leaforder(i)),2) = bigclust(i);
end

idxuni_new = unique(T_cells_tmp_new);
[~,xi] = sort(T_cells_tmp_new(:,1));
T_cells_tmp_new = T_cells_tmp_new(xi,:);
data_sorted_all = data_sorted_all(:,xi);
data_orig_all_sorted = data_orig_all_sorted(:,xi);
cellid_sorted = cellid_sorted(xi);
mapped_xy = mapped_xy(xi,:);
fc_time_sorted = fc_time_sorted(xi);
amy_flag_sorted = amy_flag_sorted(xi);
piri_flag_sorted = piri_flag_sorted(xi);
cells_bor_2 = find(diff(T_cells_tmp_new)>0)+1;
sample_sorted = sample_sorted(xi);
batch_sorted = batch_sorted(xi);
prj_sorted = prj_sorted(xi,:);

tmp = zeros(size(T_cells_tmp_new(:,2)));
tmp(1) = 1;
for i=2:length(T_cells_tmp_new)
    if T_cells_tmp_new(i,2) == T_cells_tmp_new(i-1,2)
        tmp(i) = tmp(i-1);
    else
        tmp(i) = tmp(i-1)+1;
    end
end
T_cells_tmp_new(:,2) = tmp;
% % % % % % % % % % % % % % 
T_cells_tmp = T_cells_tmp_new;
T_cells_tmp_uni = unique(T_cells_tmp);
% % % % % % % % % % % % % % 
%%
glutspatial = loadCellFile('GLUT12_spatialclusters.txt');
T_cells_tmp_spatial = zeros(length(T_cells_tmp),1);
for i=1:length(T_cells_tmp_uni)
    ii=find(T_cells_tmp(:,1)==i);
    T_cells_tmp_spatial(ii) = glutspatial{i,4};
end
T_cells_tmp_spatial_uni = unique(T_cells_tmp_spatial);
discolors = distinguishable_colors(length(T_cells_tmp_spatial_uni));
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(T_cells_tmp_spatial_uni)
    ii=find(T_cells_tmp_spatial==T_cells_tmp_spatial_uni(i)); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',discolors(i,:),'markersize',3); hold on;
end
axis tight;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
rr = 0.05*diff(xl);
for i=1:length(T_cells_tmp_uni)
    if i>=0
        in = T_cells_tmp(:,1)==i;
        th = (rand*360)*pi/180;
        xo = rr*cos(th); yo = rr*sin(th);
        plot([median(mapped_xy(in,1)),median(mapped_xy(in,1))+xo] , ...
            [median(mapped_xy(in,2)),median(mapped_xy(in,2))+yo],'-','color',0.7*[1,1,1]);
        ht = text(median(mapped_xy(in,1))+xo,median(mapped_xy(in,2))+yo,regexprep(clusteruni{leaforder(i)},'_','-'));
        set(ht,'fontsize',8,'FontName','arial');%,'BackgroundColor',0.8*[1,1,1]
    end
end
axis tight;
axis equal
axis off
title(['#C=',num2str(max(T_cells_tmp)),', perplexity=',num2str(perplexity)],'fontsize',8);
[~,ia] = unique(cell2mat(glutspatial(:,4)));
legend(glutspatial(ia,5))
if savefig_pdf==1
    eval(['export_fig tsne_final_step3_FC_Glut12_by_spatialcluster_',date,'.pdf']);
end
%%
% % % % % % % % % % % % % % 
[table1, table2] = dendrogram_split_markers(cellfun(@(x) num2str(x), m2c(T_cells_tmp_uni),'UniformOutput',0).....
    ,cellfun(@(x) num2str(x), m2c(T_cells_tmp(:,1)),'UniformOutput',0),Zpca_post,data_sorted_all,geneid);
saveCellFile(table1,['FC_final_step3_Glut2_dendrogram_junction_split_markers_',date,'.txt']);
saveCellFile(table2,['FC_final_step3_Glut2_dendrogram_junction_split_markers_by_average_',date,'.txt']);
ind = [[3:53:length(table1(:,1))];[4:53:length(table1(:,1))];[5:53:length(table1(:,1))];[6:53:length(table1(:,1))];[7:53:length(table1(:,1))]];
ind = ind(:);
treemark1 = table1(ind,[1,7]);
treemark1 = [reshape(treemark1(:,1),5,[]);reshape(treemark1(:,2),5,[])];
treemark1 = flipud(treemark1(:));
rmv = [];
for i=2:length(treemark1)
    if sum(strcmpi(treemark1{i},treemark1(1:i-1)))>0
        rmv = [rmv,i];
    end
end
treemark1(rmv) = [];


% % % % % % % % % % % % % % % % % % 
% plot by sample
% sample_uni = {'8-1','8-2','10-1','10-2','18-1','18-2','19-1','23-1','23-2','23-3','41-1','41-2','42-1','45-1','45-2','45-3','46-1','46-2','46-3'};
colors = distinguishable_colors(length(sample_uni)+1);
figure;
set(gcf,'color','w','position',[20,20,1200,800]);
[ha, pos] = tight_subplot(1, 2, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
axes(ha(1))
for i=1:length(sample_uni)
    s = scatter(mapped_xy(strcmpi(sample_sorted,sample_uni{i}),1),mapped_xy(strcmpi(sample_sorted,sample_uni{i}),2),10,colors(i,:),'filled'); hold on;
%     alpha(s,0.4);
end
axis tight
axis equal
axis off
legend(sample_uni)
fc_time_uni = [0,2,8,24,28];
colors = distinguishable_colors(length(fc_time_uni)+1);
axes(ha(2))
for i=1:length(fc_time_uni)
    s = scatter(mapped_xy(fc_time_sorted==fc_time_uni(i),1),mapped_xy(fc_time_sorted==fc_time_uni(i),2),10,colors(i,:),'filled'); hold on;
%     alpha(s,0.3);
end
axis tight
axis equal
axis off
legend('ctrl','2hr','8hr','24hr','recall')
% plot by tissue
if savefig_flag==1
    eval(['export_fig tsne_final_step3_FC_Glut2_by_annot_','_',date,'.pdf']);
    savefig(gcf,['tsne_final_step3_FC_Glut2_by_annot_perplexity_',date,'.fig'])
end

figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(4, 6, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
xy_range = [min(mapped_xy),max(mapped_xy)];
for i=1:length(sample_uni)
    axes(ha(i));
    plot(mapped_xy(strcmpi(sample_sorted,sample_uni{i}),1),mapped_xy(strcmpi(sample_sorted,sample_uni{i}),2),'.');
    title(sample_uni{i})
    set(gca,'xlim',[xy_range(1),xy_range(3)],'ylim',[xy_range(2),xy_range(4)])
    axis off
    axis equal
    %     alpha(s,0.4);
end

figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(2, 3, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
xy_range = [min(mapped_xy),max(mapped_xy)];
for i=1:length(fc_time_uni)
    axes(ha(i));
    plot(mapped_xy(fc_time_sorted==fc_time_uni(i),1),mapped_xy(fc_time_sorted==fc_time_uni(i),2),'.','markersize',3);
    title(num2str(fc_time_uni(i)))
    set(gca,'xlim',[xy_range(1),xy_range(3)],'ylim',[xy_range(2),xy_range(4)])
    axis off
    axis equal
    %     alpha(s,0.4);
end
if savefig_flag==1
    eval(['export_fig tsne_final_step3_FC_Glut2_bytimepoint_',date,'.pdf']);
    savefig(gcf,['tsne_final_step3_FC_Glut2_bytimepoint_',date,'.fig'])
end
% % 
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(2, 3, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
cmap = 0.5+0.5*flipud(colormap('gray'));
axes(ha(1))
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.8*[1,1,1],'MarkerSize',4); hold on;
xymin = 1.1*min(mapped_xy);
xymax = 1.1*max(mapped_xy);
n = 2^6;
[bandwidth,density_all,X,Y] = kde2d(mapped_xy,n,xymin,xymax);
z = density_all; z(z<0) = 0;
contour(X,Y,z,'r','LineWidth',1); hold on;axis off;axis equal
title('all')
axes(ha(2))
gr1 = fc_time_sorted==0;
plot(mapped_xy(gr1,1),mapped_xy(gr1,2),'.','color',0.8*[1,1,1],'MarkerSize',4); hold on;
[bandwidth,density,X,Y] = kde2d(mapped_xy(gr1,:),n,xymin,xymax);
z = density-density_all; z(z<0) = 0;
contour(X,Y,z,'r','LineWidth',1); hold on;axis off;axis equal
title('0hr')
axes(ha(3))
gr1 = fc_time_sorted==2;
plot(mapped_xy(gr1,1),mapped_xy(gr1,2),'.','color',0.8*[1,1,1],'MarkerSize',4); hold on;
[bandwidth,density,X,Y] = kde2d(mapped_xy(gr1,:),n,xymin,xymax);
z = density-density_all; z(z<0) = 0;
contour(X,Y,z,'r','LineWidth',1); hold on;axis off;axis equal
title('2hr')
axes(ha(4))
gr1 = fc_time_sorted==8;
plot(mapped_xy(gr1,1),mapped_xy(gr1,2),'.','color',0.8*[1,1,1],'MarkerSize',4); hold on;
[bandwidth,density,X,Y] = kde2d(mapped_xy(gr1,:),n,xymin,xymax);
z = density-density_all; z(z<0) = 0;
contour(X,Y,z,'r','LineWidth',1); hold on;axis off;axis equal
title('8hr')
axes(ha(5))
gr1 = fc_time_sorted==24;
plot(mapped_xy(gr1,1),mapped_xy(gr1,2),'.','color',0.8*[1,1,1],'MarkerSize',4); hold on;
[bandwidth,density,X,Y] = kde2d(mapped_xy(gr1,:),n,xymin,xymax);
z = density-density_all; z(z<0) = 0;
contour(X,Y,z,'r','LineWidth',1); hold on;axis off;axis equal
title('24hr')
axes(ha(6))
gr1 = fc_time_sorted==28;
plot(mapped_xy(gr1,1),mapped_xy(gr1,2),'.','color',0.8*[1,1,1],'MarkerSize',4); hold on;
[bandwidth,density,X,Y] = kde2d(mapped_xy(gr1,:),n,xymin,xymax);
z = density-density_all; z(z<0) = 0;
contour(X,Y,z,'r','LineWidth',1); hold on;axis off;axis equal
title('recall')
% if savefig_flag==1
%     eval(['export_fig tsne_final_FC_Glut2_bytimepoint_','_',date,'.pdf']);
%     savefig(gcf,['tsne_final_FC_Glut2_bytimepoint_',num2str(perplexity),'_PCAdim',num2str(initial_dims),'_',date,'.fig'])
% end
% % % % % % % % % % % % % % % % % % 
[ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp(:,1),data_sorted_all,5);
% 
cells_bor22 = find(diff(T_cells_tmp(:,2))>0);
datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
% gr_tmp_mark = gr_tmp_mark(xi);
gr_tmp_mark = geneid(ind_gr_tmp_mark);

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.02,0.88,0.87]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor)
    plot(cells_bor(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
for jj=1:length(cells_bor22)
    plot(cells_bor22(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',0.7*[1,0,0])
end
set(gca,'xtick',gr_center,'xticklabel',[1:length(gr_center)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)
colormap('summer');
freezeColors(gca);
text([cells_bor;cells_bor(end)],linspace(5,length(gr_tmp_mark),length(gr_center)),regexprep(clusteruni(leaforder_pca),'_','-'),'color','k')

% sample_uni = {'8-1','8-2','10-1','10-2','18-1','18-2','19-1','23-1','23-2','23-3','41-1','41-2','42-1','42-2','45-1','45-2','45-3','46-1','46-2','46-3'};
samples_num = false(length(sample_uni),length(sample_sorted));
for i=1:length(sample_uni)
    samples_num(i, strcmpi(sample_sorted,sample_uni{i})) = true;
end

ax2_1 = axes('position',[0.1,0.9,0.88,0.012]);
imagesc(~(fc_time_sorted'==0)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','0');

ax2_2 = axes('position',[0.1,0.9+1*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==2)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','2');

ax2_3 = axes('position',[0.1,0.9+2*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==8)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','8');

ax2_4 = axes('position',[0.1,0.9+3*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==24)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','24');

ax2_5 = axes('position',[0.1,0.9+4*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==28)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','recall');

% ax3 = axes('position',[0.1,0.81,0.88,0.01]);
% imagesc(~amy_flag_sorted'); hold on;
% colormap('gray');
% freezeColors(gca);
% set(gca,'xtick',[],'ytick',[1],'yticklabel','Amy');
% 
% ax4 = axes('position',[0.1,0.82,0.88,0.01]);
% imagesc(~piri_flag_sorted'); hold on;
% colormap('gray');
% freezeColors(gca);
% set(gca,'xtick',[],'ytick',[1],'yticklabel','Piri');

% ax5 = axes('position',[0.1,0.83,0.88,0.12]);
% imagesc(~samples_num); hold on;
% grid on
% colormap('gray');
% freezeColors(gca);
% set(gca,'xtick',[],'ytick',[1:length(sample_uni)],'yticklabel',sample_uni,'fontsize',6);
% 
% gad2 = data_orig_all_sorted(strcmpi(geneid_all,'Gad2'),:);
% slc17a7 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a7'),:);
% slc17a6 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a6'),:);
% 
% ax6 = axes('position',[0.1,0.95,0.88,0.01]);
% imagesc(~gad2); hold on;
% axes(ax6)
% colormap('gray');
% freezeColors(ax6);
% set(ax6,'xtick',[],'ytick',[1],'yticklabel','Gad2','fontsize',5);
% ax7 = axes('position',[0.1,0.95+0.01,0.88,0.01]);
% imagesc(~slc17a6); hold on;
% colormap('gray');
% freezeColors(ax7);
% set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a6','fontsize',6);
% ax8 = axes('position',[0.1,0.95+2*0.01,0.88,0.01]);
% imagesc(~slc17a7); hold on;
% colormap('gray');
% freezeColors(ax8);
% set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a7','fontsize',6);

linkaxes([ax1,ax2_1,ax2_2,ax2_3,ax2_4,ax2_5],'x');
if savefig_flag==1
    savefig(gcf,['markertable_step3_final_FC_Glut2_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
if savefig_pdf==1
    eval(['export_fig markertable_step3_final_FC_Glut2_',date,'.pdf']);
end

table1 = [cellid_sorted,m2c(T_cells_tmp)];
saveCellFile(table1,['cellidCluster_step3_Glut2_FC_',date,'.txt']);
saveCellFile([m2c([1:length(leaforder_pca)]'),clusteruni(leaforder_pca)],['ClusterOrder_step3_Glut2_FC_',date,'.txt']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% [ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp,data_sorted_all,5);
% 
[~,ind_gr_tmp_mark] = ismember(treemark1,geneid);
datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
% gr_tmp_mark = gr_tmp_mark(xi);
gr_tmp_mark = geneid(ind_gr_tmp_mark);

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.02,0.88,0.87]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor)
    plot(cells_bor(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',gr_center,'xticklabel',[1:length(gr_center)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)
colormap('summer');
freezeColors(gca);

% sample_uni = {'8-1','8-2','10-1','10-2','18-1','18-2','19-1','23-1','23-2','23-3','41-1','41-2','42-1','42-2','45-1','45-2','45-3','46-1','46-2','46-3'};
samples_num = false(length(sample_uni),length(sample_sorted));
for i=1:length(sample_uni)
    samples_num(i, strcmpi(sample_sorted,sample_uni{i})) = true;
end

ax2_1 = axes('position',[0.1,0.9,0.88,0.012]);
imagesc(~(fc_time_sorted'==0)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','0');

ax2_2 = axes('position',[0.1,0.9+1*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==2)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','2');

ax2_3 = axes('position',[0.1,0.9+2*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==8)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','8');

ax2_4 = axes('position',[0.1,0.9+3*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==24)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','24');

ax2_5 = axes('position',[0.1,0.9+4*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==28)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','recall');

linkaxes([ax1,ax2_1,ax2_2,ax2_3,ax2_4,ax2_5],'x');
if savefig_flag==1
    savefig(gcf,['junctionmarkers_final_step3_FC_Glut2_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
% % % % % % % % % 
meangr_mat_0 = zeros(length(moldata(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    meangr_mat_0(:,jjj) = mean(data_sorted_all(:,T_cells_tmp(:,1)==T_cells_tmp_uni(jjj))>0,2);
end

datamarkers = meangr_mat_0(ind_gr_tmp_mark,:);
datamarkers_cn = datamarkers;%cent_norm(datamarkers);(lo,:)

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.14,0.88,0.65]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
set(gca,'xtick',[1:length(gr_center)],'xticklabel',regexprep(clusteruni(leaforder_pca),'_','-')....
    ,'XTickLabelRotation',45,'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)
colormap('summer');
freezeColors(gca);

h2 = axes('position',[0.1,0.81,0.88,0.2]);
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','top');
axis off
set(gca,'xlim',[0.5,length(leaforder_pca)+0.5])
linkaxes([ax1,h2],'x')

if savefig_flag==1
    savefig(gcf,['junctionmarkers_final_step3_clusteraverage_FC_Glut2_',date,'.fig'])
    eval(['export_fig junctionmarkers_final_step3_clusteraverage_FC_Glut2_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
meangr_mat_0 = zeros(length(geneid_all),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    meangr_mat_0(:,jjj) = mean(data_orig_all_sorted(:,T_cells_tmp(:,1)==T_cells_tmp_uni(jjj))>0,2);
end
% % % % % % % % % % % % % % % % % % 
list = {'Calcrl';'Prkcd';'Sst';'Crh';'Tac2';'Nts';'Ppp1r1b';'Rspo2';'Drd1';'Drd2';'Adora2a';'Penk';'Tac1';'Pdyn'};
[~,ind_list] = ismember(list,geneid_all);
datamarkers = meangr_mat_0(ind_list,:);

datamarkers_cn = datamarkers;%cent_norm(datamarkers);(lo,:)

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.02,0.88,0.75]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
set(gca,'xtick',[1:length(gr_center)],'xticklabel',[1:length(gr_center)],'ytick',[1:length(list)],'yticklabel',list, 'fontsize', 10)
colormap('summer');
freezeColors(gca);

h2 = axes('position',[0.1,0.79,0.88,0.21]);
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','top');
axis off
set(gca,'xlim',[0.5,length(leaforder_pca)+0.5])
linkaxes([ax1,h2],'x')
if savefig_flag==1
    savefig(gcf,['CEA_markers_step3_clusteraverage_FC_Glut2_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
list = loadCellFile('N_Z_genelist.txt');
[~,ind_list] = ismember(list,geneid_all);
list(ind_list==0) = [];
ind_list(ind_list==0) = [];
datamarkers = meangr_mat_0(ind_list,:);

datamarkers_cn = datamarkers;%cent_norm(datamarkers);(lo,:)

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.02,0.88,0.75]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
set(gca,'xtick',[1:length(gr_center)],'xticklabel',[1:length(gr_center)],'ytick',[1:length(list)],'yticklabel',list, 'fontsize', 10)
colormap('summer');
freezeColors(gca);

h2 = axes('position',[0.1,0.79,0.88,0.21]);
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','top');
axis off
set(gca,'xlim',[0.5,length(leaforder_pca)+0.5])
linkaxes([ax1,h2],'x')
if savefig_flag==1
    savefig(gcf,['N_Z_list_step3_clusteraverage_FC_Glut2_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end

% % % % % % % % % % % % % % % % % % % % % % % % % 

[ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp(:,1),data_sorted_all,1);
list =  geneid(ind_gr_tmp_mark(1:end));
% list = {'Meis2','Maf','Zeb2','Prlr','Lhx6','Npas1','Igfbp4','Sst','Cryab','Penk','Moxd1','Prkcd','Unc5b','Ucn3','Pthlh','Vip','Chat','Htr3a','Vdr'};
% list = {'Fth1','Ftl1','Scrib','Tfr2','Tfrc','Slc11a2','Slc40a1','Aco1','Ireb2'};
figure;
set(gcf,'color','w','position',[20,20,1000,1000]);
[ha, pos] = tight_subplot(length(list), 1, [0.00,0.00], [0.01,0.01], [0.09,0.01]);
logflag = 0;
for jjj=1:length(list)
    gn = list{jjj};
    g = find(strcmpi(gn,geneid_all));
    p99 = prctile(data_orig_all_sorted(g,:),100);
    axes(ha(jjj))
    p90 = max(ceil(prctile(data_orig_all_sorted(g,:),90)),1);
    t_ed = zeros(length(T_cells_tmp_uni),1);
    t_av = zeros(length(T_cells_tmp_uni),1);
    t_75 = zeros(length(T_cells_tmp_uni),1);
    for k=1:length(T_cells_tmp_uni)
        c=k;
        gr2 = find(T_cells_tmp(:,1)==c);
        y = (data_orig_all_sorted(g,gr2));
        if logflag ==1
            y = y+1;
            yscale = 'log';
        else
            yscale = 'linear';
        end
       if length(y)>10
            [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
            fi = interp1(xi,f,y);
            fi = fi/max(fi);
            plot(k + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y'....
                ,'.','color',['#',clustername_color_hannah{k,2}],'markersize',4); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)'....
                ,'.','color',['#',clustername_color_hannah{k,2}],'markersize',4); hold on;
        end
        t_ed(k) = median(data_orig_all_sorted(g,gr2));
        t_av(k) = mean(data_orig_all_sorted(g,gr2));
        t_75(k) = prctile(data_orig_all_sorted(g,gr2),75);
    end
    axis tight
    yl = get(gca,'ylim');
    if jjj==length(list)
        set(gca,'xtick',[1:length(T_cells_tmp_uni)],'ylim',[-0.5,p99+1],'YScale',yscale);
    else
        set(gca,'xtick',[1:length(T_cells_tmp_uni)],'xticklabel',cell(size(T_cells_tmp_uni)),'ylim',[-0.5,p99+1],'YScale',yscale)
    end
    text(-3,p99/2,[gn])
    box off
end

if savefig_flag==1
    savefig(gcf,['clustermarkers_violin_final_step3_FC_Glut2_',date,'.fig'])
    eval(['export_fig clustermarkers_violin_final_step3_FC_Glut2_',date,'.pdf']);
end
% % % % % % % % % % % % % % 
[ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp(:,1),data_sorted_all,1);
list =  geneid(ind_gr_tmp_mark(1:end));
[~,ind_gr_tmp_mark] = ismember(list,geneid_all);

datamarkers = meangr_mat_0(ind_gr_tmp_mark,:);
datamarkers_cn = datamarkers;%cent_norm(datamarkers);(lo,:)

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.14,0.88,0.65]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
set(gca,'xtick',[1:length(gr_center)],'xticklabel',regexprep(clusteruni(leaforder_pca),'_','-')....
    ,'XTickLabelRotation',45,'ytick',[1:length(list)],'yticklabel',list, 'fontsize', 10)
colormap('summer');
freezeColors(gca);

h2 = axes('position',[0.1,0.81,0.88,0.2]);
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','top');
axis off
set(gca,'xlim',[0.5,length(leaforder_pca)+0.5])
linkaxes([ax1,h2],'x')

if savefig_flag==1
    savefig(gcf,['clustermarkers_final_step3_clusteraverage_FC_Glut2_',date,'.fig'])
    eval(['export_fig clustermarkers_final_step3_clusteraverage_FC_Glut2_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % % % 
% gn = 'scrib';
% g = find(strcmpi(gn,geneid_all));
% figure;
% set(gcf,'color','w','position',[20,20,1200,1000]);
% b = 1;
% if b==2;
%     [ha, pos] = tight_subplot(4, ceil(length(unique(T_cells_tmp(:,b)))/4), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% elseif b==1;
%     [ha, pos] = tight_subplot(8, ceil(length(unique(T_cells_tmp(:,b)))/8), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% end
% logflag = 0;
% p99 = prctile(data_orig_all_sorted(g,:),99.99);
% for k=1:length(unique(T_cells_tmp(:,b)))
%     axes(ha(k))
%     c=k;
%     t_ed = zeros(5,1);
%     t_av = zeros(5,1);
%     t_75 = zeros(5,1);
%     p90 = max(ceil(prctile(data_orig_all_sorted(g,T_cells_tmp(:,b)==c),90)),1);
%     for i=[1:5]
%         gr2 = find(T_cells_tmp(:,b)==c & fc_time_sorted==fc_time_uni(i));%find(fc_time_sorted==fc_time_uni(i));%
%         y = (data_orig_all_sorted(g,gr2));
%         if logflag ==1
%             y = y+1;
%             yscale = 'log';
%         else
%             yscale = 'linear';
%         end
%         if length(y)>10
%             [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
%             fi = interp1(xi,f,y);
%             fi = fi/max(fi);
%             plot(i + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.'); hold on;
%         else
%             plot(i+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)','.'); hold on;
%         end
%         
%         %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)','.'); hold on;
%         t_ed(i) = median(data_orig_all_sorted(g,gr2));
%         t_av(i) = mean(data_orig_all_sorted(g,gr2));
%         t_75(i) = prctile(data_orig_all_sorted(g,gr2),75);
%         f90(i) = sum(data_orig_all_sorted(g,gr2)>=p90);
%         pp90(i) = binocdf(f90(i),length(gr2),0.1,'upper');
%     end
%     plot([0,6],p90*[1,1],'--k');
%     %     plot(t_ed,'-sk'); hold on;
%     %     plot(t_av,'-or'); hold on;
%     %     plot(t_75,'-dg'); hold on;
%     axis tight
%     yl = get(gca,'ylim');
%     set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear','ylim',[-0.5,p99+1],'YScale',yscale)
%     text([1:5],yl(2)*0.9*ones(1,5),cellfun(@(x) num2str(x,2), m2c(-log10(pp90)),'UniformOutput',0),'fontsize',6);
%     title([gn,',c=',num2str(c)])
% end
%%
tmp = loadCellFile_turbo('Glut1-hexa.txt',1);
tmp = tmp(2:end,2:3);
clustername_color_hannah = loadCellFile_turbo('Glut2-hexa.txt',1);
clustername_color_hannah = clustername_color_hannah(2:end,2:3);
clustername_color_hannah = [tmp;clustername_color_hannah];
% 
% clustername_color_hannah = loadCellFile_turbo('Glut2-hexa.txt',1);
% clustername_color_hannah = clustername_color_hannah(2:end,2:3);
p90_mat = zeros(length(geneid_all),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    p90_mat(:,jjj) = prctile(data_orig_all_sorted(:,T_cells_tmp(:,1)==T_cells_tmp_uni(jjj)),95,2);
end
[~,loc] = ismember(clustername_color_hannah(:,1),clusteruni(leaforder_pca));
list = loadCellFile('Glut2-markers.txt');
list = list(2:end,:);
list = unique([list;treemark1]);
logflag = 0;
for jjj=1:length(list)
    jjj
    gn = list{jjj}
    g = find(strcmpi(gn,geneid_all));
    hf1 = figure;
    set(gcf,'color','w','position',[200,200,1312,30],'Visible','off');
    
    p99 = max(p90_mat(g,:));
    axes('Position',[0.1,0.1,0.85,0.9])
    p90 = max(ceil(prctile(data_orig_all_sorted(g,:),90)),1);
    t_ed = zeros(length(T_cells_tmp_uni),1);
    t_av = zeros(length(T_cells_tmp_uni),1);
    t_75 = zeros(length(T_cells_tmp_uni),1);
    for k=1:length(T_cells_tmp_uni)
        c=k;
        gr2 = find(T_cells_tmp(:,1)==c);
        y = (data_orig_all_sorted(g,gr2));
        y(y>p99) = p99;
        if logflag ==1
            y = y+1;
            yscale = 'log';
        else
            yscale = 'linear';
        end
        if length(y)>10
            [f,xi] = ksdensity(y);
%             [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
            fi = interp1(xi,f,y);
            fi = fi/max(fi);
            plot(k + fi'.*(0.7*rand(length(gr2),1)-0.35), 0.5*rand(length(gr2),1)-0.1+y',...
                '.','color',['#',clustername_color_hannah{k,2}],'markersize',4); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)'...
                ,'.','color',['#',clustername_color_hannah{k,2}],'markersize',4); hold on;
        end
        t_ed(k) = median(data_orig_all_sorted(g,gr2));
        t_av(k) = mean(data_orig_all_sorted(g,gr2));
        t_75(k) = prctile(data_orig_all_sorted(g,gr2),75);
    end
    axis tight
    plot(t_av,'sk','markersize',3,'markerfacecolor','k')
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
    set(gca,'xtick',[1:length(T_cells_tmp_uni)],'xticklabel',cell(size(T_cells_tmp_uni)),'ylim',[-0.5,p99+1],'YScale',yscale,'ytick',yt,'YTickLabel',cell(size(yt)));
    %     text(-3,p99/2,[gn])
    box off
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png(['glut12_singlegene_violin_220222/',gn,'_max',num2str(round(yt(end))),'_Glut12_violinSingle_',date],gcf,1000)
    close(hf1)
end
%%

list = loadCellFile('Glut2-markers.txt');
list = list(2:end,:);
list = unique([list;treemark1]);
% [ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for i=1:length(list)
    hf1 = figure('color','w','position',[20,20,400,400],'Visible','off');
    genePlot = list{i}
    markergene = (data_orig_all_sorted(strcmpi(geneid_all,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),99);
    if tmpthlow==tmpthhigh
        tmpthlow = 0;
    end
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    axes('Position',[0.02,0.02,0.97,0.97])
    scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),30,markergene_color(~inpos,:),'.'); hold on;
    scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),30,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
%     title(genePlot);
    axis tight
    axis equal
    axis off
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png(['tsne_singlegene_211201/',genePlot,'_Glut12_tsneSingle_',date],gcf,1000)
    close(hf1)
end
%%

list = loadCellFile('markers_for_tsne_glut12.txt');
% [ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for i=1:length(list)
    hf1 = figure('color','w','position',[20,20,400,400],'Visible','off');
    genePlot = list{i}
    markergene = (data_orig_all_sorted(strcmpi(geneid_all,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),99);
    if tmpthlow==tmpthhigh
        tmpthlow = 0;
    end
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
    %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    axes('Position',[0.02,0.02,0.97,0.97])
    scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),30,markergene_color(~inpos,:),'.'); hold on;
    scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),30,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
%     title(genePlot);
    axis tight
    axis equal
    axis off
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    eval(['export_fig ',genePlot,'_Glut12_tsneSingle_',date,'.pdf']);
    close(hf1)
end