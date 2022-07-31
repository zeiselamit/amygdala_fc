tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;

load afterloading_noneuronal_QC2000_FC_16-Nov-2021

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);
% 1=immune',2=astro',3=vascular',4=oligo',5=blood',6=GABA',7=Glut1',8=Glut2',9=doublets'
cellid_gaba_glut = loadCellFile('cellid_nonneurons_step1_FC_QC2000_16-Nov-2021.txt');%'data','geneid','cellid','sample','fc_time','batch_flag','amy_flag','piri_flag'
[~,loc]= ismember(cellid_gaba_glut(:,1),cellid);
gaba_glut = zeros(size(cellid));
gaba_glut(loc) = cell2mat(cellid_gaba_glut(:,2));
% slc17a7 = data(strcmpi(geneid,'Slc17a7'),:);
% slc17a6 = data(strcmpi(geneid,'Slc17a6'),:);
% marker_exclud_cell = {'Slc17a7','Slc17a6','C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
%     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
% [~,loc] = ismember(marker_exclud_cell,geneid);
validcells = (gaba_glut==3 ) & amy_flag ;
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
sample = sample(validcells);
fc_time = fc_time(validcells);
amy_flag = amy_flag(validcells);
piri_flag = piri_flag(validcells);
batch_flag = batch_flag(validcells);

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
% not_relevant = {'Slc17a7','Slc17a6','C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
%     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};%,'Snca','Gm9843','Gm12027','Gm12031'
in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.8 & ~ismember(geneid(:,1),sex_genes)  ....
    & ~ismember(geneid(:,1),[IEG_markers]) );%& ~ismember(geneid(:,1),[not_relevant])

corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);
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
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',geneid(in(corr_filt(leaforder))))
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
table1 = [geneid(in(corr_filt(leaforder))),m2c(idx(leaforder))];
saveCellFile(table1,['gene_correlation_vascular_',date,'.txt'])
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
D = squareform(pdist(prj,'correlation'),'tomatrix');
[Dsort,XI] = sort(D,'ascend');
per_range = 500;
x = 1:per_range;%length(D);
optk = zeros(length(D),1);
for i=1:length(D)
    y = Dsort(1:per_range,i);
    x = x(:);
    y = y(:);
    a = atan((y(end)-y(1))/(x(end)-x(1)));
    xn = x*cos(a) + y*sin(a);
    yn = -x*sin(a) + y*cos(a);
    [~,imax] = max(yn);
    optk(i) = round((x(imax)));
end

perplexity = median(optk);

options = statset('MaxIter',1000);
mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
toc

% this is just the initial tsne, can be commented later
figure;
set(gcf,'color','w','position',[20,20,900,800])
plot(mapped_xy(:,1),mapped_xy(:,2),'.'); axis tight; axis off
% plot by sample
% sample_uni = {'8-1','8-2','10-1','10-2','18-1','18-2','19-1','23-1','23-2','23-3','41-1','41-2','42-1','45-1','45-2','45-3','46-1','46-2','46-3'};
colors = distinguishable_colors(length(sample_uni)+1);
figure;
set(gcf,'color','w','position',[20,20,1800,800]);
[ha, pos] = tight_subplot(1, 3, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
axes(ha(1))
for i=1:length(sample_uni)
    s = scatter(mapped_xy(strcmpi(sample,sample_uni{i}),1),mapped_xy(strcmpi(sample,sample_uni{i}),2),10,colors(i,:),'filled'); hold on;
%     alpha(s,0.4);
end
axis tight
axis equal
axis off
legend(sample_uni);
% % 
fc_time_uni = [0,2,8,24,28];
colors = distinguishable_colors(length(fc_time_uni)+1);
axes(ha(2))
for i=1:length(fc_time_uni)
    s = scatter(mapped_xy(fc_time==fc_time_uni(i),1),mapped_xy(fc_time==fc_time_uni(i),2),10,colors(i,:),'filled'); hold on;
%     alpha(s,0.3);
end
axis tight
axis equal
axis off
legend('ctrl','2hr','8hr','24hr','recall')
% plot by tissue
colors = [1,0,0;0,0,1;1,0,0;0,0,1]; %distinguishable_colors(length(sample_uni)+1);
axes(ha(3))
s = scatter(mapped_xy(amy_flag,1),mapped_xy(amy_flag,2),10,colors(1,:),'filled'); hold on;
% alpha(s,0.3);
s = scatter(mapped_xy(piri_flag,1),mapped_xy(piri_flag,2),10,colors(2,:),'filled'); hold on;
% alpha(s,0.3);
axis tight
axis off
axis equal
legend('Amy','Piri')
% % % % % % % % % % % % % % % % % % % % % % % % 
% MinPts = 10;
% eps_prc = 90;
% [idx10, ~,~] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);
% 
% MinPts = 40;
% eps_prc = 50;
% [idx40, ~,~] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);
% m = max(idx40);
% idx = idx40;
% for i=unique(idx10)'
%     ii=find(idx10==i);
%     if sum(idx40(ii))==0 
%         idx(ii) = m+1;
%         m = m+1;
%     end
% end
MinPts = 50;
eps_prc = 80;
[idx, ~,~] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);

colors = distinguishable_colors(length(unique(idx))+1);
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=unique(idx)'
    if i==0
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'o','color',colors(i+1,:),'markersize',3); hold on;
    elseif i>0
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i+1,:),'markersize',4); hold on;
    end
end
% for i=idxuni'
%     if i>=0
%         in = idx==i;
%         ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),num2str(i));
%         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
%     end
% end
axis tight;
axis equal
axis off

title(['MinPts=',num2str(MinPts),', peplexity=',num2str(perplexity),', epsprc=',num2str(eps_prc),',#C=',num2str(max(idx)),',#out=',num2str(sum(idx==0))],'fontsize',8);
% 
% if savefig_flag==1
%     savefig(gcf,['tsne_FC_0_2_8_24_Neurons_dbscan_',date,'.fig'])
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% sort the data by the clusters and remove outliers
[idx,xi] = sort(idx);
xi(idx==0) = [];
idx(idx==0) = [];
idxuni = unique(idx);
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
xi = [1:length(idx)];
for i=1:length(idxuni)
    i
    ind = find(idx==i);
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

meangr_mat = zeros(length(moldata(:,1)),length(idxuni));
clust_cent = zeros(length(idxuni),2);
for jjj=1:length(idxuni)
    jjj
    meangr_mat(:,jjj) = mean(log2(data_sorted_all(:,idx==idxuni(jjj))+1),2);
    clust_cent(jjj,:) = [median(mapped_xy(idx==idxuni(jjj),1)),median(mapped_xy(idx==idxuni(jjj),2))];
end

meangr_mat1 = meangr_mat;

[prj,m,D,V,Q] = pca_wis(meangr_mat1',size(meangr_mat1,2)-1);
Zpca = linkage(prj,'ward','correlation');
Dpca = pdist(prj,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
figure;
set(gcf,'position',[100,100,1000,1000],'color','w')
h1 = axes('position',[0.03,0.03,0.3,0.93]);
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
h2 = axes('position',[0.35,0.03,0.63,0.93]);
x=squareform(Dpca); imagesc(x(leaforder_pca,leaforder_pca));
colormap('summer')
set(gca,'ytick',[1:length(leaforder_pca)],'xtick',[],'fontsize',8,'ydir','normal')
linkaxes([h1,h2],'y');
% if savefig_flag==1
%     savefig(gcf,['tree_FC_GABA_',date,'.fig'])    
% end
% if savefig_pdf==1
%     eval(['export_fig tree_FC_GABA_',date,'.pdf']);
% end
leaforder = leaforder_pca;

T_cells_tmp_new = zeros(length(idx),2);
for i=1:length(leaforder)
    T_cells_tmp_new(idx==idxuni(leaforder(i)),1) = i;
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

% % % % % % % % % % % % % % 
T_cells_tmp = T_cells_tmp_new;
T_cells_tmp_uni = unique(T_cells_tmp);
% % % % % % % % % % % % % % % 
% [table1, table2] = dendrogram_split_markers(cellfun(@(x) num2str(x), m2c(T_cells_tmp_uni),'UniformOutput',0).....
%     ,cellfun(@(x) num2str(x), m2c(T_cells_tmp(:,1)),'UniformOutput',0),Zpca_post,data_sorted_all,geneid);
% saveCellFile(table1,['FC_GABA_dendrogram_junction_split_markers_',date,'.txt']);
% saveCellFile(table2,['FC_GABA_dendrogram_junction_split_markers_by_average_',date,'.txt']);
% ind = [[3:53:length(table1(:,1))];[4:53:length(table1(:,1))]];
% ind = ind(:);
% treemark1 = table1(ind,[1,7]);
% treemark1 = [reshape(treemark1(:,1),2,[]);reshape(treemark1(:,2),2,[])];
% treemark1 = flipud(treemark1(:));
% rmv = [];
% for i=2:length(treemark1)
%     if sum(strcmpi(treemark1{i},treemark1(1:i-1)))>0
%         rmv = [rmv,i];
%     end
% end
% treemark1(rmv) = [];
% % % % % % % % % % % % % % % 
colors = distinguishable_colors(length(T_cells_tmp_uni)+1);
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(T_cells_tmp_uni)
    if i==-1
        ii=find(T_cells_tmp(:,1)==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',[0,0,1],'markersize',3); hold on;
    else
        ii=find(T_cells_tmp(:,1)==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i+1,:),'markersize',5); hold on;
    end
end
for i=1:length(T_cells_tmp_uni)
    in = T_cells_tmp(:,1)==i;
    ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),num2str(i));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',8)
end
axis tight;
axis equal
axis off
title(['MinPts=',num2str(MinPts),', epsprc=',num2str(eps_prc),',#C=',num2str(max(T_cells_tmp(:,1))),',#out=',num2str(sum(T_cells_tmp(:,1)==0))],'fontsize',8);
% if savefig_flag==1
%     savefig(gcf,['tsne_FC_GABA_by_cluster_',date,'.fig'])
%     % eval(['export_fig tsne_AmyPiri_FC_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);
% end
% if savefig_pdf==1
%     eval(['export_fig tsne_FC_GABA_by_cluster_',date,'.pdf']);
% end

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
% if savefig_flag==1
%     % eval(['export_fig tsne_AmyPiri_FC_by_sampleannot_','_',date,'.pdf']);
%     savefig(gcf,['tsne_FC_GABA_by_annot_perplexity_',num2str(perplexity),'_PCAdim',num2str(initial_dims),'_',date,'.fig'])
% end

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
ax1 = axes('position',[0.1,0.02,0.88,0.73]);
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

% sample_uni = {'8-1','8-2','10-1','10-2','18-1','18-2','19-1','23-1','23-2','23-3','41-1','41-2','42-1','42-2','45-1','45-2','45-3','46-1','46-2','46-3'};
samples_num = false(length(sample_uni),length(sample_sorted));
for i=1:length(sample_uni)
    samples_num(i, strcmpi(sample_sorted,sample_uni{i})) = true;
end

ax2_1 = axes('position',[0.1,0.75,0.88,0.012]);
imagesc(~(fc_time_sorted'==0)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','0');

ax2_2 = axes('position',[0.1,0.75+1*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==2)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','2');

ax2_3 = axes('position',[0.1,0.75+2*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==8)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','8');

ax2_4 = axes('position',[0.1,0.75+3*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==24)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','24');

ax2_5 = axes('position',[0.1,0.75+4*0.012,0.88,0.012]);
imagesc(~(fc_time_sorted'==28)); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','recall');

ax3 = axes('position',[0.1,0.81,0.88,0.01]);
imagesc(~amy_flag_sorted'); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Amy');

ax4 = axes('position',[0.1,0.82,0.88,0.01]);
imagesc(~piri_flag_sorted'); hold on;
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Piri');

ax5 = axes('position',[0.1,0.83,0.88,0.12]);
imagesc(~samples_num); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1:length(sample_uni)],'yticklabel',sample_uni,'fontsize',6);

Aqp4 = data_orig_all_sorted(strcmpi(geneid_all,'Aqp4'),:);
Gja1 = data_orig_all_sorted(strcmpi(geneid_all,'Gja1'),:);
Foxj1 = data_orig_all_sorted(strcmpi(geneid_all,'Foxj1'),:);

ax6 = axes('position',[0.1,0.95,0.88,0.01]);
imagesc(~Aqp4); hold on;
axes(ax6)
colormap('gray');
freezeColors(ax6);
set(ax6,'xtick',[],'ytick',[1],'yticklabel','Aqp4','fontsize',5);
ax7 = axes('position',[0.1,0.95+0.01,0.88,0.01]);
imagesc(~Gja1); hold on;
colormap('gray');
freezeColors(ax7);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Gja1','fontsize',6);
ax8 = axes('position',[0.1,0.95+2*0.01,0.88,0.01]);
imagesc(~Foxj1); hold on;
colormap('gray');
freezeColors(ax8);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Foxj1','fontsize',6);

linkaxes([ax1,ax2_1,ax2_2,ax2_3,ax2_4,ax2_5,ax3,ax4,ax5,ax6,ax7,ax8],'x');
% if savefig_flag==1
%     savefig(gcf,['markertable_FC_GABA_',date,'.fig'])
%     % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
% end
% if savefig_pdf==1
%     eval(['export_fig markertable_FC_GABA_',date,'.pdf']);
% end

table1 = [cellid_sorted,m2c(T_cells_tmp)];
saveCellFile(table1,['cellidCluster_step2_vascular_FC_',date,'.txt']);

