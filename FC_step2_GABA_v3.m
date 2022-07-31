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

cellid_gaba_glut = loadCellFile('cellid_gaba_glut_FC_QC3000_06-Oct-2021.txt');%'data','geneid','cellid','sample','fc_time','batch_flag','amy_flag','piri_flag'
[~,loc]= ismember(cellid_gaba_glut(:,1),cellid);
gaba_glut = zeros(size(cellid));
gaba_glut(loc) = cell2mat(cellid_gaba_glut(:,2));
% slc17a7 = data(strcmpi(geneid,'Slc17a7'),:);
% slc17a6 = data(strcmpi(geneid,'Slc17a6'),:);
marker_exclud_cell = {'Slc17a7','Slc17a6','C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
    ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
[~,loc] = ismember(marker_exclud_cell,geneid);
validcells = (gaba_glut==1 ) & amy_flag & (sum(data(loc,:)>1)')<2;
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
not_relevant = {'Slc17a7','Slc17a6','C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
    ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};%,'Snca','Gm9843','Gm12027','Gm12031'
in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),sex_genes)  ....
    & ~ismember(geneid(:,1),[IEG_markers]) & ~ismember(geneid(:,1),[not_relevant]));

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
saveCellFile(table1,['gene_correlation_GABA_',date,'.txt'])
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
MinPts = 40;
eps_prc = 50;
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

[prj,m,D,V,Q] = pca_wis(meangr_mat1',initial_dims);
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
if savefig_flag==1
    savefig(gcf,['tree_FC_GABA_',date,'.fig'])    
end
if savefig_pdf==1
    eval(['export_fig tree_FC_GABA_',date,'.pdf']);
end
leaforder = leaforder_pca;
Zpca_post = linkage(prj(leaforder_pca,:),'ward','correlation');

bigclust = cluster(Zpca_post,'maxclust',10);
T_cells_tmp_new = zeros(length(idx),2);
for i=1:length(leaforder)
    T_cells_tmp_new(idx==idxuni(leaforder(i)),1) = i;
    T_cells_tmp_new(idx==idxuni(leaforder(i)),2) = bigclust(i);
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
[table1, table2] = dendrogram_split_markers(cellfun(@(x) num2str(x), m2c(T_cells_tmp_uni),'UniformOutput',0).....
    ,cellfun(@(x) num2str(x), m2c(T_cells_tmp(:,1)),'UniformOutput',0),Zpca_post,data_sorted_all,geneid);
saveCellFile(table1,['FC_GABA_dendrogram_junction_split_markers_',date,'.txt']);
saveCellFile(table2,['FC_GABA_dendrogram_junction_split_markers_by_average_',date,'.txt']);
ind = [[3:53:length(table1(:,1))];[4:53:length(table1(:,1))]];
ind = ind(:);
treemark1 = table1(ind,[1,7]);
treemark1 = [reshape(treemark1(:,1),2,[]);reshape(treemark1(:,2),2,[])];
treemark1 = flipud(treemark1(:));
rmv = [];
for i=2:length(treemark1)
    if sum(strcmpi(treemark1{i},treemark1(1:i-1)))>0
        rmv = [rmv,i];
    end
end
treemark1(rmv) = [];
% % % % % % % % % % % % % % 
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
if savefig_flag==1
    savefig(gcf,['tsne_FC_GABA_by_cluster_',date,'.fig'])
    % eval(['export_fig tsne_AmyPiri_FC_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);
end
if savefig_pdf==1
    eval(['export_fig tsne_FC_GABA_by_cluster_',date,'.pdf']);
end

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
    % eval(['export_fig tsne_AmyPiri_FC_by_sampleannot_','_',date,'.pdf']);
    savefig(gcf,['tsne_FC_GABA_by_annot_perplexity_',num2str(perplexity),'_PCAdim',num2str(initial_dims),'_',date,'.fig'])
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

gad2 = data_orig_all_sorted(strcmpi(geneid_all,'Gad2'),:);
slc17a7 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a7'),:);
slc17a6 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a6'),:);

ax6 = axes('position',[0.1,0.95,0.88,0.01]);
imagesc(~gad2); hold on;
axes(ax6)
colormap('gray');
freezeColors(ax6);
set(ax6,'xtick',[],'ytick',[1],'yticklabel','Gad2','fontsize',5);
ax7 = axes('position',[0.1,0.95+0.01,0.88,0.01]);
imagesc(~slc17a6); hold on;
colormap('gray');
freezeColors(ax7);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a6','fontsize',6);
ax8 = axes('position',[0.1,0.95+2*0.01,0.88,0.01]);
imagesc(~slc17a7); hold on;
colormap('gray');
freezeColors(ax8);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a7','fontsize',6);

linkaxes([ax1,ax2_1,ax2_2,ax2_3,ax2_4,ax2_5,ax3,ax4,ax5,ax6,ax7,ax8],'x');
if savefig_flag==1
    savefig(gcf,['markertable_FC_GABA_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
if savefig_pdf==1
    eval(['export_fig markertable_FC_GABA_',date,'.pdf']);
end

table1 = [cellid_sorted,m2c(T_cells_tmp)];
saveCellFile(table1,['cellidCluster_step2_GABA_FC_',date,'.txt']);

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
ax1 = axes('position',[0.1,0.02,0.88,0.73]);
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

gad2 = data_orig_all_sorted(strcmpi(geneid_all,'Gad2'),:);
slc17a7 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a7'),:);
slc17a6 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a6'),:);

ax6 = axes('position',[0.1,0.95,0.88,0.01]);
imagesc(~gad2); hold on;
axes(ax6)
colormap('gray');
freezeColors(ax6);
set(ax6,'xtick',[],'ytick',[1],'yticklabel','Gad2','fontsize',5);
ax7 = axes('position',[0.1,0.95+0.01,0.88,0.01]);
imagesc(~slc17a6); hold on;
colormap('gray');
freezeColors(ax7);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a6','fontsize',6);
ax8 = axes('position',[0.1,0.95+2*0.01,0.88,0.01]);
imagesc(~slc17a7); hold on;
colormap('gray');
freezeColors(ax8);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a7','fontsize',6);

linkaxes([ax1,ax2_1,ax2_2,ax2_3,ax2_4,ax2_5,ax3,ax4,ax5,ax6,ax7,ax8],'x');
if savefig_flag==1
    savefig(gcf,['junctionmarkers_FC_GABA_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
% % % % % % % % % 
meangr_mat_0 = zeros(length(moldata(:,1)),length(idxuni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    meangr_mat_0(:,jjj) = mean(data_sorted_all(:,T_cells_tmp(:,1)==T_cells_tmp_uni(jjj))>0,2);
end
% meangr_mat1 = meangr_mat1(:,leaforder);
datamarkers = meangr_mat_0(ind_gr_tmp_mark,:);
% d = squareform(pdist(datamarkers,'correlation'));
% z = linkage(datamarkers,'ward','correlation');
% lo = optimalleaforder(z,d);
datamarkers_cn = datamarkers;%cent_norm(datamarkers);(lo,:)
% for i=1:length(ind_gr_tmp_mark)
%     t = sort(datamarkers_cn(i,:),'descend');
%     c(i) = t(findknee(t));
% end
% datamarkers_cn =  datamarkers_cn >= repmat(c',1,length(datamarkers_cn(1,:)));

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.02,0.88,0.75]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
set(gca,'xtick',[1:length(gr_center)],'xticklabel',[1:length(gr_center)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)
colormap('summer');
freezeColors(gca);

h2 = axes('position',[0.1,0.79,0.88,0.21]);
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','top');
axis off
set(gca,'xlim',[0.5,length(leaforder_pca)+0.5])
linkaxes([ax1,h2],'x')

if savefig_flag==1
    savefig(gcf,['junctionmarkers_clusteraverage_FC_GABA_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
meangr_mat_0 = zeros(length(geneid_all),length(idxuni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    meangr_mat_0(:,jjj) = mean(data_orig_all_sorted(:,T_cells_tmp(:,1)==T_cells_tmp_uni(jjj))>0,2);
end
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
    savefig(gcf,['CEA_markers_clusteraverage_FC_GABA_',date,'.fig'])
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
    savefig(gcf,['N_Z_list_clusteraverage_FC_GABA_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
    ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
[~,loc] = ismember(exclude_markers,geneid_all);
nonneuro = sum(data_orig_all_sorted(loc,:));

gabaglut = zeros(size(T_cells_tmp_uni));
gabaglut_sc = zeros(size(T_cells_tmp(:,1)));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    [~,gabaglut(jjj)] = max([mean(gad2(:,T_cells_tmp(:,1)==jjj)>0);mean(slc17a7(:,T_cells_tmp(:,1)==jjj)>0);....
        mean(slc17a6(:,T_cells_tmp(:,1)==jjj)>0);mean(nonneuro(:,T_cells_tmp(:,1)==jjj)>0)]);
    gabaglut_sc(T_cells_tmp(:,1)==jjj) = gabaglut(jjj);
end
figure;
set(gcf,'color','w','position',[20,20,800,800])
colors = distinguishable_colors(4);
for idx=unique(gabaglut_sc)'
    ii=find(gabaglut_sc==idx); 
    h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(idx,:),'markersize',3); hold on;
end
axis tight
axis equal
axis off
title('GABA/Glut1/Glut2')
legend('GABA','Glut1','Glut2','non-neurons')
% % % % % % % % % % % % % % % % % % % % % % % % % %
% markers = {'fos','arc','npas4','rrad'};
% t1 = [];
% t2 = [];
% t3 = [];
% for ii=1:length(markers)
%     tmpmark = markers{ii};
%     marker = data_orig_all_sorted(strcmpi(geneid_all,tmpmark),:);
%     marker_percent = zeros(1, length(T_cells_tmp_uni));
%     marker_ave = zeros(1, length(T_cells_tmp_uni));
%     p_marker = zeros(1, length(T_cells_tmp_uni));
%     M = length(marker);
%     K = sum(marker>0);
%     for j=1:length(T_cells_tmp_uni)
%         c1 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j));
%         c2 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j) & marker'>0);
%         marker_percent(j) = 100*c2/c1;
%         marker_ave(j) = mean(marker(T_cells_tmp(:,1)==T_cells_tmp_uni(j)));
%         p_marker(j) = hygecdf(c2,M,K,c1,'upper');
%     end
%     marker_percent(isnan(marker_percent)) = 0;
%     marker_ave(isnan(marker_ave)) = 0;
%     eval([markers{ii},'_percent = marker_percent;'])
%     eval([markers{ii},'_ave = marker_ave;'])
%     eval([markers{ii},'_p = p_marker;']);
%     eval(['t1 = [t1;',markers{ii},'_percent];'])
%     eval(['t2 = [t2;',markers{ii},'_ave];'])
%     eval(['t3 = [t3;',markers{ii},'_p];'])
% end
% 
% figure;
% set(gcf,'color','w','position',[20,20,900,1800])
% [ha, pos] = tight_subplot(1, 5, [0.02,0.00], [0.02,0.02], [0.05,0.02]);
% axes(ha(1));
% barh(t1');
% set(gca,'ytick',[1:length(idxuni)])
% title('%>0')
% 
% axes(ha(2));
% barh(t2');
% set(gca,'ytick',[1:length(idxuni)],'YTickLabel',cell(length(idxuni),1),'XLim',[0,16])
% title('Average')
% 
% axes(ha(3));
% barh(-log10(t3)');
% set(gca,'ytick',[1:length(idxuni)],'YTickLabel',cell(length(idxuni),1),'XLim',[0,10])
% legend(markers);
% title('-log10(p)')
% % % % % % % % % % % % % % %  % % % % % % % % % %
% fc_time_uni = [0,2,8,24,28];
% t1 = [];
% t2 = [];
% for ii=1:length(fc_time_uni)
%     marker = fc_time_sorted==fc_time_uni(ii);
%     marker_percent = zeros(1, length(T_cells_tmp_uni));
%     p_marker = zeros(1, length(T_cells_tmp_uni));
%     M = length(marker);
%     K = sum(marker>0);
%     for j=1:length(T_cells_tmp_uni)
%         c1 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j));
%         c2 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j) & marker>0);
%         marker_percent(j) = 100*c2/c1;
%         p_marker(j) = hygecdf(c2,M,K,c1,'upper');
%     end
%     marker_percent(isnan(marker_percent)) = 0;
% %     eval([markers{ii},'_percent = marker_percent;'])
% %     eval([markers{ii},'_p = p_marker;']);
%     t1 = [t1;marker_percent];
%     t2 = [t2;p_marker];
% end
% 
% axes(ha(4));
% barh(t1');
% set(gca,'ytick',[1:length(idxuni)],'YTickLabel',cell(length(idxuni),1),'XLim',[0,100])
% title('%>0')
% 
% axes(ha(5));
% barh(-log10(t2)');
% legend({'0','2','8','24','recall'});
% set(gca,'ytick',[1:length(idxuni)],'YAxisLocation','right','XLim',[0,10])
% title('-log10(p)')
% 
% linkaxes([ha(1),ha(2),ha(3),ha(4),ha(5)],'y');
% if savefig_flag==1
%     % eval(['export_fig tsne_markers2_AmyPiri_FC_',date,'.pdf']);
%     savefig(gcf,['timepoints_FC_Glut1_enrichment_percluster_',date,'.fig'])
% end
% 
% 
% 
% markers = IEG_markers;
% t1 = [];
% t2 = [];
% for ii=1:length(markers)
%     tmpmark = markers{ii};
%     marker = data_orig_all_sorted(strcmpi(geneid_all,tmpmark),:);
%     marker_percent = zeros(1, length(T_cells_tmp_uni));
%     p_marker = zeros(1, length(T_cells_tmp_uni));
%     M = length(marker);
%     K = sum(marker>0);
%     for j=1:length(T_cells_tmp_uni)
%         c1 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j));
%         c2 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j) & marker'>0);
%         marker_percent(j) = 100*c2/c1;
%         p_marker(j) = hygecdf(c2,M,K,c1,'upper');
%     end
%     marker_percent(isnan(marker_percent)) = 0;
%     eval([markers{ii},'_percent = marker_percent;'])
%     eval([markers{ii},'_p = p_marker;']);
%     eval(['t1 = [t1;',markers{ii},'_percent];'])
%     eval(['t2 = [t2;',markers{ii},'_p];'])
% end
% figure;
% set(gcf,'color','w','position',[20,20,900,1200])
% [ha, pos] = tight_subplot(1, 2, [0.05,0.00], [0.05,0.05], [0.05,0.02]);
% axes(ha(1));
% imagesc(cent_norm(t1)');
% colormap('summer');
% set(gca,'ytick',[1:length(idxuni)],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
% axes(ha(2));
% imagesc(-log10(t2)',[0,10]);
% colormap('summer');
% set(gca,'ytick',[],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
% linkaxes([ha(1),ha(2)],'y');
% 
% if savefig_flag==1
%     % eval(['export_fig tsne_markers2_AmyPiri_FC_',date,'.pdf']);
%     savefig(gcf,['immidiateEarly_FC_Glut1_enrichment_percluster_',date,'.fig'])
% end
% 
% 
% % % % % % % % % % % % % %
% list = {'Slc17a6','Slc17a7','Sim1','Arhgdib','Pou3f1','Ppp1r1b','Erbb4','Col11a1','Cxcl14','Cbln2','C1ql1','Rorb','Nxph1','Rims3','Tmem132c','Arhgap36'};
% figure;
% set(gcf,'color','w','position',[20,20,1100,960])
% [ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% for i=1:length(list)
%     genePlot = list{i};
%     markergene = (data_orig_all_sorted(strcmpi(geneid_all,genePlot),:));
%     inpos = markergene>0;
%     tmpthlow = prctile(markergene(markergene>0),1);
%     tmpthhigh = prctile(markergene(markergene>0),99);
%     if tmpthlow==tmpthhigh
%         tmpthlow = 0;
%     end
%     markergene(markergene>tmpthhigh) = tmpthhigh;
%     markergene(markergene<tmpthlow) = tmpthlow;
%     c_rgb = [1,0,0];rand([1,3]);
%     %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%     %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
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
%     % eval(['export_fig tsne_markers2_AmyPiri_FC_',date,'.pdf']);
%     savefig(gcf,['tsne_FC_Glut1_markers1_',date,'.fig'])
% end
% if savefig_pdf==1
%     eval(['export_fig tsne_FC_Glut1_markers1_',date,'.pdf']);
% end
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %
%     list = {'Trhr','Sertm1','Neurod6','Cyp26b1','Htr2c','Trh','Itgav','Sox5','Reln','Met','Grik3','Rorb','Slc1a3','Sema3e','Sim1','Nos1'};
%     figure;
%     set(gcf,'color','w','position',[20,20,1100,960])
%     [ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
%     for i=1:length(list)
%         genePlot = list{i};
%         markergene = (data_orig_all_sorted(strcmpi(geneid_all,genePlot),:));
%         inpos = markergene>0;
%         tmpthlow = prctile(markergene(markergene>0),1);
%         tmpthhigh = prctile(markergene(markergene>0),99);
%         if tmpthlow==tmpthhigh
%             tmpthlow = 0;
%         end
%         markergene(markergene>tmpthhigh) = tmpthhigh;
%         markergene(markergene<tmpthlow) = tmpthlow;
%         c_rgb = [1,0,0];rand([1,3]);
%         %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%         %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
%         markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
%             interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
%             ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
%         axes(ha(i));
%         scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
%         scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
%         set(gca,'xlim',[-150,150],'ylim',[-150,150])
%         title(genePlot);
%         axis tight
%         axis equal
%         axis off
%     end
% if savefig_flag==1
%     % eval(['export_fig tsne_markers2_AmyPiri_FC_',date,'.pdf']);
%     savefig(gcf,['tsne_FC_Glut1_markers1_',date,'.fig'])
% end
% if savefig_pdf==1
%     eval(['export_fig tsne_FC_Glut1_markers1_',date,'.pdf']);
% end
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %
% list = {'Rspo2','Rspo3','Tshz1','Tshz2','Trh','Ecel1','Tgfb2','Foxp2','Neurod6','Lypd1','Meis2','Bdnf','Gad2','Arc','Jun','Fos'};
% figure;
% set(gcf,'color','w','position',[20,20,1100,960])
% [ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% for i=1:length(list)
%     genePlot = list{i};
%     markergene = (data_orig_all_sorted(strcmpi(geneid_all,genePlot),:));
%     inpos = markergene>0;
%     tmpthlow = prctile(markergene(markergene>0),1);
%     tmpthhigh = prctile(markergene(markergene>0),99);
%     if tmpthlow==tmpthhigh
%         tmpthlow = 0;
%     end
%     markergene(markergene>tmpthhigh) = tmpthhigh;
%     markergene(markergene<tmpthlow) = tmpthlow;
%     c_rgb = [1,0,0];rand([1,3]);
%     %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%     %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
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
%     % eval(['export_fig tsne_markers2_AmyPiri_FC_',date,'.pdf']);
%     savefig(gcf,['tsne_FC_Glut1_markers1_',date,'.fig'])
% end
% if savefig_pdf==1
%     eval(['export_fig tsne_FC_Glut1_markers1_',date,'.pdf']);
% end
% % % % % % % % % % % % % %
% 
% fc_time_uni = [0,2,8,24,28];
% data_mean_tp = zeros(length(geneid_all),length(fc_time_uni));
% data_frac_tp = zeros(length(geneid_all),length(fc_time_uni));
% for ii=1:length(fc_time_uni)
%     data_mean_tp(:,ii) = mean(log2(data_orig_all_sorted(:, fc_time_sorted==fc_time_uni(ii))+1),2);
%     data_frac_tp(:,ii) = mean(data_orig_all_sorted(:, fc_time_sorted==fc_time_uni(ii))>0,2);
% end
% 
% % 0 vs 2
% top_g = 50;
% x1 = data_mean_tp(:,1);
% x2 = data_mean_tp(:,2);
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
% title(['cluster ',num2str(c2),' vs ',num2str(c1)])
% xlabel(['mean (0hr)'])
% ylabel(['mean (2hr)'])
% axis tight
% x1 = data_frac_tp(:,1);
% x2 = data_frac_tp(:,2);
% d = x1-x2 ;
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
% title(['cluster ',num2str(c2),' vs ',num2str(c1)])
% xlabel(['frac (0hr)'])
% ylabel(['frac (2hr)'])
% axis tight
% % % % % % % % % % % % % % % % % % % % % 
% 
% % in = find(sum(data_orig_all_sorted>0,2)>100 & sum(data_orig_all_sorted>0,2)<length(data_orig_all_sorted(1,:))*0.5 & ~ismember(geneid_all(:,1),sex_genes)  ....
% %     & ~ismember(geneid_all(:,1),[batch_genes]));
% in = find(sum(data_orig_all_sorted>0,2)>100 & sum(data_orig_all_sorted>0,2)<length(data_orig_all_sorted(1,:))*0.9);
% ave_flag = 1;
% b = 1;
% c_top = cell(length(unique(T_cells_tmp(:,b))),1);
% c_top2 = cell(length(unique(T_cells_tmp(:,b))),1);
% c_top3 = cell(length(unique(T_cells_tmp(:,b))),1);
% c_top_rnd = cell(length(unique(T_cells_tmp(:,b))),1);
% fc_time_rnd = fc_time_sorted;
% for c1=1:length(unique(T_cells_tmp(:,b)))
%     c1
%     fc_time_uni = [0,2,8,24,28];
%     top_g = 20;
%     t = fc_time_rnd(T_cells_tmp(:,b)==c1);
%     fc_time_rnd(T_cells_tmp(:,b)==c1) = t(randperm(length(t)));
%     gr1 = find(T_cells_tmp(:,b)==c1 & fc_time_sorted==fc_time_uni(1));
%     gr1_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(1));
%     x1 = mean(log2(data_orig_all_sorted(:,gr1)+1),2);
%     x1_rnd = mean(log2(data_orig_all_sorted(:,gr1_rnd)+1),2);
%     p90gene = prctile(data_orig_all_sorted(:,T_cells_tmp(:,b)==c1),90,2);
%     k = 0;
%     qmat = ones(length(geneid_all),4);
%     qmat2 = ones(length(geneid_all),4);
%     qmat3 = ones(length(geneid_all),4);
%     qmat_rnd = ones(length(geneid_all),4);
%     for i=[2,3,4,5]
%         k = k+1;
%         gr2 = find(T_cells_tmp(:,b)==c1 & fc_time_sorted==fc_time_uni(i));
%         gr2_not = find(T_cells_tmp(:,b)==c1 & fc_time_sorted~=fc_time_uni(i));
%         gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(i));
%         x2 = mean(log2(data_orig_all_sorted(:,gr2)+1),2);
%         x2_rnd = mean(log2(data_orig_all_sorted(:,gr2_rnd)+1),2);
%         d = x1-x2 ;
%         x90 = sum(data_orig_all_sorted(:,gr2)>repmat(p90gene,1,length(gr2)),2);
%         x90(p90gene==0) = nan;
%         pbino = binocdf(x90,length(gr2),0.1,'upper');
%         pbino(isnan(pbino)) = 1;
%         %         [~,ptt] = ttest2(data_orig_all_sorted(in,gr1)',data_orig_all_sorted(in,gr2)','tail','both');
%         ptt = ones(length(in),1);
%         ptt2 = ones(length(in),1);
%         ptt_rnd = ones(length(in),1);
%         if length(gr1)>10 & length(gr2)>10
%             for s=1:length(in)
%                 %             s
%                 %             ptt(s) = ranksum(data_orig_all_sorted(in(s),gr1),data_orig_all_sorted(in(s),gr2),'tail','left');
%                 [~,ptt(s)] = kstest2(data_orig_all_sorted(in(s),gr1),data_orig_all_sorted(in(s),gr2),'tail','larger');
%                 [~,ptt2(s)] = kstest2(data_orig_all_sorted(in(s),gr2_not),data_orig_all_sorted(in(s),gr2),'tail','larger');
%                 [~,ptt_rnd(s)] = kstest2(data_orig_all_sorted(in(s),gr1_rnd),data_orig_all_sorted(in(s),gr2_rnd),'tail','larger');
%             end
%         end
%         ptt(isnan(ptt)) = 1;
%         ptt2(isnan(ptt2)) = 1;
%         ptt_rnd(isnan(ptt_rnd)) = 1;
%         %         in2 = (in(ptt<1));
%         %         q = qval_from_pval(ptt(ptt<1));
%         %         q(q<1e-10) = 1e-10;
%         qmat(in,k) = ptt;
%         qmat2(in,k) = ptt2;
%         qmat3(in,k) = pbino(in);
%         qmat_rnd(in,k) = ptt_rnd;
%     end
%     c_av_tp = zeros(length(geneid_all),length(fc_time_uni));
%     c_av_tp_rnd = zeros(length(geneid_all),length(fc_time_uni));
%     rank_tp = zeros(length(geneid_all),length(fc_time_uni));
%     [~,xii] = sort(data_orig_all_sorted(:,T_cells_tmp(:,b)==c1)');
%     tp_c1 = fc_time_sorted(T_cells_tmp(:,b)==c1);
%     xii = xii';
%     grsize = zeros(5,1);
%     for i=[1,2,3,4,5]
%         gr2 = find(T_cells_tmp(:,b)==c1 & fc_time_sorted==fc_time_uni(i));
%         gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(i));
%         grsize(i) = length(gr2);
%         if ave_flag==1
%             c_av_tp(:,i) = median(data_orig_all_sorted(:,gr2),2);
%             c_av_tp_rnd(:,i) = median(data_orig_all_sorted(:,gr2_rnd),2);
%         else
%             c_av_tp(:,i) = mean(data_orig_all_sorted(:,gr2)>0,2);
%         end
%         rank_tp(:,i) = mean(xii(:,tp_c1==fc_time_uni(i)),2)/length(tp_c1);
%     end
%     
%     [qsort,xi] = sort(qmat);
%     [qsort2,xi2] = sort(qmat2);
%     [qsort3,xi3] = sort(qmat3);
%     [qsort_rnd,xi_rnd] = sort(qmat_rnd);
%     meddiff = c_av_tp(:,2:5)-repmat(c_av_tp(:,1),1,4);
%     meddiff_rnd = c_av_tp_rnd(:,2:5)-repmat(c_av_tp_rnd(:,1),1,4);
%     c_top{c1} = [xi(1:200,:),qsort(1:200,:),[meddiff(xi(1:200,1),1),meddiff(xi(1:200,2),2),meddiff(xi(1:200,3),3),meddiff(xi(1:200,4),4)]];
%     c_top2{c1} = [xi2(1:200,:),qsort2(1:200,:),[meddiff(xi(1:200,1),1),meddiff(xi(1:200,2),2),meddiff(xi(1:200,3),3),meddiff(xi(1:200,4),4)]];
%     c_top3{c1} = [xi3(1:200,:),qsort3(1:200,:),[meddiff(xi(1:200,1),1),meddiff(xi(1:200,2),2),meddiff(xi(1:200,3),3),meddiff(xi(1:200,4),4)]];
%     c_top_rnd{c1} = [xi_rnd(1:200,:),qsort_rnd(1:200,:),[meddiff_rnd(xi(1:200,1),1),meddiff_rnd(xi(1:200,2),2),meddiff_rnd(xi(1:200,3),3),meddiff_rnd(xi(1:200,4),4)]];
%     %     xi = xi(1:10,:);
%     %     xi = xi(:);
%     %     hh = figure('position',[200,200,1600,1000],'color','w','Name',['cluster',num2str(c1)],'Visible','on');
%     %     [ha, pos] = tight_subplot(4, 10, [0.03,0.03], [0.1,0.05], [0.05,0.05]);
%     %     for i=1:length(ha)
%     % %         axes(ha(i));
%     %         set(hh,'CurrentAxes',ha(i))
%     %         t_ed = zeros(5,1);
%     %         t_av = zeros(5,1);
%     %         t_75 = zeros(5,1);
%     %         for j=[1:5]
%     %             gr2 = find(T_cells_tmp(:,b)==c1 & fc_time_sorted==fc_time_uni(j));
%     %             plot(j+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(xi(i),gr2)','.'); hold on;
%     %             t_ed(j) = median(data_orig_all_sorted(xi(i),gr2));
%     %             t_av(j) = mean(data_orig_all_sorted(xi(i),gr2));
%     %             t_75(j) = prctile(data_orig_all_sorted(xi(i),gr2),75);
%     %         end
%     %         plot(t_ed,'-sk'); hold on;
%     %         plot(t_av,'-or'); hold on;
%     %         plot(t_75,'-dg'); hold on;
%     %         set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear')
%     %         axis tight
%     %         title(geneid_all(xi(i)));
%     %     end
%     %     axes('Position',[0.1,0.02,0.2,0.1]);
%     %     text(0,0,['0(',num2str(grsize(1)),'),2(',num2str(grsize(2)),'),8(',.....
%     %         num2str(grsize(3)),'),24(',num2str(grsize(4)),'),recall(',num2str(grsize(4)),')'])
%     %     axis off
%     %     eval(['export_fig GABA_fine_cluster',num2str(c1),'_topDEG_',date,'.pdf']);
%     %     pause
%     %     xi_rnd = xi_rnd(1:10,:);
%     %     xi_rnd = xi_rnd(:);
%     %     hh_rnd = figure('position',[200,200,1600,1000],'color','w','Name',['cluster',num2str(c1)],'Visible','on');
%     %     [ha, pos] = tight_subplot(4, 10, [0.03,0.03], [0.1,0.05], [0.05,0.05]);
%     %     for i=1:length(ha)
%     % %         axes(ha(i));
%     %         set(hh_rnd,'CurrentAxes',נha(i))
%     %         t_ed = zeros(5,1);
%     %         t_av = zeros(5,1);
%     %         t_75 = zeros(5,1);
%     %         for j=[1:5]
%     %             gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(j));
%     %             plot(j+0.5*rand(length(gr2_rnd),1)-0.1,0.5*rand(length(gr2_rnd),1)-0.1+data_orig_all_sorted(xi(i),gr2_rnd)','.'); hold on;
%     %             t_ed(j) = median(data_orig_all_sorted(xi(i),gr2_rnd));
%     %             t_av(j) = mean(data_orig_all_sorted(xi(i),gr2_rnd));
%     %             t_75(j) = prctile(data_orig_all_sorted(xi(i),gr2_rnd),75);
%     %         end
%     %         plot(t_ed,'-sk'); hold on;
%     %         plot(t_av,'-or'); hold on;
%     %         plot(t_75,'-dg'); hold on;
%     %         set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear')
%     %         axis tight
%     %         title(geneid_all(xi_rnd(i)));
%     %     end
%     %     close(hh)
% end
% 
% % % % % % % % % % % % % % % % % 
% save('Glut1_gene_cluster_timepoint_scatter','geneid_all','T_cells_tmp','data_orig_all_sorted','fc_time_sorted','fc_time_uni','-v7.3');
% 
% gn = 'ecel1';
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
% %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)','.'); hold on;
%         t_ed(i) = median(data_orig_all_sorted(g,gr2));
%         t_av(i) = mean(data_orig_all_sorted(g,gr2));
%         t_75(i) = prctile(data_orig_all_sorted(g,gr2),75);
%         f90(i) = sum(data_orig_all_sorted(g,gr2)>=p90);
%         pp90(i) = binocdf(f90(i),length(gr2),0.1,'upper');
%     end
%     plot([0,6],p90*[1,1],'--k');
% %     plot(t_ed,'-sk'); hold on;
% %     plot(t_av,'-or'); hold on;
% %     plot(t_75,'-dg'); hold on;
%     axis tight
%     yl = get(gca,'ylim');
%     set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear','ylim',[-0.5,p99+1],'YScale',yscale)
%     text([1:5],yl(2)*0.9*ones(1,5),cellfun(@(x) num2str(x,2), m2c(-log10(pp90)),'UniformOutput',0),'fontsize',6);
%     title([gn,',c=',num2str(c)])
% end
% 
% % % % % % % % % for all the data
% gn = 'bdnf';
% g = find(strcmpi(gn,geneid_all));
% 
% figure;
% set(gcf,'color','w','position',[20,20,1000,900]);
% [ha, pos] = tight_subplot(2, 1, [0.05,0.05], [0.05,0.05], [0.05,0.05]);
% t_ed = zeros(5,1);
% t_av = zeros(5,1);
% t_75 = zeros(5,1);
% set(gcf,'CurrentAxes',ha(1))
% f = zeros(5,max(data_orig_all_sorted(g,:))+1);
% for i=[1:5]
%     gr2 = find( fc_time_sorted==fc_time_uni(i));%find(fc_time_sorted==fc_time_uni(i));%
%     set(gcf,'CurrentAxes',ha(2))
%     plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)','.'); hold on;
%     t_ed(i) = median(data_orig_all_sorted(g,gr2));
%     t_av(i) = mean(data_orig_all_sorted(g,gr2));
%     t_75(i) = prctile(data_orig_all_sorted(g,gr2),75);
%     [f(i,:)] = histcounts(data_orig_all_sorted(g,gr2),[0:(max(data_orig_all_sorted(g,:))+1)]);
%     set(gcf,'CurrentAxes',ha(1))
% %     plot(f,'s','markersize',16); hold on;
% end
% set(gcf,'CurrentAxes',ha(1))
% bar(f')
% legend({'0','2','8','24','recall'})
% set(gca,'yscale','linear')
% set(gcf,'CurrentAxes',ha(2))
% plot(t_ed,'-sk'); hold on;
% plot(t_av,'-or'); hold on;
% plot(t_75,'-dg'); hold on;
% set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear')
% axis tight
% title([gn])
% 
% % % % % % % % % % % % % % % % % % % 
% top_g = 150;
% gr1 = find(T_cells_tmp(:,1)>=3 & T_cells_tmp(:,1)<=13);
% gr2 = find(T_cells_tmp(:,1)>=14 & T_cells_tmp(:,1)<=17);
% x1 = mean(log2(data_orig_all_sorted(:,gr1)+1),2);
% x2 = mean(log2(data_orig_all_sorted(:,gr2)+1),2);
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
% title(['Sex comparison (average)'])
% xlabel(['mean (Females)'])
% ylabel(['mean (Males)'])
% axis tight
% 
% top_g = 250;
% x1 = mean(data_orig_all_sorted(:,gr1)>0,2);
% x2 = mean(data_orig_all_sorted(:,gr2)>0,2);
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
% title(['Sex comparison (precentage)'])
% xlabel(['mean>0 (Females)'])
% ylabel(['mean>0 (Males)'])
% axis tight
% if savefig_flag==1
%     savefig(gcf,['scatter_cluster_',num2str(c1),'_vs_',num2str(c2),'_',date,'.fig'])
% end
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% markers = {'Sema4c' ,'Sema6d' ,'Sema4a' ,'Sema6c' ,'Sema3d' ,'Sema3a' ,'Sema3e' ,'Sema3c' ...
%     ,'Sema4f' ,'Sema4b' ,'Sema3g' ,'Sema7a' ,'Sema3b' ,'Sema3f' ,'Sema4d' ,'Sema5a' ,'Sema5b' ...
%     ,'Sema6b' ,'Sema6a' ,'Sema4g' ,'Creb1' ,'Creb3l1' ,'Creb3l4' ,'Creb3' ,'Creb3l2' ,'Creb5' ...
%     ,'Crebl2' ,'Crebzf' ,'Creb3l3' ,'Crebbp' ,'Crebrf' ,'Nrxn3' ,'Nrxn1' ,'Nrxn2' ,'Nlgn3' ,'Nlgn1' ...
%     ,'Nlgn2' ,'Nxph2' ,'Nxpe5' ,'Nxph1' ,'Nxph4' ,'Nxpe2' ,'Nxpe4' ,'Nxph3' ,'Nxpe3'};
% % markers = IEG_markers;
% t1 = [];
% t2 = [];
% for ii=1:length(markers)
%     tmpmark = markers{ii};
%     marker = data_orig_all_sorted(strcmpi(geneid_all,tmpmark),:);
%     marker_percent = zeros(1, length(T_cells_tmp_uni));
%     p_marker = zeros(1, length(T_cells_tmp_uni));
%     M = length(marker);
%     K = sum(marker>0);
%     for j=1:length(T_cells_tmp_uni)
%         c1 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j));
%         c2 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j) & marker'>0);
%         marker_percent(j) = 100*c2/c1;
%         p_marker(j) = hygecdf(c2,M,K,c1,'upper');
%     end
%     marker_percent(isnan(marker_percent)) = 0;
%     eval([markers{ii},'_percent = marker_percent;'])
%     eval([markers{ii},'_p = p_marker;']);
%     eval(['t1 = [t1;',markers{ii},'_percent];'])
%     eval(['t2 = [t2;',markers{ii},'_p];'])
% end
% figure;
% set(gcf,'color','w','position',[20,20,900,1200])
% [ha, pos] = tight_subplot(1, 2, [0.05,0.00], [0.05,0.05], [0.05,0.02]);
% axes(ha(1));
% imagesc(t1');
% colormap('summer');
% set(gca,'ytick',[1:length(idxuni)],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
% axes(ha(2));
% imagesc(-log10(t2)',[0,10]);
% colormap('summer');
% set(gca,'ytick',[],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
% linkaxes([ha(1),ha(2)],'y');
% 
% 
% % % % % % % % % % % % % % % % % 
% gn = 'dlk1';
% g = find(strcmpi(gn,geneid_all));
% c = 26;
% figure;
% set(gcf,'color','w','position',[20,20,1200,1000]);
% [ha, pos] = tight_subplot(4, ceil(length(unique(T_cells_tmp(:,2)))/4), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% for k=1:length(unique(T_cells_tmp(:,2)))
%     axes(ha(k))
%     c=k;
%     for i=[1:5]
%         gr2 = find(T_cells_tmp(:,2)==c & fc_time_sorted==fc_time_uni(i));%find(fc_time_sorted==fc_time_uni(i));%
%         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)','.'); hold on;
%         plot(i,median(data_orig_all_sorted(g,gr2)),'sk'); hold on;
%     end
%     set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear')
%     axis tight
%     title([gn,',c=',num2str(c)])
% end
% 
% % % % % % % % % % % 
% table3 = [cellid_sorted,m2c(T_cells_tmp(:,1))];
% saveCellFile(table3,['cellid_cluster_FC_Glut1_preclean_',date,'.txt']);
% 
% % % % % % % % % % % % % % % % % % % 
% per_cluster_insample = zeros(length(sample_uni),length(T_cells_tmp_uni));
% for i=1:length(sample_uni)
%     in = find(strcmpi(sample_sorted,sample_uni{i}));
%     per_cluster_insample(i,:) = histcounts(T_cells_tmp(in,1),[1:length(T_cells_tmp_uni)+1]);
% end
% per_cluster_insample = per_cluster_insample./repmat(sum(per_cluster_insample,2),1,length(T_cells_tmp_uni));
% % condition_sorted_uni = unique(condition_sorted);
% colorvec = distinguishable_colors(length(fc_time_uni));
% figure;
% set(gcf,'color','w','position',[20,20,1200,1000]);
% [ha, pos] = tight_subplot(8, ceil(length(T_cells_tmp_uni)/8), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% for j=1:length(T_cells_tmp_uni)
%     axes(ha(j))
%     for i=1:length(fc_time_uni)
%         s = unique(sample_sorted(fc_time_sorted==fc_time_uni(i)));
%         [~,loc] = ismember(s,sample_uni);
%         t = per_cluster_insample(loc,j);
%         plot(i,100*t,'.','color',colorvec(i,:)); hold on;
%         plot(i,100*median(t),'sk','markersize',12)
%     end
%     box off
%     axis tight
%     yl = get(gca,'ylim');
%     set(ha(j),'xlim',[0.5,length(fc_time_uni)+0.5])
%     title(num2str(j))
% end
% % % % % % % % % % % % % % % % % % % 
% per_cluster_insample2 = zeros(length(sample_uni),length(unique(T_cells_tmp(:,2))));
% for i=1:length(sample_uni)
%     in = find(strcmpi(sample_sorted,sample_uni{i}));
%     per_cluster_insample2(i,:) = histcounts(T_cells_tmp(in,2),[1:length(unique(T_cells_tmp(:,2)))+1]);
% end
% per_cluster_insample2 = per_cluster_insample2./repmat(sum(per_cluster_insample2,2),1,length(unique(T_cells_tmp(:,2))));
% % condition_sorted_uni = unique(condition_sorted);
% colorvec = distinguishable_colors(length(fc_time_uni));
% figure;
% set(gcf,'color','w','position',[20,20,1200,1000]);
% [ha, pos] = tight_subplot(4, ceil(length(unique(T_cells_tmp(:,2)))/4), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% for j=1:length(unique(T_cells_tmp(:,2)))
%     axes(ha(j))
%     for i=1:length(fc_time_uni)
%         s = unique(sample_sorted(fc_time_sorted==fc_time_uni(i)));
%         [~,loc] = ismember(s,sample_uni);
%         t = per_cluster_insample(loc,j);
%         plot(i,100*t,'.','color',colorvec(i,:)); hold on;
%         plot(i,100*median(t),'sk','markersize',12)
%     end
%     box off
%     axis tight
%     yl = get(gca,'ylim');
%     set(ha(j),'xlim',[0.5,length(fc_time_uni)+0.5])
%     title(num2str(j))
% end
% 
% % % % % % % % % % % % % 
% % in = find(sum(data_orig_all_sorted>0,2)>100 & sum(data_orig_all_sorted>0,2)<length(data_orig_all_sorted(1,:))*0.5 & ~ismember(geneid_all(:,1),sex_genes)  ....
% %     & ~ismember(geneid_all(:,1),[batch_genes]));
% % tic
% % idx = clusterdata(data_orig_all_sorted(in,:),'distance','correlation','linkage','ward','maxclust',100);
% % toc
% % idxuni = unique(idx);
% % % n = histcounts(idx,length(idxuni));
% % d = corr_mat(data_orig_all_sorted(in,:)');
% % dsum = sum(d>0.2,2);
% % 
% % [~,xi] = sort(idx);
% % figure;
% % set(gcf,'color','w','position',[20,20,900,800])
% % h1 = axes('position',[0.2,0.1,0.7,0.85]);
% % imagesc(d(xi,xi),[0,0.3]);
% % set(gca,'ytick',[1:length(xi)],'YTickLabel',geneid_all(in(xi)))
% % freezeColors(h1)
% % h2 = axes('position',[0.91,0.1,0.04,0.85]);
% % imagesc(idx(xi));
% % colormap('prism');
% % h3 = axes('position',[0.95,0.1,0.05,0.85]);
% % plot(sum(d(xi,xi)>0.2,2),[1:length(xi)],'.');
% % set(h3,'ydir','reverse');
% % linkaxes([h1,h2,h3],'y')
% % % % % % % % % % % % % % 
% % options = statset('MaxIter',1000);
% % tmp = log2(data_orig_all_sorted(in,:)+1);
% % tmp = tmp./repmat(sum(tmp,2),1,length(tmp(1,:)));
% % tmp = cent_norm(tmp')';
% % genes_xy = tsne(tmp,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',50,.....
% %     'Perplexity',50,'Standardize',true,'LearnRate',500,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
% % figure;
% % set(gcf,'color','w','position',[20,20,900,800])
% % plot(genes_xy(:,1),genes_xy(:,2),'.'); axis tight; axis off
% % 
% % MinPts = 30;
% % eps_prc = 90;
% % [idx, isnoise] = dbscan_epsprc_mipts(genes_xy,eps_prc,MinPts);
% % idxuni = unique(idx);
% % 
% % colors = distinguishable_colors(length(unique(idx))+1);
% % figure;
% % set(gcf,'color','w','position',[20,20,900,800])
% % for i=unique(idx)'
% %     if i>=0
% %         ii=find(idx==i); h=plot(genes_xy(ii,1),genes_xy(ii,2),'o','color',colors(i+1,:),'markersize',3); hold on;
% %     elseif i>0
% %         ii=find(idx==i); h=plot(genes_xy(ii,1),genes_xy(ii,2),'.','color',colors(i+1,:),'markersize',4); hold on;
% %     end
% % end
% % % for i=idxuni'
% % %     if i>=0
% % %         in = idx==i;
% % %         ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),num2str(i));
% % %         set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',6)
% % %     end
% % % end
% % axis tight;
% % axis equal
% % axis off
% % 
% % title(['MinPts=',num2str(MinPts),', peplexity=',num2str(perplexity),', epsprc=',num2str(eps_prc),',#C=',num2str(max(idx)),',#out=',num2str(sum(idx==0))],'fontsize',8);
% % % 