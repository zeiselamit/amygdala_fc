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
sam_folder = {'out10X08_1_210219','out10X08_2_190709','out10X10_1_210219','out10X10_2_210219'....
    ,'out10X18_1_210219','out10X18_2_210219','out10X19_1_210219','out10X23_1_210219','out10X23_2_210219'...
    ,'out10X23_3_210219','out10X41_1_210220','out10X41_2_200318','out10X42_1_200319','out10X42_2_210220',...
    'out10X45_1_200319','out10X45_2_210220','out10X45_3_200319','out10X46_1_210220','out10X46_2_200319',....
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
    | strcmpi(sample,'75-2') | strcmpi(sample,'75-3') | strcmpi(sample,'76-1') |....
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
% load afterloading_QC3000_FC_06-Oct-2021
tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);

figure('position',[100,100,800,800],'color','w');
subplot(2,2,1);
[f,xi] = ksdensity(tot_mol);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_mol,98)]);
xlabel('total mol')
subplot(2,2,3);
[f,xi] = ecdf(tot_mol);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_mol,98)]);
xlabel('total mol')
subplot(2,2,2);
[f,xi] = ksdensity(tot_genes);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_genes,98)]);
xlabel('total genes')
subplot(2,2,4);
[f,xi] = ecdf(tot_mol);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_genes,98)]);
xlabel('total genes')

data = normalize(data,'norm');
data = round(data./repmat(sum(data),length(data(:,1)),1)*30e3);
[median(sum(data(:,batch_flag==1))),median(sum(data(:,batch_flag==2))),median(sum(data(:,batch_flag==3)))]
[median(sum(data(:,batch_flag==1)>0)),median(sum(data(:,batch_flag==2)>0)),median(sum(data(:,batch_flag==3)>0))]

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};
in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),sex_genes));
corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);

data_orig_all = data;
geneid_all = geneid;
data = data(in(corr_filt),:);
geneid = geneid(in(corr_filt));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
moldata = data;

data_tsne = cent_norm([(log2(moldata+1))]);
initial_dims = length(corr_filt);
[prj,~,D,~,~] = pca_wis(data_tsne',initial_dims);
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
% % % % % % % % % % % % % %  harmony py
s = prj;
batchid = batch_flag;
usepylib = 1;
[sout]=harmonypy(s,batchid,usepylib);
prj = sout;
% % % % % % % % % % % % % % 
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
% sample_uni = {'8-1','10-1','18-1','19-1','23-1','23-3','41-1','41-2','42-1','45-1','45-2','45-3','46-1','46-2','46-3'};
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
legend(sample_uni)
% % % % % % % % % % % % 
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
legend({'ctrl','2hr','8hr','24hr','recall'})
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
legend({'Amy','Piri'})
if savefig_flag==1
    % eval(['export_fig tsne_AmyPiri_FC_by_sampleannot_','_',date,'.pdf']);
    savefig(gcf,['tsne_FC_0_2_8_24_Neurons_by_annot_perplexity_',num2str(perplexity),'_PCAdim',num2str(initial_dims),'_',date,'.fig'])
end
% % % % % 
% Z = linkage(mapped_xy,'ward','euclidean');
% c = cluster(Z,'Maxclust',120);
% idx = c;
% clustering with dbscan
MinPts = 40;
eps_prc = 85;
[idx, isnoise] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);
idxuni = unique(idx);

colors = distinguishable_colors(length(unique(idx))+1);
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=unique(idx)'
    if i>=0
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

title(['peplexity=',num2str(perplexity),', epsprc=',num2str(eps_prc),',#C=',num2str(max(idx)),',#out=',num2str(sum(idx==0))],'fontsize',8);

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
% for i=1:length(idxuni)
%     i
%     ind = find(idx==i);
%     if length(ind)>20
%         tmp1d = fast_tsne((data_tsne(:,ind))', no_dims, initial_dims, perplexity,theta, rand_seed);
%         [~,xitmp] = sort(tmp1d);
%         xi(ind) = xi(ind((xitmp)));
%     end
% end

data_sorted_all = data_sorted_all(:,xi);
data_orig_all_sorted = data_orig_all_sorted(:,xi);
cellid_sorted = cellid_sorted((xi));
sample_sorted = sample_sorted((xi));
mapped_xy = mapped_xy(xi,:);
fc_time_sorted = fc_time_sorted(xi);
amy_flag_sorted = amy_flag_sorted(xi);
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
% meangr_mat1(loc,:) = [];

% meangr_mat1 = cent_norm(meangr_mat1(:,leaforder));
[prj,m,D,V,Q] = pca_wis(meangr_mat1',initial_dims);
Zpca = linkage(prj,'ward','correlation');
Dpca = pdist(prj,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
figure;
set(gcf,'position',[100,100,1000,1000],'color','w')
axes('position',[0.03,0.03,0.3,0.93])
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
axes('position',[0.35,0.03,0.63,0.93])
x=squareform(Dpca); imagesc(x(leaforder_pca,leaforder_pca));
colormap('summer')
set(gca,'ytick',[1:length(leaforder_pca)],'xtick',[],'fontsize',8,'ydir','normal')

if savefig_flag==1
    savefig(gcf,['tree_FC_',date,'.fig'])    
end
leaforder = leaforder_pca;

T_cells_tmp_new = zeros(size(idx));
for i=1:length(leaforder)
    T_cells_tmp_new(idx==idxuni(leaforder(i))) = i;
end
idxuni_new = unique(T_cells_tmp_new);
[~,xi] = sort(T_cells_tmp_new);
T_cells_tmp_new = T_cells_tmp_new(xi);
data_sorted_all = data_sorted_all(:,xi);
data_orig_all_sorted = data_orig_all_sorted(:,xi);
cellid_sorted = cellid_sorted(xi);
% % tissue = tissue(xi);
mapped_xy = mapped_xy(xi,:);
fc_time_sorted = fc_time_sorted(xi);
amy_flag_sorted = amy_flag_sorted(xi);
piri_flag_sorted = piri_flag_sorted(xi);
cells_bor_2 = find(diff(T_cells_tmp_new)>0)+1;
sample_sorted = sample_sorted(xi);
prj_sorted = prj_sorted(xi,:);
T_cells_tmp = T_cells_tmp_new;
T_cells_tmp_uni = unique(T_cells_tmp);
% % % %
idx = T_cells_tmp_new;
idxuni = idxuni_new;


colors = distinguishable_colors(length(unique(idx))+1);
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=unique(idx)'
    if i==-1
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',[0,0,1],'markersize',3); hold on;
    else
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i+1,:),'markersize',5); hold on;
    end
end
for i=idxuni'
    in = idx==i;
    ht = text(median(mapped_xy(in,1)),median(mapped_xy(in,2)),num2str(i));
    set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',8)
end
axis tight;
axis equal
axis off
title(['MinPts=',num2str(MinPts),', epsprc=',num2str(eps_prc),',#C=',num2str(max(idx)),',#out=',num2str(sum(idx==0))],'fontsize',8);
if savefig_flag==1
    savefig(gcf,['tsne_FC_0_2_8_24_Neurons_by_cluster_',date,'.fig'])
    % eval(['export_fig tsne_AmyPiri_FC_perplexity_',num2str(perplexity),'_Ngenes=',num2str(length(data_tsne(:,1))),'_NPCA=',num2str(initial_dims),'_',date,'.pdf']);
end
if savefig_pdf==1
    eval(['export_fig tsne_FC_0_2_8_24_Neurons_by_cluster_',date,'.pdf']);
end

% % % % % % % % % % % % % 
[ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp,data_sorted_all,5);
% % % % % % % % % 
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

% sample_uni = {'8-1','10-1','18-1','19-1','23-1','23-3','41-1','41-2','42-1','45-1','45-2','45-3','46-1','46-2','46-3'};
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

ax5 = axes('position',[0.1,0.83,0.88,0.135]);
imagesc(~samples_num); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1:length(sample_uni)],'yticklabel',sample_uni,'fontsize',6);

gad2 = data_orig_all_sorted(strcmpi(geneid_all,'Gad2'),:);
slc17a7 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a7'),:);
slc17a6 = data_orig_all_sorted(strcmpi(geneid_all,'Slc17a6'),:);

ax6 = axes('position',[0.1,0.965,0.88,0.01]);
imagesc(~gad2); hold on;
axes(ax6)
colormap('gray');
freezeColors(ax6);
set(ax6,'xtick',[],'ytick',[1],'yticklabel','Gad2','fontsize',5);
ax7 = axes('position',[0.1,0.965+0.01,0.88,0.01]);
imagesc(~slc17a6); hold on;
colormap('gray');
freezeColors(ax7);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a6','fontsize',6);
ax8 = axes('position',[0.1,0.965+2*0.01,0.88,0.01]);
imagesc(~slc17a7); hold on;
colormap('gray');
freezeColors(ax8);
set(gca,'xtick',[],'ytick',[1],'yticklabel','Slc17a7','fontsize',6);

linkaxes([ax1,ax2_1,ax2_2,ax2_3,ax2_4,ax2_5,ax3,ax4,ax5,ax6,ax7,ax8],'x');
if savefig_flag==1
    savefig(gcf,['markertable_all_onNeurons_Rep4_',date,'.fig'])
    % eval(['export_fig markertable_AmyPiri_FC_',date,'.pdf']);
end
if savefig_pdf==1
    eval(['export_fig markertable_FC_0_2_8_24_Neurons_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plxnb3','Cldn11'....
    ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc25a18','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2','Ttr'};

[~,loc] = ismember(exclude_markers,geneid_all);
nonneuro = sum(data_orig_all_sorted(loc,:));

gabaglut = zeros(size(T_cells_tmp_uni));
gabaglut_sc = zeros(size(T_cells_tmp));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    %     tmp = [mean(gad2(:,T_cells_tmp==jjj)>0);mean(slc17a7(:,T_cells_tmp==jjj)>0);....
    %         mean(slc17a6(:,T_cells_tmp==jjj)>0);mean(nonneuro(:,T_cells_tmp==jjj)>0)];
    tmp = [mean(gad2(:,T_cells_tmp==jjj));mean(slc17a7(:,T_cells_tmp==jjj));....
        mean(slc17a6(:,T_cells_tmp==jjj));mean(nonneuro(:,T_cells_tmp==jjj))];
    tmpsort = [tmp(1),max(tmp(2:3)),tmp(4)];
    tmpsort= sort(tmpsort,'descend');
    if tmpsort(1)>2*tmpsort(2)
        [~,gabaglut(jjj)] = max(tmp);
    else
        gabaglut(jjj) = 5;
    end
    gabaglut_sc(T_cells_tmp==jjj) = gabaglut(jjj);
end
figure;
set(gcf,'color','w','position',[20,20,800,800])
colors = distinguishable_colors(length(unique(gabaglut_sc)));
for idx=unique(gabaglut_sc)'
    ii=find(gabaglut_sc==idx); 
    h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(idx,:),'markersize',3); hold on;
end
axis tight
axis equal
axis off
title('GABA/Glut1/Glut2')
legend('GABA','Glut1','Glut2','non-neurons','doublets')
% cid = cellfun(@(x,y) [x,'_',y], cellid_sorted, sample_sorted,'UniformOutput',0);
figure('color','w','position',[20,20,800,800])
markergene = nonneuro;
inpos = markergene>0;
tmpthlow = prctile(markergene(markergene>0),1);
tmpthhigh = prctile(markergene(markergene>0),90);
markergene(markergene>tmpthhigh) = tmpthhigh;
markergene(markergene<tmpthlow) = tmpthlow;
c_rgb = [1,0,0];rand([1,3]);
markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
    interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
    ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
title('nonneurons sum');
axis tight
axis equal
axis off

table1 = [cellid_sorted,m2c(gabaglut_sc)];
saveCellFile(table1,['cellid_gaba_glut_FC_QC3000_',date,'.txt']);
% % % % % % % % % % % % % %
list = {'Snap25','Stmn2','Gad2','Slc32a1','Slc17a7','Slc17a6','Sst','Sim1','Foxj1','Pdgfra','Mog','C1qc','Flt1','Cldn5','Aqp4','Plp1'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for i=1:length(list)
    genePlot = list{i};
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
    axes(ha(i));
    scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    axis tight
    axis equal
    axis off
end


toc