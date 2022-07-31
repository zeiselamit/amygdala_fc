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
    
data = normalize(data,'norm');
data = ceil(data./repmat(sum(data),length(data(:,1)),1)*10e3);
[median(sum(data(:,batch_flag==1))),median(sum(data(:,batch_flag==2))),median(sum(data(:,batch_flag==3)))]
[median(sum(data(:,batch_flag==1)>0)),median(sum(data(:,batch_flag==2)>0)),median(sum(data(:,batch_flag==3)>0))]

IEG_markers = loadCellFile('IEG_cluster.txt');
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};%,'Snca','Gm9843','Gm12027','Gm12031'
in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.8 & ~ismember(geneid(:,1),sex_genes)  ....
    & ~ismember(geneid(:,1),[IEG_markers]));

corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);
g = unique([geneid(in(corr_filt));geneselection_oligos(:,1);geneselection_astro(:,1);geneselection_vascular(:,1);geneselection_immune(:,1)]);
[~,corr_filt] = ismember(g,geneid);
corr_filt(sum(data(corr_filt,:)>0,2)<10) = [];
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
% % % % % % % 
% % % % % % % % % % % % % % % % % % 
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

perplexity = 100;%median(optk);

options = statset('MaxIter',1000);
mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',5);
toc
% 


% % % % % % % % % % % % % % % % % % % 
clusteruni = unique(cellid_clusters(:,4));
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,4),clusteruni{i}));
    T_cells_tmp(ind,1) = i;
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp(:,1));
saveCellFile([clusteruni,m2c(clustersize)],'clustername_all_nonneurons.txt')
% % % % % % % % % % % % % % % % % % % 
clustername_color_hannah = loadCellFile_turbo('non-neuronal-hexa.txt',1);
clustername_color_hannah = clustername_color_hannah(2:end,2:3);
[~,loc] = ismember(upper(clusteruni),upper(regexprep(clustername_color_hannah(:,1),'-','_')));
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
    savefig(gcf,['tsne_final_step3_FC_nonneuronal_by_cluster_',date,'.fig'])    
end
if savefig_pdf==1
    eval(['export_fig tsne_final_step3_FC_nonneuronal_by_cluster_',date,'.pdf']);
end
% % 

figure;
set(gcf,'color','w','position',[20,20,1350,500])
subplot(1,5,1);
dscatter(mapped_xy(fc_time==0,1),mapped_xy(fc_time==0,2));
axis tight;
axis equal
axis off
title('Control')
subplot(1,5,2);
dscatter(mapped_xy(fc_time==2,1),mapped_xy(fc_time==2,2));
axis tight;
axis equal
axis off
title('2hr')
subplot(1,5,3);
dscatter(mapped_xy(fc_time==8,1),mapped_xy(fc_time==8,2));
axis tight;
axis equal
axis off
title('8hr')
subplot(1,5,4);
dscatter(mapped_xy(fc_time==24,1),mapped_xy(fc_time==24,2));
axis tight;
axis equal
axis off
title('24hr')
subplot(1,5,5);
dscatter(mapped_xy(fc_time==28,1),mapped_xy(fc_time==28,2));
axis tight;
axis equal
axis off
title('recall')
if savefig_pdf==1
    eval(['export_fig tsne_density_bytimepoint_final_step3_FC_nonneuronal_',date,'.pdf']);
end

% % % % % % % % % % % % % % % % % 

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
save(['all_nonneurons_cluster_ave'],'meangr_mat_all','clusteruni','geneid_all','geneid');

[prj,m,D,V,Q] = pca_wis(meangr_mat1',10);
% prj = prj_mean';
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
h2 = axes('position',[0.37,0.03,0.6,0.93]);
x = squareform(Dpca); imagesc(x(leaforder_pca,leaforder_pca));
colormap('summer')
set(gca,'ytick',[1:length(clusteruni)],'YTickLabel',regexprep(clusteruni(leaforder_pca),'_','-'),'xtick',[],'fontsize',8,'ydir','normal')
linkaxes([h1,h2],'y');
if savefig_pdf==1
    eval(['export_fig tree_corrmat_nonneurons_per_cluster_',date,'.pdf']);
end
fc_time_uni = unique(fc_time);
%%
data = data_orig_all;
geneid = geneid_all;
%%
gn = 'fth1';
g = find(strcmpi(gn,geneid));
figure;
set(gcf,'color','w','position',[20,20,1200,1000],'Name',gn);
b = 1;
[ha, pos] = tight_subplot(3, ceil(length(unique(T_cells_tmp(:,b)))/3), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
logflag = 0;
p99 = prctile(data(g,:),99.99);
colorvec = distinguishable_colors(3);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    c=k;
    t_ed = zeros(5,1);
    t_av = zeros(5,1);
    t_75 = zeros(5,1);
    p90 = prctile(data(g,T_cells_tmp(:,b)==c),90);
    for i=[1:5]
        gr2 = find(T_cells_tmp(:,b)==c & fc_time==fc_time_uni(i));%find(fc_time_sorted==fc_time_uni(i));%
        bfgr2 = batch_flag(gr2);
        y = (data(g,gr2));
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
            x1 = i + fi'.*(0.9*rand(length(gr2),1)-0.45);
            y1 = 0.5*rand(length(gr2),1)-0.1+y';
            for jj=1:3
                v = bfgr2==jj;
            plot(x1(v),y1(v) ,'.','color',colorvec(jj,:)); hold on;
            end
        else
            plot(i+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
        end
        
        %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
        t_ed(i) = median(data(g,gr2));
        t_av(i) = mean(data(g,gr2));
        t_75(i) = prctile(data(g,gr2),75);
        f90(i) = sum(data(g,gr2)>p90);
        pp90(i) = binocdf(f90(i),length(gr2),0.1,'upper');
    end
    plot([0,6],p90*[1,1],'--k');
    %     plot(t_ed,'-sk'); hold on;
    %     plot(t_av,'-or'); hold on;
    %     plot(t_75,'-dg'); hold on;
    axis tight
    yl = get(gca,'ylim');
    if k>=(length(T_cells_tmp_uni)-ceil(length(unique(T_cells_tmp(:,b)))/8))
        set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear','ylim',[-0.5,p99+1],'YScale',yscale)
    else
        set(gca,'xtick',[1:5],'XTickLabel',cell(1,5),'yscale','linear','ylim',[-0.5,p99+1],'YScale',yscale)
    end
    for ii=1:5
        if -log10(pp90(ii))>1
            fs = 10;
        else
            fs = 6;
        end
        text(ii,p99*0.7,cellfun(@(x) num2str(x,2), m2c(-log10(pp90(ii))),'UniformOutput',0),'fontsize',fs);
    end
    title(regexprep(clusteruni{c},'_','-'));%[gn,',c=',num2str(c)])
end

if 1==100
    %     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
    eval(['export_fig ',gn,'_violin_FC_nonneurons_',date,'.pdf']);    
end
%%

%%
nsim = 100;
t_sample_uni = [0,0,0,0,0,0,0,2,2,2,2,2,8,8,8,8,24,24,24,24,28,28,28];
subset_p = [randi(7,[nsim,6]),randi([8,12],[nsim,4]),randi([13,16],[nsim,3]),randi([17,20],[nsim,3]),randi([21,23],[nsim,2])];
pbino_cell = cell(nsim,1);
generank = zeros(length(geneid),nsim);
for jjj=1:nsim
    jjj
    s = ismember(sample, sample_uni(subset_p(jjj,:)));
    k = 0;
    b = 1;
    pbino_mat = ones(length(geneid),5*length(T_cells_tmp_uni));
    for c1=1:length(unique(T_cells_tmp(:,b)))
        c1;
        fc_time_uni = [0,2,8,24,28];
        p90gene = prctile(data(:,s & T_cells_tmp(:,b)==c1),90,2);
        for i=[1,2,3,4]
            k = k+1;
            gr2 = find(s & T_cells_tmp(:,b)==c1 & fc_time==fc_time_uni(i));
            x90 = sum(data(:,gr2)>repmat(p90gene,1,length(gr2)),2);
            x90(p90gene==0) = nan;
            pbino = binocdf(x90,length(gr2),0.1,'upper');
            pbino(isnan(pbino)) = 1;
            pbino_mat(:,k) = pbino;
        end
    end
    x = sum(pbino_mat(:,:)<1e-2,2);
    [~,xi] = sort(x,'descend');
    generank(xi,jjj) = [1:length(xi)]';
    pbino_cell{jjj} = pbino_mat;
end


f = sum(generank<300,2) - sum(generank>1000,2);
f = sum(generank,2) ;
f = mean(generank<500,2);
[~,xi] = sort(f,'descend');
[geneid(xi(1:100)),m2c(f(xi(1:100)))];
hitgenes = [geneid(f>0.5),m2c(f(f>0.5))]
saveCellFile([hitgenes],['hitgenes_nonneurons_',date,'.txt'])

in = find(sum(data>0,2)>100 & sum(data>0,2)<length(data(1,:))*0.9);
k = 0;
b = 1;
pbino_mat = ones(length(geneid),5*length(T_cells_tmp_uni));
pbino_mat_rnd = ones(length(geneid),5*length(T_cells_tmp_uni));
fc_time_rnd = fc_time;
p90_mat = zeros(size(pbino_mat));
for c1=1:length(unique(T_cells_tmp(:,b)))
    c1
    fc_time_uni = [0,2,8,24,28];
    t = fc_time_rnd(T_cells_tmp(:,b)==c1);
    fc_time_rnd(T_cells_tmp(:,b)==c1) = t(randperm(length(t)));
    p90gene = prctile(data(:,T_cells_tmp(:,b)==c1),50,2);
    for i=[1,2,3,4,5]
        k = k+1
        gr2 = find(T_cells_tmp(:,b)==c1 & fc_time==fc_time_uni(i));
        gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(i));
        x90 = sum(data(:,gr2)>repmat(p90gene,1,length(gr2)),2);
        x90(p90gene==0) = nan;
        pbino = binocdf(x90,length(gr2),0.5,'upper');
        pbino(isnan(pbino)) = 1;
        pbino_mat(:,k) = pbino;
        x90_rnd = sum(data(:,gr2_rnd)>repmat(p90gene,1,length(gr2_rnd)),2);
        x90_rnd(p90gene==0) = nan;
        pbino_rnd = binocdf(x90_rnd,length(gr2_rnd),0.1,'upper');
        pbino_rnd(isnan(pbino_rnd)) = 1;
        pbino_mat_rnd(:,k) = pbino_rnd;
        p90_mat(:,k) = prctile(data(:,gr2),90,2);
    end
end
save(['pbino_nonneuron_',date],'pbino_mat','geneid','clusteruni')
save(['p90mat_nonneuron_',date],'p90_mat','geneid','clusteruni')
x = sum(pbino_mat<1e-2,2);
[~,xi] = sort(x,'descend');
[geneid(xi(1:50)),m2c(x(xi(1:50)))]

[~,inhit] = ismember(hitgenes(:,1),geneid);
m = -log10(pbino_mat(inhit,:));
z = linkage(m,'ward','correlation');
d = corr_mat(m');
leaforder = optimalleaforder(z,squareform(1-d,'tovector'));
tp_order = reshape([1:5*length(T_cells_tmp_uni)],5,length(T_cells_tmp_uni))';
tp_order = tp_order(:);
figure('position',[200,200,1200,1000],'color','w');
set(gca,'Position',[0.1,0.03,0.85,0.9])
imagesc((m(leaforder,tp_order)),[0,3]); hold on;
colormap('summer');

freezeColors(gca);
yl = get(gca,'ylim');
for i=1:length(fc_time_uni)-1
    plot(i*length(T_cells_tmp_uni)*[1,1]+0.5,yl,'-k');
end
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',hitgenes(leaforder,1),'fontsize',10,'xtick'...
    ,[length(T_cells_tmp_uni)/2:length(T_cells_tmp_uni):5*length(T_cells_tmp_uni)],'XTickLabel',{'0','2','8','24','recall'})

mnew = m(leaforder,tp_order);
c_axis = [];
for i=1:5
    ind = (i-1)*length(T_cells_tmp_uni)+1:i*length(T_cells_tmp_uni);
    mtmp = mnew(:,ind)';
    z = linkage(mtmp,'ward','correlation');
    d = corr_mat(mtmp');
    leaforder_c = optimalleaforder(z,squareform(1-d,'tovector'));
    mnew(:,ind) = mnew(:,ind(leaforder_c));
    c_axis = [c_axis,leaforder_c];
end
figure('position',[200,200,1200,1000],'color','w');
set(gca,'Position',[0.1,0.03,0.85,0.9])
imagesc((mnew),[0,3]); hold on;
colormap('summer');
freezeColors(gca);
yl = get(gca,'ylim');
for i=1:length(fc_time_uni)-1
    plot(i*length(T_cells_tmp_uni)*[1,1]+0.5,yl,'-k');
end
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',hitgenes(leaforder,1),'fontsize',10,'xtick'...
    ,[length(T_cells_tmp_uni)/2:length(T_cells_tmp_uni):5*length(T_cells_tmp_uni)],'XTickLabel',{'0','2','8','24','recall'})


%%


