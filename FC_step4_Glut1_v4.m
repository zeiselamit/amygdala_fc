tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 1;
savefig_pdf = 1;


load afterloading_QC3000_FC_31-Mar-2022
IEG_markers = loadCellFile('IEG_cluster.txt');
[~,in_ieg] = ismember(IEG_markers,geneid);
tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);
batch_flag(strcmpi(sample,'76-2')) = 3;

cellidclsuter_gaba = loadCellFile('cellidCluster_step3_Glut1_FC_24-Oct-2021.txt');
% cellidclsuter_gaba = [cellidclsuter_gaba, cellfun(@(x,y) ['gaba_',num2str(x),'_',num2str(y)],cellidclsuter_gaba(:,2),cellidclsuter_gaba(:,3),'UniformOutput',0)]; 

gaba_clustername = loadCellFile('ClusterOrder_step3_Glut1_FC_24-Oct-2021.txt');
cluster_info = loadCellFile('Clusters_info_hannah.txt');
cluster_info = cluster_info(58:89,:);

cellidclsuter_gaba = [cellidclsuter_gaba,cell(length(cellidclsuter_gaba),1)];
c = cell2mat(cellidclsuter_gaba(:,2));
for i=1:length(gaba_clustername)
    ind = find(c==gaba_clustername{i});
    cellidclsuter_gaba(ind,end) = repmat(gaba_clustername(i,2),length(ind),1);
end


cellid_clusters = cell(length(cellid),4);
[~,loc]= ismember(cellidclsuter_gaba(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_gaba;

validcells = cell2mat(cellfun(@(x) ~isempty(x), cellid_clusters(:,1),'UniformOutput',0));
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
cellid_clusters = cellid_clusters(validcells,:);
sample = sample(validcells);
fc_time = fc_time(validcells);
amy_flag = amy_flag(validcells);
batch_flag = batch_flag(validcells);

[~,xi] = sort(cell2mat(cellid_clusters(:,2)));
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
mtgenes = contains(geneid,'mt-');
tot_mol = sum(data);
sum_mtgenes = sum(data(mtgenes,:));
figure; loglog(sum_mtgenes,tot_mol,'.');
figure; loglog(sum_mtgenes(batch_flag==1),tot_mol(batch_flag==1),'.r');hold on;
loglog(sum_mtgenes(batch_flag==2),tot_mol(batch_flag==2),'.g');hold on;
loglog(sum_mtgenes(batch_flag==3),tot_mol(batch_flag==3),'.b');hold on;

mtratio = sum_mtgenes./tot_mol;
figure;
[f,xi] = ksdensity(mtratio(batch_flag==1));hold on
plot(xi,f);
[f,xi] = ksdensity(mtratio(batch_flag==2));
plot(xi,f);
[f,xi] = ksdensity(mtratio(batch_flag==3));
plot(xi,f);
legend('b1','b2','b3');

figure;
for i=1:length(sample_uni)
    [f,xi] = ksdensity(mtratio(strcmpi(sample,sample_uni{i})));hold on
    plot(xi,f);
end
legend(sample_uni)



data = normalize(data,'norm');
data = round(data./repmat(sum(data),length(data(:,1)),1)*20e3);
[median(sum(data(:,batch_flag==1))),median(sum(data(:,batch_flag==2))),median(sum(data(:,batch_flag==3)))]
[median(sum(data(:,batch_flag==1)>0)),median(sum(data(:,batch_flag==2)>0)),median(sum(data(:,batch_flag==3)>0))]

% % % % % % % % % % % % % % % % % % % 
clusteruni = gaba_clustername(:,2);
cluster_location = cluster_info(:,6);
cluster_loc_flag = cell(size(cellid));
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,4),clusteruni{i}));
    cluster_loc_flag(ind) = repmat(cluster_location(i),length(ind),1);
    T_cells_tmp(ind) = i;
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp);
fc_time_uni = unique(fc_time);
% % % % % % % % % % % % % % % % % % % 
%%
b = 1;
p50_batch2_cluster_tp = repmat({zeros(length(geneid),3)},1,length(T_cells_tmp_uni));
p90_batch2_cluster_tp = repmat({zeros(length(geneid),3)},1,length(T_cells_tmp_uni));
grsz_batch2 = zeros(length(T_cells_tmp_uni),3);
for k=1:length(T_cells_tmp_uni)
    k
    jj = 0;
    for i=[1,3,4]
        jj = jj+1;
        gr2 = find(T_cells_tmp(:,b)==k & fc_time==fc_time_uni(i) & batch_flag==2);%find(fc_time_sorted==fc_time_uni(i));%
        p50_batch2_cluster_tp{k}(:,jj) = prctile(data(:,gr2),50,2);
        p90_batch2_cluster_tp{k}(:,jj) = prctile(data(:,gr2),90,2);
        grsz_batch2(k,jj) = length(gr2);
    end
end

p50_batch3_cluster_tp = repmat({zeros(length(geneid),3)},1,length(T_cells_tmp_uni));
p90_batch3_cluster_tp = repmat({zeros(length(geneid),3)},1,length(T_cells_tmp_uni));
grsz_batch3 = zeros(length(T_cells_tmp_uni),3);
for k=1:length(T_cells_tmp_uni)
    k
    jj = 0;
    for i=[1,2,5]
        jj = jj+1;
        gr2 = find(T_cells_tmp(:,b)==k & fc_time==fc_time_uni(i) & batch_flag==3);%find(fc_time_sorted==fc_time_uni(i));%
        p50_batch3_cluster_tp{k}(:,jj) = prctile(data(:,gr2),50,2);
        p90_batch3_cluster_tp{k}(:,jj) = prctile(data(:,gr2),90,2);
        grsz_batch3(k,jj) = length(gr2);
    end
end

p50_batch23_cluster_tp = repmat({zeros(length(geneid),5)},1,length(T_cells_tmp_uni));
p90_batch23_cluster_tp = repmat({zeros(length(geneid),5)},1,length(T_cells_tmp_uni));
grsz_batch23 = zeros(length(T_cells_tmp_uni),5);
for k=1:length(T_cells_tmp_uni)
    k
    jj = 0;
    for i=[1,2,3,4,5]
        jj = jj+1;
        gr2 = find(T_cells_tmp(:,b)==k & fc_time==fc_time_uni(i) & batch_flag>=2);%find(fc_time_sorted==fc_time_uni(i));%
        p50_batch23_cluster_tp{k}(:,jj) = prctile(data(:,gr2),50,2);
        p90_batch23_cluster_tp{k}(:,jj) = prctile(data(:,gr2),90,2);
        grsz_batch23(k,jj) = length(gr2);
    end
end

top_g = 20;
topgene_b2 = [];
for k=1:length(T_cells_tmp_uni)
    k
    if min(grsz_batch2(k,:))>30
        tmp = p50_batch2_cluster_tp{k};
        in = find(max(tmp,[],2)>10);
        tmp2 = p90_batch2_cluster_tp{k};
        tmp3 = tmp2-tmp; %(tmp2+1)./(tmp+1);
        fc50 = (max(tmp(in,:),[],2)+1)./(min(tmp(in,:),[],2)+1);
%         fc90 = tmp2(in,2)-tmp2(in,1);%(max(tmp2(in,:),[],2)+1)./(min(tmp2(in,:),[],2)+1);
        fc90 = tmp3(in,3)-tmp3(in,1);
        [~,xi] = sort(fc90,'descend');%
        topgene_b2 = [topgene_b2; [geneid(in(xi(1:top_g))),m2c([fc90(xi(1:top_g)),fc50(xi(1:top_g)),repmat(k,top_g,1)])]];
    end
end
    
%%
gn = 'trh';
g = find(strcmpi(gn,geneid));
figure;
set(gcf,'color','w','position',[20,20,1200,1000],'Name',gn);
b = 1;
if b==2;
    [ha, pos] = tight_subplot(4, ceil(length(unique(T_cells_tmp(:,b)))/4), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
elseif b==1;
    [ha, pos] = tight_subplot(8, ceil(length(unique(T_cells_tmp(:,b)))/8), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
end
logflag = 0;
colorvec = distinguishable_colors(3);
p99 = prctile(data(g,:),99);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    t_90 =[];
    c=k;
    t_ed = zeros(5,1);
    t_av = zeros(5,1);
    t_75 = zeros(5,1);
    p90 = prctile(data(g,T_cells_tmp(:,b)==c),90);
    for i=[1:5]
        gr2 = find(T_cells_tmp(:,b)==c & fc_time==fc_time_uni(i) & batch_flag>=2);%find(fc_time_sorted==fc_time_uni(i));%
        bfgr2 = batch_flag(gr2);
        y = (data(g,gr2));
        if logflag ==1
            y = y+1;
            yscale = 'log';
        else
            yscale = 'linear';
        end
        if length(y)>10
            %             [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
            [f,xi] = ksdensity(y);
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
        t_90(i) = prctile(data(g,gr2),90);
        f90(i) = sum(data(g,gr2)>p90);
        pp90(i) = binocdf(f90(i),length(gr2),0.1,'upper');
    end
    plot([0,6],p90*[1,1],'--k');
    %     plot(t_ed,'-sk'); hold on;
    %     plot(t_av,'-or'); hold on;
    %     plot(t_75,'-dg'); hold on;
        plot(t_90,'-dg'); hold on;
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
    eval(['export_fig ',gn,'_violin_FC_Glut1_',date,'.pdf']);    
end
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
saveCellFile([hitgenes],['hitgenes_glut1_',date,'.txt'])


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
    p90gene = prctile(data(:,T_cells_tmp(:,b)==c1),90,2);
    for i=[1,2,3,4,5]
        k = k+1
        gr2 = find(T_cells_tmp(:,b)==c1 & fc_time==fc_time_uni(i));
        gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(i));
        x90 = sum(data(:,gr2)>repmat(p90gene,1,length(gr2)),2);
        x90(p90gene==0) = nan;
        pbino = binocdf(x90,length(gr2),0.1,'upper');
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
save(['pbino_glut1_',date],'pbino_mat','geneid','clusteruni')
save(['p90mat_glut1_',date],'p90_mat','geneid','clusteruni')
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
%%
in = find(sum(data>0,2)>100 & sum(data>0,2)<length(data(1,:))*0.9);
k = 0;
b = 1;
pbino_mat = ones(length(geneid),5*length(T_cells_tmp_uni));
pbino_mat_rnd = ones(length(geneid),5*length(T_cells_tmp_uni));
fc_time_rnd = fc_time;
for c1=1:length(unique(T_cells_tmp(:,b)))
    c1
    fc_time_uni = [0,2,8,24,28];
    t = fc_time_rnd(T_cells_tmp(:,b)==c1);
    fc_time_rnd(T_cells_tmp(:,b)==c1) = t(randperm(length(t)));
    p90gene = prctile(data(:,T_cells_tmp(:,b)==c1),90,2);
    for i=[1,2,3,4]
        k = k+1
        gr2 = find(T_cells_tmp(:,b)==c1 & fc_time==fc_time_uni(i));
        gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(i));
        x90 = sum(data(:,gr2)>repmat(p90gene,1,length(gr2)),2);
        x90(p90gene==0) = nan;
        pbino = binocdf(x90,length(gr2),0.1,'upper');
        pbino(isnan(pbino)) = 1;
        pbino_mat(:,k) = pbino;
        x90_rnd = sum(data(:,gr2_rnd)>repmat(p90gene,1,length(gr2_rnd)),2);
        x90_rnd(p90gene==0) = nan;
        pbino_rnd = binocdf(x90_rnd,length(gr2_rnd),0.1,'upper');
        pbino_rnd(isnan(pbino_rnd)) = 1;
        pbino_mat_rnd(:,k) = pbino_rnd;
    end
end
x = sum(pbino_mat(in,:)<1e-2,2);
x_rnd = sum(pbino_mat_rnd(in,:)<1e-2,2);
[~,xi] = sort(x,'descend');
g = [geneid(in(xi(1:100))),m2c(x(xi(1:100)))]
m = -log10(pbino_mat(in(xi(1:100)),:)); m = m(:,1:k);
m_ieg = -log10(pbino_mat(in_ieg,:));
percluster_tp = cell(length(T_cells_tmp_uni),1);
for i=1:length(T_cells_tmp_uni)
    [z,zmin] = min(pbino_mat(:,(i-1)*5+1:i*5),[],2);
    [~,zi] = sort(z);
    percluster_tp{i}  = [geneid(zi(1:20)),m2c([fc_time_uni(zmin(zi(1:20)))',-log10(z(zi(1:20)))])];
end
x = cat(1,percluster_tp{:});
[~,xi] = sort(cell2mat(x(:,2)),'ascend');
x(xi,:)

saveCellFile(g,'top300_mostdiff_FC_Glut1.txt')

msum = reshape(sum(m>=2),length(T_cells_tmp_uni),5);
figure('position',[200,200,1200,1000],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.05,0.05], [0.02,0.02], [0.05,0.05]);
for i=1:length(fc_time_uni)
    s = msum(:,i);
    [ssort,sxi] = sort(s,'descend');
    axes(ha(i));
    scatter(ssort,[1:length(T_cells_tmp_uni)],clustersize(sxi));
    set(ha(i),'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',regexprep(clusteruni(sxi),'_','-'),'fontsize',6)
    title(['t=',num2str(fc_time_uni(i))])
end

msum = reshape(sum(m_ieg>1),length(T_cells_tmp_uni),5);
figure('position',[200,200,1200,1000],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.05,0.05], [0.02,0.02], [0.05,0.05]);
for i=1:length(fc_time_uni)
    s = msum(:,i);
    [ssort,sxi] = sort(s,'descend');
    axes(ha(i));
    scatter(ssort,[1:length(T_cells_tmp_uni)],clustersize(sxi));
    set(ha(i),'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',regexprep(clusteruni(sxi),'_','-'),'fontsize',6)
    title(['t=',num2str(fc_time_uni(i))])
end


figure;
set(gcf,'color','w','position',[20,20,1200,1000]);
[ha, pos] = tight_subplot(8, ceil(length(unique(T_cells_tmp(:,1)))/8), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    bar(msum(k,:));
    title(regexprep(clusteruni(k),'_','-'))
end
    

z = linkage(m,'ward','correlation');
d = corr_mat(m');
leaforder = optimalleaforder(z,squareform(1-d,'tovector'));
tp_order = reshape([1:4*length(T_cells_tmp_uni)],4,length(T_cells_tmp_uni))';
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
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',g(leaforder),'fontsize',10,'xtick'...
    ,[length(T_cells_tmp_uni)/2:length(T_cells_tmp_uni):4*length(T_cells_tmp_uni)],'XTickLabel',{'0','2','8','24','recall'})
if savefig_flag==100
    %     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
    eval(['export_fig top100_mostdiff_FC_Glut1_',date,'.pdf']);    
end


inrmv = (sum(m_ieg,2)==0);
g = IEG_markers(~inrmv);
m_ieg(inrmv,:) = [];
z = linkage(m_ieg,'ward','correlation');
d = corr_mat(m_ieg');
leaforder = optimalleaforder(z,squareform(1-d,'tovector'));
tp_order = reshape([1:5*length(T_cells_tmp_uni)],5,length(T_cells_tmp_uni))';
tp_order = tp_order(:);
figure('position',[200,200,1200,1000],'color','w');
set(gca,'Position',[0.1,0.03,0.85,0.9])
imagesc((m_ieg(leaforder,tp_order)),[0,3]); hold on;
colormap('summer');
freezeColors(gca);
yl = get(gca,'ylim');
for i=1:length(fc_time_uni)-1
    plot(i*length(T_cells_tmp_uni)*[1,1]+0.5,yl,'-k');
end
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',g(leaforder),'fontsize',10,'xtick'...
    ,[length(T_cells_tmp_uni)/2:length(T_cells_tmp_uni):5*length(T_cells_tmp_uni)],'XTickLabel',{'0','2','8','24','recall'})



%%
% [ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp(:,1),data_sorted_all,1);
list =  {'Lynx1','Diras2','Slit2','Cdh8'....
    ,'Neurod6','Bmp3','Rspo2','Man1a','Sim1','C1ql1','Tac1','Meis2','Cartpt','Prox1'};



% list = {'Fth1','Ftl1','Scrib','Tfr2','Tfrc','Slc11a2','Slc40a1','Aco1','Ireb2'};
figure;
set(gcf,'color','w','position',[20,20,1000,1000]);
[ha, pos] = tight_subplot(length(list), 1, [0.01,0.01], [0.02,0.01], [0.15,0.01]);
logflag = 0;
for jjj=1:length(list)
    gn = list{jjj};
    g = find(strcmpi(gn,geneid));
    p99 = prctile(data(g,:),100);
    axes(ha(jjj))
    p90 = max(ceil(prctile(data(g,:),90)),1);
    t_ed = zeros(length(T_cells_tmp_uni),1);
    t_av = zeros(length(T_cells_tmp_uni),1);
    t_75 = zeros(length(T_cells_tmp_uni),1);
    for k=1:length(T_cells_tmp_uni)
        c=k;
        gr2 = find(T_cells_tmp(:,1)==c);
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
            plot(k + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.'); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
        end
        t_ed(k) = median(data(g,gr2));
        t_av(k) = mean(data(g,gr2));
        t_75(k) = prctile(data(g,gr2),75);
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
if savefig_flag==100
    %     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
    eval(['export_fig markers1_violin_final_step3_FC_Glut1_',date,'.pdf']);
    
end


%%
% [ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp(:,1),data_sorted_all,1);
% list =  {'Maf','Nfib','Zeb2','Htr3a','Sncg','Vip','Chat','Reln','Ndnf','Lamp5','Pvalb','Sst','Chodl'};%geneid(ind_gr_tmp_mark(1:end));
list =  {'Mpped2','Krt12','Gpr88','Trhr','Cd24a','Igfbpl1','Pde11a','Trh','Rab3b','Kctd1','Gsg1l'};
% list = {'Meis2','Maf','Zeb2','Prlr','Lhx6','Npas1','Igfbp4','Sst','Cryab','Penk','Moxd1','Prkcd','Unc5b','Ucn3','Pthlh','Vip','Chat','Htr3a','Vdr'};
% list = {'Fth1','Ftl1','Scrib','Tfr2','Tfrc','Slc11a2','Slc40a1','Aco1','Ireb2'};
figure;
set(gcf,'color','w','position',[20,20,1000,1000]);
[ha, pos] = tight_subplot(length(list), 1, [0.01,0.01], [0.02,0.01], [0.15,0.01]);
logflag = 0;
for jjj=1:length(list)
    gn = list{jjj};
    g = find(strcmpi(gn,geneid));
    p99 = prctile(data(g,:),100);
    axes(ha(jjj))
    p90 = max(ceil(prctile(data(g,:),90)),1);
    t_ed = zeros(length(T_cells_tmp_uni),1);
    t_av = zeros(length(T_cells_tmp_uni),1);
    t_75 = zeros(length(T_cells_tmp_uni),1);
    for k=1:length(T_cells_tmp_uni)
        c=k;
        gr2 = find(T_cells_tmp(:,1)==c);
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
            plot(k + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.'); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
        end
        t_ed(k) = median(data(g,gr2));
        t_av(k) = mean(data(g,gr2));
        t_75(k) = prctile(data(g,gr2),75);
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
if savefig_flag==100
%     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
    eval(['export_fig markers2_violin_final_step3_FC_Glut1_',date,'.pdf']);
end


%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
top_g = 50;
gr1name = '0';
gr2name = '2';
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d'};%,'Snca','Gm9843','Gm12027','Gm12031'
gr1 = find((T_cells_tmp==2 & fc_time==0 & batch_flag==2));
gr2 = find((T_cells_tmp==2 & fc_time==0 & batch_flag==3));
% in = find(mean(data_orig_all_sorted(:,gr1)>0,2)>0.03 | mean(data_orig_all_sorted(:,gr2)>0,2)<0.7);
in = find(mean(data>0,2)>0.05 & mean(data>0,2)<0.5 & ~ismember(geneid(:,1),sex_genes));
ptt = zeros(length(in),1);
for s=1:length(in)
    if mod(s,100)==0
        s
    end
    ptt(s) = ranksum(data(in(s),gr1),data(in(s),gr2),'tail','both');
end
ptt(isnan(ptt)) = 1;
ptt(ptt<1) = qval_from_pval(ptt(ptt<1));
ptt(ptt<1e-300) = 1e-300;
x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
d = x1-x2 ;
figure('position',[200,200,1400,580],'color','w');
[ha, pos] = tight_subplot(1, 4, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
axes(ha(1))
plot(d(in),-log10(ptt),'.'); hold on;
[~,xi] = sort(ptt);
plot(d(in(xi(1:100))),-log10(ptt(xi(1:100))),'or');
text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid(in(xi(1:100))))
xlabel([gr1name,'-',gr2name]);
ylabel(['-log10(q)']);
grid on;
title(['comparing ',gr1name,' to ',gr2name])
[~,xi] = sort(d);
axes(ha(2))
plot(x1, x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,xmax],'-k'); grid on
plot([0,xmax],[0,xmax]+1,'--k'); grid on
plot([1,xmax],[0,xmax-1],'--k'); grid on
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',8);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',8);
xlabel(gr1name)
ylabel(gr2name)
axis tight

top_g = 50;
x1 = mean(data(:,gr1)>0,2);
x2 = mean(data(:,gr2)>0,2);
d = x1-x2 ;
[~,xi] = sort(d);
axes(ha(3))
plot(x1, x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,xmax],'-k'); grid on
plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',8);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',8);
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight

top_g = 50;
x1 = log2(prctile(data(:,gr1),90,2)+1);
x2 = log2(prctile(data(:,gr2),90,2)+1);
d = x1-x2 ;
d(max([x1,x2],[],2)<2) = 0;
[~,xi] = sort(d);
axes(ha(4))
plot(x1, x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,xmax],'-k'); grid on
plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',8);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',8);
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight


figure;
top_g = 50;
x1 = log2(prctile(data(:,gr1),90,2)+1);
x2 = log2(prctile(data(:,gr2),90,2)+1);
d = x1-x2 ;
d(max([x1,x2],[],2)<2) = 0;
[~,xi] = sort(d);

plot((x1 + x2)/2, x1-x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,0],'-k'); grid on
plot([0,xmax],[1,1]+0.4,'--k'); grid on
plot([0,xmax],[1,1]*0.4,'--k'); grid on
plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,'.r'); hold on;
text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
xi = flipud(xi);
plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,'.r'); hold on;
text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight
% eval(['export_fig scatter_volcano_',gr1name,'_',gr2name,'_',clusteruni((j)),'_',date,'.pdf']);

%%
gr1 = find((T_cells_tmp==2 & (fc_time==0) & batch_flag==3));
gr2 = find((T_cells_tmp==2 & fc_time==2 & batch_flag==3));
x1 = log2(prctile(data(:,gr1),90,2)+1);
x2 = log2(prctile(data(:,gr2),90,2)+1);
p90 = prctile(data(:,[gr1;gr2]),90,2);
f90 = sum(data(:,gr2)>repmat(p90,1,length(gr2)),2);
pp90 = binocdf(f90,length(gr2),0.1,'upper');
pp90(max([x1,x2],[],2)<4) = 1;
[~,xi] = sort(pp90);
geneid(xi(1:100))

%%
for c=1:length(clusteruni)
    c
    % c = 3;
    list = {'Arc','Fos','Fosb','Fosl2','Npas4','Btg2','Bdnf','Junb','Egr1','Egr2','Egr4','Nr4a1','Nr4a3','Dusp1'};
    hf = figure('color','w','position',[20,20,1200,200],'Name',clusteruni{c});
    [ha, pos] = tight_subplot(1, length(list), [0.02,0.02], [0.05,0.12], [0.03,0.03]);
    for jjj=1:length(list)
        
        gn = list{jjj};
        g = find(strcmpi(gn,geneid));
        
        if min(grsz_batch23(c,1))>30
            axes(ha(jjj))
            b23_50 = p50_batch23_cluster_tp{c}(g,:)+1;
            b23_50 = b23_50/b23_50(1);
            b23_90 = p90_batch23_cluster_tp{c}(g,:)+1;
            b23_90 = b23_90/b23_90(1);
            
            plot([1,2,3,4,5],b23_50,'--sb'); hold on;
            plot([1,2,3,4,5],b23_90,'-sc'); hold on;
%             text([1,2,3,4,5],7.5*[1,1,1,1,1],cellfun(@num2str, m2c(grsz_batch23(c,:)),'UniformOutput',0),'fontsize',6)
            set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear','YScale',yscale,'xlim',[1,5],'ylim',[0,4])
            
            %[gn,',c=',num2str(c)])
            if jjj==1
                ylabel(clusteruni{c},'Interpreter','none','fontsize',7)
                title(gn);
%                 title([clusteruni{c},',',gn],'Interpreter','none');
            else
                title(gn);
            end
            if jjj==length(list)
                text([1,2,3,4,5],4-[linspace(0,1,5)],cellfun(@num2str, m2c(grsz_batch23(c,:)),'UniformOutput',0),'fontsize',6)
            end
        end
    end
    linkaxes(ha,'xy');
    eval(['export_fig Glut1_IEG_FC_220411/',clusteruni{c},'_IEG_FC',date,'.pdf']);
%     pause
    close(hf)
end
%%
%%
list = {'Arc','Fos','Fosb','Fosl2','Npas4','Btg2','Bdnf','Junb','Egr1','Egr2','Egr4','Nr4a1','Nr4a3','Dusp1'};
[~,loc] = ismember(list,geneid);
maxind_ieg = zeros(length(list),length(T_cells_tmp_uni));
for i=1:length(T_cells_tmp_uni)
    if min(grsz_batch23(i,1))>30
        tmp = p90_batch23_cluster_tp{i}(loc,:)+1;
        tmp = tmp./repmat(tmp(:,1),1,5);
        maxind_ieg(:,i) = max(tmp,[],2);
    end
end
hf = figure('color','w','position',[20,20,800,1000],'Name',clusteruni{c});
[ha, pos] = tight_subplot(1, length(list), [0.01,0.002], [0.05,0.02], [0.12,0.03]);
for i=1:length(list)
    axes(ha(i))
    barh(maxind_ieg(i,:),'BarWidth',0.3);
    set(ha(i),'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',cell(1,length(T_cells_tmp_uni))....
        ,'xlim',[1,4],'TickLabelInterpreter','none','fontsize',6)
    title(list{i})
end
set(ha(1),'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',clusteruni,'TickLabelInterpreter','none','fontsize',6)

hf = figure('color','w','position',[20,20,800,1000],'Name',clusteruni{c});
axes('position',[0.15,0.1,0.72,0.85])
imagesc(maxind_ieg'>2)
colormap('gray')
set(gca,'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',clusteruni,'TickLabelInterpreter','none','fontsize',6,...
    'xtick',[1:length(list)],'XTickLabel',list,'ylim',[0.5,length(T_cells_tmp_uni)+0.5])
axes('position',[0.9,0.1,0.07,0.85])
barh(sum(maxind_ieg'>2,2))
axis tight
set(gca,'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',cell(1,length(T_cells_tmp_uni))....
    ,'fontsize',6,'ydir','reverse','ylim',[0.5,length(T_cells_tmp_uni)+0.5])


%%

ieggenes = {'Egr1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'};
[~,loc] = ismember(ieggenes,geneid);
pcth = 95;
d = data(loc,:);
d95 = d>repmat(prctile(d,95,2),1,length(d));
iegfram_mat = cell(length(T_cells_tmp_uni),1);
for k=1:length(T_cells_tmp_uni)
    k
    tmp =zeros(length(ieggenes)+1,5);    
    for jj=1:5
        gr1 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag>=2));
        tmp(:,jj) = [mean(d95(:,gr1),2);mean(sum(d95(:,gr1))>0)];
    end
    iegfram_mat{k} = tmp;
end


c = cell2mat(iegfram_mat);
n = 56*5;
% c = c(:,[1:5:n,2:5:n,3:5:n,4:5:n,5:5:n]);
figure('position',[100,100,600,900],'color','w'); 
imagesc(cent_norm(c)); hold on;
for i=1:length(T_cells_tmp_uni)
    plot([0.5,5.5],i*(length(ieggenes)+1)*[1,1]+0.5,'k');
    text([1:5],i*length(ieggenes)*[1,1,1,1,1]-2, cellfun(@num2str, m2c(grsz_batch23(i,:)),'UniformOutput',0),'fontsize',7);
end
text(4*ones(length(T_cells_tmp_uni),1),[4:(length(ieggenes)+1):(length(ieggenes)+1)*length(T_cells_tmp_uni)]...
    ,regexprep(clusteruni,'_','-'),'fontsize',7);
set(gca,'ytick',[1:length(c)],'yticklabel',repmat([ieggenes,{'any'}]',length(T_cells_tmp_uni),1));   

tmp = [];
for i=1:length(T_cells_tmp_uni)
    tmp = [tmp;[sum(cent_norm(iegfram_mat{i}(1:7,:)))]];
end
s = [(max(tmp(:,2:5),[],2)-tmp(:,1))];
[~,xi] = sort(s,'descend');
sgrsz = sum(grsz_batch23,2);
xi = xi(sgrsz(xi)>120);
s = s(xi);
sgrsz = sgrsz(xi);
figure('position',[100,100,600,900],'color','w');
axes('Position',[0.1,0.1,0.6,0.85])
scatter([1:length(s)],s,sgrsz*0.3,'filled');
text([1:length(s)]+1,s,clusteruni(xi),'Interpreter','none','fontsize',6)
axis tight
ylabel('IEG score');
xlabel('rank');
eval(['export_fig Glut1_IEGfrac_score_scatter',date,'.pdf']);

c = cell2mat(iegfram_mat(xi));
figure('position',[100,100,600,900],'color','w'); 
imagesc(cent_norm(c)); hold on;
for i=1:length(xi)
    plot([0.5,5.5],i*(length(ieggenes)+1)*[1,1]+0.5,'k');
    text([1:5],i*(length(ieggenes)+1)*[1,1,1,1,1]-2, cellfun(@num2str, m2c(grsz_batch23(xi(i),:)),'UniformOutput',0),'fontsize',7);
end
text(4*ones(length(xi),1),[4:(length(ieggenes)+1):(length(ieggenes)+1)*length(xi)]...
    ,regexprep(clusteruni(xi),'_','-'),'fontsize',7);
set(gca,'ytick',[1:length(c)],'yticklabel',repmat([ieggenes,{'any'}]',length(T_cells_tmp_uni),1));   
set(gca,'xtick',[1:5],'xticklabel',{'0','2','8','24','28'})


figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    axes(ha(i))
    for j=[2,3,4,6];%length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    title(regexprep(clusteruni{xi(i)},'_','-'))
end
legend(ieggenes([2,3,4,6]))
eval(['export_fig top5_clusters_Glut1_IEGfrac_',date,'.pdf']);

