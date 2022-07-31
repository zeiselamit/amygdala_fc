tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;

load afterloading_QC3000_FC_31-Mar-2022
IEG_markers = loadCellFile('IEG_cluster.txt');
[~,in_ieg] = ismember(IEG_markers,geneid);

gaba_clustername = loadCellFile('ClusterOrder_step3_Glut1_FC_24-Oct-2021.txt');
cluster_info = loadCellFile('Clusters_info_hannah.txt');
cluster_info = cluster_info(58:end,:);

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

fc_time2 = fc_time;
fc_time2(fc_time==0 & batch_flag==2)= 1;
fc_time_uni2 = unique(fc_time2);
% % % % 
clusteruni = unique(cellid_clusters(:,5));
tmp = cellfun(@(x) x(7:end),clusteruni,'UniformOutput',0);
tmp = regexprep(tmp,'Cd63','Cd36');
[~,loc] = ismember(cluster_info(:,5),tmp);
clusteruni = clusteruni(loc);
clusteruni = [clusteruni,cellfun(@(x,y) [x,'-',y] , cluster_info(:,4),cluster_info(:,5),'UniformOutput',0)];
cluster_location = cluster_info(:,6);
cluster_loc_flag = cell(size(cellid));
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,5),clusteruni{i}));
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
gn = 'snca';
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
fc_time2 = fc_time;
fc_time2(fc_time==0 & batch_flag==2)= 1;
fc_time_uni2 = unique(fc_time2);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    t_90 =[];
    c=k;
    t_ed = zeros(5,1);
    t_av = zeros(5,1);
    t_75 = zeros(5,1);
    p90 = prctile(data(g,T_cells_tmp(:,b)==c & batch_flag>=2),90);
    for i=[1:6]
        gr2 = find(T_cells_tmp(:,b)==c & fc_time2==fc_time_uni2(i) & batch_flag>=2);%find(fc_time_sorted==fc_time_uni(i));%
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
% ieggenes = {'Bdnf','Fos','Fosl2','Fosb','Jun','Junb','Egr1','Egr2','Egr4'...
%     ,'Dusp1','Dusp5','Dusp6','Arc','Btg1','Btg2','Npas4','Nr4a1','Homer1','Ptgs2','Atf3','Cdkn1a'};
ieggenes = {'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'};
colorvec = distinguishable_colors(length(ieggenes));
[~,loc] = ismember(ieggenes,geneid);
pcth = 90;
d = data(loc,:);
d95 = d>repmat(prctile(d,pcth,2),1,length(d));
iegfram_mat = cell(length(T_cells_tmp_uni),1);
iegfram_mat_fold = cell(length(T_cells_tmp_uni),1);
for k=1:length(T_cells_tmp_uni)
    k
    if min(grsz_batch23(k,:),[],2)>0
        tmp =zeros(length(ieggenes)+1,5);
        for jj=1:5
            gr1 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag>=2));
            tmp(:,jj) = [mean(d95(:,gr1),2);mean(sum(d95(:,gr1))>0)];
        end
        iegfram_mat{k} = tmp;
        tmp1 = tmp(:,1);
        tmp1(tmp1==0) = 1/grsz_batch23(k,1);
        iegfram_mat_fold{k} = tmp./repmat(tmp1,1,length(fc_time_uni));
    end
end


% c = cell2mat(iegfram_mat);
% n = 56*5;
% % c = c(:,[1:5:n,2:5:n,3:5:n,4:5:n,5:5:n]);
% figure('position',[100,100,600,900],'color','w'); 
% imagesc(cent_norm(c)); hold on;
% for i=1:length(T_cells_tmp_uni)
%     plot([0.5,5.5],i*(length(ieggenes)+1)*[1,1]+0.5,'k');
%     text([1:5],i*length(ieggenes)*[1,1,1,1,1]-2, cellfun(@num2str, m2c(grsz_batch23(i,:)),'UniformOutput',0),'fontsize',7);
% end
% text(4*ones(length(T_cells_tmp_uni),1),[4:(length(ieggenes)+1):(length(ieggenes)+1)*length(T_cells_tmp_uni)]...
%     ,regexprep(clusteruni,'_','-'),'fontsize',7);
% set(gca,'ytick',[1:length(c)],'yticklabel',repmat([ieggenes,{'any'}]',length(T_cells_tmp_uni),1));   

tmp = zeros(size(grsz_batch23));
for i=1:length(T_cells_tmp_uni)
    if min(grsz_batch23(i,:),[],2)>0
        tmp(i,:) = [sum((iegfram_mat{i}(1:length(ieggenes),:)))];
    end
end
% tmp2 = zeros(size(grsz_batch23));
% for i=1:length(T_cells_tmp_uni)
%     if min(grsz_batch23(i,:),[],2)>20
%         tmp2(i,:) = [sum(iegfram_mat_fold{i}(1:length(ieggenes),:))];
%     end
% end
v = (min(grsz_batch23,[],2))>20;
s = [(max(tmp(:,2:5),[],2)-tmp(:,1))];
% s = max(tmp2,[],2);
[~,xi] = sort(s,'descend');
sgrsz = sum(grsz_batch23,2);
min_sgrsz = min(grsz_batch23,[],2);
% xi = xi(s(xi)~=0);
s = s(xi);
v = v(xi);
sgrsz = sgrsz(xi);
min_sgrsz = min_sgrsz(xi);
figure('position',[100,100,400,800],'color','w');
axes('Position',[0.1,0.1,0.6,0.85])
% scatter([1:length(s(v))],s(v),sgrsz(v)*0.3,repmat([1,0,0],sum(v),1),'filled');hold on;
% scatter([1:length(s(v==0))],s(v==0),sgrsz(v==0)*0.3,repmat([1,1,1]*.5,sum(v==0),1),'filled');
markergene = min_sgrsz;
markergene(markergene<10) = 10;
markergene(markergene>60) = 60;
c_rgb = [0,0.4,0.2];
valid = sgrsz>60;
marker_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene)...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene)];
scatter([1:length(s(valid))],s(valid),sgrsz(valid)*0.3,marker_color(valid,:),'filled');hold on;
text([1:length(s(valid))]+2,s(valid),clusteruni(xi(valid),2),'Interpreter','none','fontsize',6)

axis tight
ylabel('IEG score');
xlabel('rank');
set(gca,'xlim',[-5,length(xi(valid))+8])
eval(['export_fig Glut12_IEGfrac_score_scatter_p90_grsz20',date,'.pdf']);

table1 = [clusteruni(xi,2),m2c([s,sgrsz,min_sgrsz])];
saveCellFile(table1,'Glut12_IEGfrac_score_p90_grsz20.txt')

c = cell2mat(iegfram_mat(xi));
c_fold = cell2mat(iegfram_mat_fold(xi));
% figure('position',[100,100,600,900],'color','w'); 
% imagesc(cent_norm(c)); hold on;
% for i=1:length(xi)
%     plot([0.5,5.5],i*(length(ieggenes)+1)*[1,1]+0.5,'k');
%     text([1:5],i*(length(ieggenes)+1)*[1,1,1,1,1]-2, cellfun(@num2str, m2c(grsz_batch23(xi(i),:)),'UniformOutput',0),'fontsize',7);
% end
% text(4*ones(length(xi),1),[4:(length(ieggenes)+1):(length(ieggenes)+1)*length(xi)]...
%     ,regexprep(clusteruni(xi),'_','-'),'fontsize',7);
% set(gca,'ytick',[1:length(c)],'yticklabel',repmat([ieggenes,{'any'}]',length(T_cells_tmp_uni),1));   
% set(gca,'xtick',[1:5],'xticklabel',{'0','2','8','24','28'})

glut12_cluster_by_iegscore = clusteruni(xi(v));
figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
intop = xi(v);
for i=1:5
    axes(ha(i))
    for j=1:length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat{intop(i)}(j,:),'-s','color',colorvec(j,:)); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    set(ha(i),'xtick',[0,2,8,24,26])
    title(regexprep(clusteruni{intop(i),2},'_','-'))
end
legend(ieggenes)
eval(['export_fig top5_clusters_Glut12_IEGfrac_p90_grsz20_',date,'.pdf']);

figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    axes(ha(i))
    for j=1:length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat_fold{intop(i)}(j,:),'-s','color',colorvec(j,:)); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    set(ha(i),'xtick',[0,2,8,24,26])
    title(regexprep(clusteruni{intop(i),2},'_','-'))
end
legend(ieggenes)
eval(['export_fig top5_clusters_Glut12_IEG_fold_frac__p90_grsz20_',date,'.pdf']);


figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
intop = flipud(xi(v));
for i=1:5
    axes(ha(i))
    for j=1:length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat{intop(i)}(j,:),'-s','color',colorvec(j,:)); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    set(ha(i),'xtick',[0,2,8,24,26])
    title(regexprep(clusteruni{intop(i),2},'_','-'))
end
legend(ieggenes)
eval(['export_fig low5_clusters_Glut12_IEGfrac_p90_grsz20_',date,'.pdf']);

figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    axes(ha(i))
    for j=1:length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat_fold{intop(i)}(j,:),'-s','color',colorvec(j,:)); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    set(ha(i),'xtick',[0,2,8,24,26])
    title(regexprep(clusteruni{intop(i),2},'_','-'))
end
legend(ieggenes)
eval(['export_fig low5_clusters_Glut12_IEG_fold_frac_p90_grsz20_',date,'.pdf']);
% figure('position',[100,100,1400,1000],'color','w');
% [ha, pos] = tight_subplot(6, ceil(length(xi)/6), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
% m = [];
% for i=1:length(xi)
%     axes(ha(i))
%     for j=1:length(ieggenes)
%         plot([0,2,8,24,26],iegfram_mat{xi(i)}(j,:),'-s','color',colorvec(j,:)); hold on;
% %         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
%     end
%     set(ha(i),'xtick',[0,2,8,24,26])
%     title(regexprep(clusteruni{xi(i)},'_','-'))
%     axis tight
%     m = max([m,max(get(ha(i),'ylim'))]);
% end
% linkaxes(ha,'y')
% legend(ieggenes)
% set(ha(1),'ylim',[0,m])
% eval(['export_fig top5_clusters_Glut12_IEGfrac_',date,'.pdf']);


figure('position',[100,100,1400,1400],'color','w');
[ha, pos] = tight_subplot(2, 4, [0.05,0.13], [0.04,0.04], [0.1,0.02]);
tp_ptp_gene = zeros(length(ieggenes),length(fc_time_uni));
for i=1:length(ieggenes)
    crowt = repmat([ieggenes,{'any'}]',length(xi),1);
    c_fos = c_fold(strcmpi(crowt,ieggenes{i}),:);
    [~,imax] = max(c_fos,[],2);
    [~,fosxi] = sort(imax);
    tp_ptp_gene(i,:) = [sum(imax==1),sum(imax==2),sum(imax==3),sum(imax==4),sum(imax==5)];
%     fosxi = 1:length(xi);
    % figure('position',[200,200,280,480],'color','w');
    axes(ha(i))
    imagesc((c_fos(fosxi,:)),[1,5])
    set(gca,'ytick',[1:length(xi)],'yticklabel',clusteruni(xi(fosxi),2),'ticklabelInterpreter','none','fontsize',7);
    set(gca,'xtick',[1:5],'xticklabel',{'0','2','8','24','28'})
    title(ieggenes{i})
    colormap('summer')
end
eval(['export_fig IEG_pergene_heatmap_Glut12_p90_grsz20_',date,'.pdf']);
%% IEG fraction per time point
ieggenes = {'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'};
[~,loc] = ismember(ieggenes,geneid);
pcth = 90;
d = data(loc,:);
d95 = d>repmat(prctile(d,pcth,2),1,length(d));
iegpos = sum(d95>0)>2;

t = xi(v);
top5types = false(size(fc_time));
for i=1:5
    top5types(T_cells_tmp==t(i)) = 1;
end
t = flipud(t);
low5types = false(size(fc_time));
for i=1:5
    low5types(T_cells_tmp==t(i)) = 1;
end

for i=1:length(fc_time_uni)
    frac_t(i) = sum(fc_time_uni(i)==fc_time & top5types);
    frac_t_ieg(i) = sum(fc_time_uni(i)==fc_time & iegpos' & top5types);
end
sum_t = round(frac_t*mean(iegpos(top5types)));
frac_t = frac_t/sum(frac_t);
sum_t_ieg = frac_t_ieg;
frac_t_ieg  = frac_t_ieg/sum(frac_t_ieg);
top5ratio = frac_t_ieg./frac_t

figure('Position',[100,100,1400,300],'Color','w');
subplot(1,4,1);
sb = bar([frac_t;frac_t_ieg],'stacked'); hold on;
labelStackedBar(sb,strsplit(num2str([sum_t,sum_t_ieg])))
set(gca,'xtick',[1,2],'xticklabel',{'all','IEG+'})
title('top5clusters')

for i=1:length(fc_time_uni)
    frac_t(i) = sum(fc_time_uni(i)==fc_time & low5types);
    frac_t_ieg(i) = sum(fc_time_uni(i)==fc_time & iegpos' & low5types);
end
sum_t = round(frac_t*mean(iegpos(low5types)));
frac_t = frac_t/sum(frac_t);
sum_t_ieg = frac_t_ieg;
frac_t_ieg  = frac_t_ieg/sum(frac_t_ieg);
low5ratio = frac_t_ieg./frac_t

subplot(1,4,2);
sb = bar([frac_t;frac_t_ieg],'stacked'); hold on;
labelStackedBar(sb,strsplit(num2str([sum_t,sum_t_ieg])))
set(gca,'xtick',[1,2],'xticklabel',{'all','IEG+'})
title('low5clusters')

for i=1:length(fc_time_uni)
    frac_t(i) = sum(fc_time_uni(i)==fc_time);
    frac_t_ieg(i) = sum(fc_time_uni(i)==fc_time & iegpos');
end
sum_t = round(frac_t*mean(iegpos));
frac_t = frac_t/sum(frac_t);
sum_t_ieg = frac_t_ieg;
frac_t_ieg  = frac_t_ieg/sum(frac_t_ieg);
allratio = frac_t_ieg./frac_t

subplot(1,4,3);
sb = bar([frac_t;frac_t_ieg],'stacked'); hold on;
labelStackedBar(sb,strsplit(num2str([sum_t,sum_t_ieg])))
set(gca,'xtick',[1,2],'xticklabel',{'all','IEG+'})
title('Allclusters')

subplot(1,4,4);
plot(allratio,'s-'); hold on;
plot(top5ratio,'o-')
plot(low5ratio,'d-')
legend({'all','top10','low10'})
set(gca,'xtick',[1,2,3,4,5],'xticklabel',{'0','2','8','24','28'})
title(['ratio IEG+/all, prct=',num2str(pcth),',IEG+=',num2str(round(100*mean(iegpos))),'%'])
eval(['export_fig ratio_IEG_pertimepoint_Glut12_p90_3IEG_grsz20_top5_',date,'.pdf']);
%% IEG +/- scatters
ieggenes = {'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'};
[~,loc] = ismember(ieggenes,geneid);
pcth = 90;
d = data(loc,:);
d95 = d>repmat(prctile(d,pcth,2),1,length(d));
iegpos = sum(d95>0,1)>0;

sum(iegpos)
indcorrg = [];
indcorrg_d = [];
for k=1:length(T_cells_tmp_uni)
    gr1 = find((T_cells_tmp==k & fc_time>=2 & batch_flag>=2 & ~iegpos'));
    gr2 = find((T_cells_tmp==k & fc_time>=2 & batch_flag>=2 & iegpos'));
    if min([length(gr1),length(gr2)])>30
        x1 = mean(log2(data(:,gr1)+1),2);
        x2 = mean(log2(data(:,gr2)+1),2);
        d = x2-x1 ;
        d(max([x1,x2],[],2)<1) = 0;
        [~,xi] = sort(d,'descend');
        indcorrg = [indcorrg;xi(1:50)];
        [~,xi] = sort(d,'ascend');
        indcorrg_d = [indcorrg_d;xi(1:50)];
    end
end
indcorrg = unique(indcorrg);
indcorrg_d = unique(indcorrg_d);

% r = corr(iegpos',data');
% r(isnan(r)) = 0;
% [~,xi] = sort(r,'descend');
% indcorrg = xi(1:500);
r_top_pertype = zeros(length(indcorrg),length(T_cells_tmp_uni));
d_top_pertype = zeros(length(indcorrg),length(T_cells_tmp_uni));
p_top_pertype = zeros(length(indcorrg),length(T_cells_tmp_uni));

d_down_pertype = zeros(length(indcorrg_d),length(T_cells_tmp_uni));
p_down_pertype = zeros(length(indcorrg_d),length(T_cells_tmp_uni));

top_g = 150;
gr1name = 'gr1-IEG(-)';
gr2name = 'gr2-IEG(+)';
for k=1:length(T_cells_tmp_uni)
    k
    %     if min(grsz_batch23(k,1))>30
    gr1 = find((T_cells_tmp==k & fc_time>=2 & batch_flag>=2));
    if gr1>100
        tmp = corr(iegpos(gr1)',data(indcorrg,gr1)');
        r_top_pertype(:,k) = tmp';
    end
    sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d','Erdr1','Mid1',};%,'Snca','Gm9843','Gm12027','Gm12031'
    gr1 = find((T_cells_tmp==k & fc_time>=2 & batch_flag>=2 & ~iegpos'));
    gr2 = find((T_cells_tmp==k & fc_time>=2 & batch_flag>=2 & iegpos'));
    if min([length(gr1),length(gr2)])>30
        in = find(mean(data>0,2)>0.05 & mean(data>0,2)<0.5 & ~ismember(geneid(:,1),sex_genes));
        ptt = zeros(length(in),1);
        for s=1:length(in)
            %             if mod(s,100)==0
            %                 s
            %             end
                ptt(s) = ranksum(data(in(s),gr1),data(in(s),gr2),'tail','both');
%             [~,ptt(s)] = kstest2(data(in(s),gr1),data(in(s),gr2));
        end
        ptt(isnan(ptt)) = 1;
        ptt(ptt<1) = qval_from_pval(ptt(ptt<1));
        ptt(ptt<1e-300) = 1e-300;
        ptmp = ones(length(geneid),1);
        ptmp(in) = ptt;
        p_top_pertype(:,k) = -log10(ptmp(indcorrg));
        p_down_pertype(:,k) = -log10(ptmp(indcorrg_d));
        x1 = mean(log2(data(:,gr1)+1),2);
        x2 = mean(log2(data(:,gr2)+1),2);
        d = x2-x1 ;
        hf = figure('position',[200,200,1400,580],'color','w','name',clusteruni{k});
        [ha, pos] = tight_subplot(1, 3, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
        axes(ha(1))
        plot(d(in),-log10(ptt),'.'); hold on;
        [~,xi] = sort(ptt);
        plot(d(in(xi(1:100))),-log10(ptt(xi(1:100))),'or');
        text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid(in(xi(1:100))))
        xlabel([gr1name,'-',gr2name]);
        ylabel(['-log10(q)']);
        grid on;
        title(['comparing ',gr1name,' to ',gr2name,', gr1=',num2str(length(gr1))....
            ,', gr2=',num2str(length(gr2)),', IEGfrac=',num2str(length(gr2)/length([gr1;gr2]))])
        [~,xi] = sort(d);
        axes(ha(2))
        
        %
        d(max([x1,x2],[],2)<1) = 0;
        [~,xi] = sort(d);
        d_top_pertype(:,k) = d(indcorrg);
        d_down_pertype(:,k) = d(indcorrg_d);
        plot((x1 + x2)/2, x2-x1, '.');hold on;
        xmax = max(x1);
        plot([0,xmax],[0,0],'-k'); grid on
        plot([0,xmax],[0,0]+1,'--k'); grid on
        plot([0,xmax],[0,0]-1,'--k'); grid on
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
        xi = flipud(xi);
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);

        xlabel('(gr1+gr2)/2')
        ylabel('gr1-gr2')
        axis tight
        yl = get(gca,'ylim');
        yl(2) = max(yl(2),1);yl(1) = min(yl(1),-1);
        set(gca,'ylim',yl);
        
        top_g = 100;
        x1 = mean(data(:,gr1)>0,2);
        x2 = mean(data(:,gr2)>0,2);
        d = x1-x2 ;
        d(max([x1,x2],[],2)<0.2) = 0;
        [~,xi] = sort(d);
        axes(ha(3))
        %
        plot((x1 + x2)/2, x2-x1, '.');hold on;
        xmax = max(x1);
        plot([0,xmax],[0,0],'-k'); grid on
        plot([0,xmax],[0,0]+0.3,'--k'); grid on
        plot([0,xmax],[0,0]-0.3,'--k'); grid on
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
        xi = flipud(xi);
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x2(xi(1:top_g)) - x1(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);        
        %
         xlabel('(gr1+gr2)/2')
        ylabel('gr1-gr2')
        axis tight
        yl = get(gca,'ylim');
        yl(2) = max(yl(2),0.5);yl(1) = min(yl(1),-0.5);
        set(gca,'ylim',yl);
%         eval(['export_fig Glut12_scatter_IEG_pos_neg_FC_220423/',clusteruni{k},'_scatter_IEG_FC',date,'.pdf']);
        pause
        close(hf)
    end
    
end
d_top_pertype(isnan(d_top_pertype)) = 0;
% o = get1Dorderlinkage(d_top_pertype,1);
[z,o] = sort(sum(d_top_pertype>0.5,2),'descend');
o = o(z>1);
z = z(z>1);
v = find(~(sum(d_top_pertype==0)==length(indcorrg)));
[~,gg] = ismember(glut12_cluster_by_iegscore,clusteruni(v,1));
v = v(gg(gg>0));
[~,imax] = max(d_top_pertype(o,v),[],2);
[~,oo] = sortrows([-z,imax]);
o = o(oo);


hf = figure('position',[200,200,1400,1200],'color','w');
ax1 = axes('Position',[0.1,0.1,0.4,0.87]);
imagesc(d_top_pertype(o,v),[-1,1])
set(gca,'ytick',[1:length(indcorrg)],'YTickLabel',geneid(indcorrg(o))....
    ,'xtick',[1:length(v)],'xTickLabel',clusteruni(v,2),'xTickLabelrotation',45,.....
    'fontsize',7,'ticklabelinterpreter','none')
title('upregulation diff')
cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);
colormap((cm))
colorbar
list = geneid(indcorrg(o));
saveCellFile(list,['Glut12_up_d05_2_IEGposneg_',date,'.txt']);
% 
d_down_pertype(isnan(d_down_pertype)) = 0;
% o = get1Dorderlinkage(d_top_pertype,1);
[z,o] = sort(sum(d_down_pertype<-0.5,2),'descend');
o = o(z>1);
z = z(z>1);
v = find(~(sum(d_down_pertype==0)==length(indcorrg_d)));
[~,gg] = ismember(glut12_cluster_by_iegscore,clusteruni(v,1));
v = v(gg(gg>0));
[~,imax] = min(d_down_pertype(o,v),[],2);
[~,oo] = sortrows([-z,imax]);
o = o(oo);
ax2 = axes('Position',[0.6,0.1,0.4,0.87]);
imagesc(d_down_pertype(o,v),[-1,1])
set(gca,'ytick',[1:length(indcorrg)],'YTickLabel',geneid(indcorrg_d(o))....
    ,'xtick',[1:length(v)],'xTickLabel',clusteruni(v,2),'xTickLabelrotation',45,.....
    'fontsize',7,'ticklabelinterpreter','none')
title('downregulation diff')
cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);
colormap((cm))
colorbar
list = geneid(indcorrg_d(o));
saveCellFile(list,['Glut12_down_d05_2_IEGposneg_',date,'.txt']);
eval(['export_fig Glut12_updown_d05_2_IEGposneg_',date,'.pdf']);
% 
% ax2 = axes('Position',[0.6,0.1,0.4,0.87]);
% % o = get1Dorderlinkage(p_top_pertype,1);
% [~,o] = sort(sum(p_top_pertype>1,2),'descend');
% % figure;
% imagesc(p_top_pertype(o(1:150),v),[0,2])
% set(gca,'ytick',[1:length(indcorrg)],'YTickLabel',geneid(indcorrg(o(1:150)))....
%     ,'xtick',[1:length(v)],'xTickLabel',clusteruni(v),'xTickLabelrotation',45,.....
%     'fontsize',7,'ticklabelinterpreter','none')
% 
% set(gca,'ytick',[1:length(indcorrg)],'YTickLabel',geneid(indcorrg(o)))
% list = geneid(indcorrg(o(1:150)));
% saveCellFile(list,['Glut12_topgenes_p150_IEGposneg_',date,'.txt']);
% eval(['export_fig Glut12_top150genes_d_p_IEGposneg_',date,'.pdf']);

%% 


for k=1:length(T_cells_tmp_uni)
    k
    if min(grsz_batch23(k,:),[],2)>20
        ieggenes = {'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'};
        [~,loc] = ismember(ieggenes,geneid);
        pcth = 95;
        d = data(loc,:);
        d95 = d>repmat(prctile(d,pcth,2),1,length(d));
        iegpos = sum(d95>0)>0;
        
        
        gr1 = find((T_cells_tmp==k &  batch_flag>=2));
        t = fc_time(gr1);
        b = batch_flag(gr1);
        iegpos = iegpos(gr1);
        d = data(:,gr1);
        in = find(mean(d>0,2)>0.05 & mean(d>0,2)<0.5);
        d = d(in,:);
        
        data_tsne = cent_norm([(log2(d+1))]);
        initial_dims = 20;%length(corr_filt);
        [prj,m,D,V,Q] = pca_wis(data_tsne',initial_dims);
        D = diag(D);
        % % %
        % s = prj;
        % batchid = b;
        % usepylib = 1;
        % [sout]=harmonypy(s,batchid,usepylib);
        % prj = sout;
        % % %
        perplexity = 50;
        options = statset('MaxIter',1000);
        init = prj(:,1:2)/std(prj(:,1))*1e-4;
        init = init-repmat(mean(init),length(init),1);
        mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
            'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',5);
        
        % figure;
        % plot(mapped_xy(:,1),mapped_xy(:,2),'.')
        
        figure('color','w','name',clusteruni{k},'position',[400,10,800,1400]);
        subplot(3,1,1)
        for i=1:length(fc_time_uni)
            ind = find(t==fc_time_uni(i));
            plot(mapped_xy(ind,1),mapped_xy(ind,2),'.'); hold on;
        end
        legend('0','2','8','24','26')
        subplot(3,1,2)
        for i=2:3
            ind = find(b==i);
            plot(mapped_xy(ind,1),mapped_xy(ind,2),'.'); hold on;
        end
        legend('b2','b3')
        subplot(3,1,3)
        for i=0:1
            ind = find(iegpos==i);
            plot(mapped_xy(ind,1),mapped_xy(ind,2),'.'); hold on;
        end
        legend('-','+')
        pause
        close all
    end
end
%% 2hr vs recall
diffmat_2_recall = zeros(length(geneid),length(T_cells_tmp_uni));
for k=1:length(T_cells_tmp_uni)
    k
    gr1 = find((T_cells_tmp==k & fc_time==2 & batch_flag>=3));
    gr2 = find((T_cells_tmp==k & fc_time==28 & batch_flag>=3));
    if min([length(gr1),length(gr2)])>10

        hf = figure('Position',[100,100,800,800],'Color','w','Name',clusteruni{k,2});
        top_g = 50;
        x1 = log2(prctile(data(:,gr1),90,2)+1);
        x2 = log2(prctile(data(:,gr2),90,2)+1);

        d = x1-x2 ;
        diffmat_2_recall(:,k) = d;
        s = (x1+x2)/2 ;
        d(max([x1,x2],[],2)<3) = 0;
        [~,xi] = sort(d);

        plot((x1 + x2)/2, x1-x2, '.');hold on;
        xmax = max(x1);
        plot([0,xmax],[0,0],'-k'); grid on
        plot([0,xmax],[1,1]*0.4,'--k'); grid on
        plot([0,xmax],[1,1]*(-0.4),'--k'); grid on
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
        xi = flipud(xi);
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
        title(['gr1=',num2str(length(gr1)),',gr2=',num2str(length(gr2))])
        xlabel('(gr1+gr2)/2')
        ylabel(['2hr-recall'])
        axis tight
%         eval(['export_fig Glut12_scatter_2hr_recall_90p_220525/',clusteruni{k,2},'_scatter_2hr_recall_90p_',date,'.pdf']);
%             pause
        close(hf)
    end
end
top_g = 100;
s = sum(diffmat_2_recall>1,2);
[~,sxi] = sort(s,'descend');
topup = geneid(sxi(1:top_g))
figure('Position',[100,100,500,770],'color','w')
subplot(1,2,1)
barh(s(sxi(1:top_g)));
set(gca,'ytick',[1:top_g],'YTickLabel',geneid(sxi(1:top_g)),'fontsize',7)
title('upregulated, 2hr>recall')
ylabel('genes')
xlabel('#celltypes d>1')

s = sum(diffmat_2_recall<-1,2);
[~,sxi] = sort(s,'descend');
todown = geneid(sxi(1:top_g))
subplot(1,2,2)
barh(s(sxi(1:top_g)));
set(gca,'ytick',[1:top_g],'YTickLabel',geneid(sxi(1:top_g)),'fontsize',7)
title('downregulated, 2hr<recall')
ylabel('genes')
xlabel('#celltypes d<-1')
eval(['export_fig Glut12_topgenes_2hr_recall_90p_',date,'.pdf']);
table1 = [ {'upreg','downreg'};[topup,todown] ];
saveCellFile(table1,['topgenes_2hr_recall_90p_Glut12_',date,'.txt']);
%% 8hr vs 24hr
diffmat_2_recall = zeros(length(geneid),length(T_cells_tmp_uni));
for k=1:length(T_cells_tmp_uni)
    k
    gr1 = find((T_cells_tmp==k & fc_time==8 & batch_flag==2));
    gr2 = find((T_cells_tmp==k & fc_time==24 & batch_flag==2));
    if min([length(gr1),length(gr2)])>10
        hf = figure('Position',[100,100,800,800],'Color','w','Name',clusteruni{k,2});
        top_g = 50;
        x1 = log2(prctile(data(:,gr1),90,2)+1);
        x2 = log2(prctile(data(:,gr2),90,2)+1);
        % x1 = (prctile(data(:,gr1),90,2));
        % x2 = (prctile(data(:,gr2),90,2));
        d = x1-x2 ;
        diffmat_2_recall(:,k) = d;
        s = (x1+x2)/2 ;
        d(max([x1,x2],[],2)<3) = 0;
        [~,xi] = sort(d);

        plot((x1 + x2)/2, x1-x2, '.');hold on;
        xmax = max(x1);
        plot([0,xmax],[0,0],'-k'); grid on
        plot([0,xmax],[1,1]*0.4,'--k'); grid on
        plot([0,xmax],[1,1]*(-0.4),'--k'); grid on
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
        xi = flipud(xi);
        plot( (x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,'.r'); hold on;
        text((x1(xi(1:top_g)) + x2(xi(1:top_g)))/2 , x1(xi(1:top_g)) - x2(xi(1:top_g)) ,geneid(xi(1:top_g)),'fontsize',8);
        title(['gr1=',num2str(length(gr1)),',gr2=',num2str(length(gr2))])
        xlabel('(gr1+gr2)/2')
        ylabel(['2hr-recall'])
        axis tight
        eval(['export_fig Glut12_scatter_8hr_24hr_90p_220602/',clusteruni{k,2},'_scatter_8hr_24hr_90p_',date,'.pdf']);
        %     pause
        close(hf)
    end
end
top_g = 100;
s = sum(diffmat_2_recall>1,2);
[~,sxi] = sort(s,'descend');
topup = geneid(sxi(1:top_g))
figure('Position',[100,100,500,770],'color','w')
subplot(1,2,1)
barh(s(sxi(1:top_g)));
set(gca,'ytick',[1:top_g],'YTickLabel',geneid(sxi(1:top_g)),'fontsize',7)
title('upregulated, 2hr>recall')
ylabel('genes')
xlabel('#celltypes d>1')

s = sum(diffmat_2_recall<-1,2);
[~,sxi] = sort(s,'descend');
todown = geneid(sxi(1:top_g))
subplot(1,2,2)
barh(s(sxi(1:top_g)));
set(gca,'ytick',[1:top_g],'YTickLabel',geneid(sxi(1:top_g)),'fontsize',7)
title('downregulated, 2hr<recall')
ylabel('genes')
xlabel('#celltypes d<-1')
eval(['export_fig Glut12_topgenes_8hr_24hr_90p_',date,'.pdf']);
table1 = [ {'upreg','downreg'};[topup,todown] ];
saveCellFile(table1,['topgenes_8hr_24hr_90p_Glut12_',date,'.txt']);

%% gene correlation per cluster
list1 = {'Egr1'
    'Egr4'
    'Junb'
    'Scg2'
    'Ntrk2'
    'Vgf'
    'Dclk1'
    'Syt4'
    'Clstn3'
    'Pde10a'
    'Ptprn'
    'Zdbf2'
    'Syt4'
    'Nptxr'
    'Arpp21'
    'Sv2c'
    'Tiparp'
    'Wipf3'
    'Ptpn5'
    'Plppr4'
    'Sorcs3'
    'Chgb'
    'Nrxn3'
    'Nptx2'
    'Kcnv1'
    'Kcnf1'
    'Lgi2'};
ieggenes = [{'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'},list1'];
tmp = cell(length(ieggenes));
for i=1:length(ieggenes)
    for j=i+1:length(ieggenes)
        tmp{i,j} = [ieggenes{i},'-',ieggenes{j}];
    end
end
yt = tmp(~cellfun(@isempty,tmp));
[~,loc] = ismember(ieggenes,geneid);
pcth = 95;
d = data(loc,:);
rmat = zeros(length(T_cells_tmp_uni),length(ieggenes)*(length(ieggenes)-1)/2);
hf = figure('position',[100,100,800,800],'color','w');
ha = tight_subplot(8,10,0.02,0.02,0.02);
for k=1:length(T_cells_tmp_uni)
    gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time==0);
    r = corr_mat(d(:,gr1)');
    t = r.*sqrt((length(gr1)-2)./(1-r.^2));
    p = 2*(1-tcdf(abs(t),length(gr1)-2));
    
%     hf = figure('position',[100,100,800,800],'color','w');
    axes(ha(k))
%     imagesc(r,[0,0.5]);title(clusteruni{k,2},'Interpreter','none')
    imagesc(-log10(p+1e-5),[0,5]);title(clusteruni{k,2},'Interpreter','none')
    axis off
%     set(gca,'xtick',[1:length(ieggenes)],'xticklabel',ieggenes,'ytick',[1:length(ieggenes)],'yticklabel',ieggenes)
%     pause
%     close(hf)
%     r = triu(r,1);
%     r = r(r~=0);
%     rmat(k,:) = r';
%     hf = figure('Position',[200,200,800,800],'Color','w');
%     imagesc(r);
%     set(gca,'YTick',[1:length(ieggenes)],'YTickLabel',ieggenes,'XTick',[1:length(ieggenes)],'XTickLabel',ieggenes)
%     title(clusteruni{k,2},'Interpreter','none')
%     pause
%     close (hf)
end
 
hf = figure('Position',[200,200,800,1200],'Color','w');
imagesc(rmat)
set(gca,'YTick',[1:length(clusteruni)],'YTickLabel',clusteruni(:,2),'XTick',[1:length(yt)],'XTickLabel',yt,'TickLabelInterpreter','none')
%%
list1 = {'Egr1'
    'Egr4'
    'Junb'
    'Scg2'
    'Ntrk2'
    'Vgf'
    'Dclk1'
    'Syt4'
    'Clstn3'
    'Pde10a'
    'Ptprn'
    'Zdbf2'
    'Syt4'
    'Nptxr'
    'Arpp21'
    'Sv2c'
    'Tiparp'
    'Wipf3'
    'Ptpn5'
    'Plppr4'
    'Sorcs3'
    'Chgb'
    'Nrxn3'
    'Nptx2'
    'Kcnv1'
    'Kcnf1'
    'Lgi2'};
ieggenes = list1;
colorvec = distinguishable_colors(length(ieggenes));
[~,loc] = ismember(ieggenes,geneid);
pcth = 90;
d = data(loc,:);
d95 = d>repmat(prctile(d,pcth,2),1,length(d));
iegfram_mat = cell(length(T_cells_tmp_uni),1);
iegfram_mat_fold = cell(length(T_cells_tmp_uni),1);
for k=1:length(T_cells_tmp_uni)
    k
    if min(grsz_batch23(k,:),[],2)>0
        tmp =zeros(length(ieggenes)+1,5);
        for jj=1:5
            gr1 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag>=2));
            tmp(:,jj) = [mean(d95(:,gr1),2);mean(sum(d95(:,gr1))>0)];
        end
        iegfram_mat{k} = tmp;
        tmp1 = tmp(:,1);
        tmp1(tmp1==0) = 1/grsz_batch23(k,1);
        iegfram_mat_fold{k} = tmp./repmat(tmp1,1,length(fc_time_uni));
    end
end

figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    k = find(strcmpi(glut12_cluster_by_iegscore{i},clusteruni));
    axes(ha(i))
    for j=1:10%length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat_fold{k}(j,:),'-s','color',colorvec(j,:)); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    set(ha(i),'xtick',[0,2,8,24,26])
    title(regexprep(clusteruni{k,2},'_','-'))
end
legend(ieggenes(1:10))
figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    k = find(strcmpi(glut12_cluster_by_iegscore{i},clusteruni));
    axes(ha(i))
    for j=11:20%length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat_fold{k}(j,:),'-s','color',colorvec(j,:)); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    set(ha(i),'xtick',[0,2,8,24,26])
    title(regexprep(clusteruni{k,2},'_','-'))
end
legend(ieggenes(11:20))
figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    k = find(strcmpi(glut12_cluster_by_iegscore{i},clusteruni));
    axes(ha(i))
    for j=21:length(ieggenes)
        plot([0,2,8,24,26],iegfram_mat_fold{k}(j,:),'-s','color',colorvec(j,:)); hold on;
%         plot([1:5],iegfram_mat{xi(i)}(j,:),'-s'); hold on;
    end
    set(ha(i),'xtick',[0,2,8,24,26])
    title(regexprep(clusteruni{k,2},'_','-'))
end
legend(ieggenes(21:end))

%%
list1 = {'ACTR2','ABI2','CDK5','BAIAP2','PDLIM5','CFL1','LRRK2','ZDHHC15','DBN1','DHX36','EFNA1','EPHA4'....
    ,'EPHB2','NLGN1','ARC','CUX2','NGEF','SIPA1L1','DNM3','TANC2','STAU2','DBNL',....
    'EEF2K','ITPKA','CAPRIN1','PAFAH1B1','PAK3','ABI3','KIF1A','RELN','PTEN','CAPRIN2'....
    ,'UBE3A','LRP8','SRCIN1','BHLHB9','CAMK2B','DTNBP1','PPFIA2','SHANK3','CDK5R1','LZTS3','ARHGAP44'};
list1 = {'AVIL','CFL1','CFL2','DSTN','FLII','SRGAP2','GMFB','GSN','MYH9','VILL','SVIL','VIL1','FMNL1','CAPG','SCIN','GMFG'};
list1 = {'CDH2','ARMCX5-GPRASP2','ACTR2','ABI2','CDK5','LILRB2','BAIAP2','PDLIM5','CFL1','IL1RAPL1','LRRK2','ZDHHC15','DBN1','DHX36','EFNA1','EPHA4','EPHA7','EPHB2','FCGR2B','NLGN1','ARC','CUX2','FYN','NGEF','SIPA1L1','LRRTM2','DNM3','TANC2','GHSR','STAU2','GRID2','DBNL','GRIN2B','EEF2K','HSPA8','APOE','IL1RAP','INS','ITPKA','ARF4','MIR30B','CAPRIN1','NPTX1','NRCAM','NTRK3','PAFAH1B1','PAK3','ABHD17B','ABI3','GHRL','KIF1A','SRGN','VPS35','PRNP','RELN','PTEN','LRFN2','PTPRD','ABHD17C','ARHGAP22','SEMA3F','CAPRIN2','CDKL5','UBE3A','WNT5A','WNT7A','LRP8','LRFN4','SRCIN1','BHLHB9','CAMK2B','ABHD17A','DTNBP1','PPFIA2','SHANK3','CBLN1','NRP2','CDK5R1','ZNF804A','LRRC4B','CRIPT','HOMER1','LZTS3','ARHGAP44'};
ieggenes = list1;
colorvec = distinguishable_colors(length(ieggenes));
[~,loc] = ismember(upper(ieggenes),upper(geneid));
ieggenes = ieggenes(loc>0);
loc = loc(loc>0);
pcth = 90;
d = data(loc,:);
d95 = d>repmat(prctile(d,pcth,2),1,length(d));
iegfram_mat = cell(length(T_cells_tmp_uni),1);
iegfram_mat_fold = cell(length(T_cells_tmp_uni),1);
for k=1:length(T_cells_tmp_uni)
    k
    if min(grsz_batch23(k,:),[],2)>0
        tmp =zeros(length(ieggenes)+1,5);
        for jj=1:5
            gr1 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag>=2));
            tmp(:,jj) = [mean(d95(:,gr1),2);mean(sum(d95(:,gr1))>0)];
        end
        iegfram_mat{k} = tmp;
        tmp1 = tmp(:,1);
        tmp1(tmp1==0) = 1/grsz_batch23(k,1);
        iegfram_mat_fold{k} = tmp./repmat(tmp1,1,length(fc_time_uni));
    end
end

figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    k = find(strcmpi(glut12_cluster_by_iegscore{i},clusteruni));
    axes(ha(i))
    imagesc(iegfram_mat_fold{k},[0,3])
    if i==1
        set(gca,'xtick',[1:5],'xticklabel',[0,2,8,24,26],'ytick',[1:length(ieggenes)],'yticklabel',ieggenes);
    else
        set(gca,'xtick',[1:5],'xticklabel',[0,2,8,24,26],'ytick',[1:length(ieggenes)],'yticklabel',cell(size(ieggenes)));
    end
    colormap('summer')
    title(regexprep(clusteruni{k,2},'_','-'))
end

figure('position',[100,100,1400,400],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.02,0.02], [0.05,0.1], [0.05,0.02]);
for i=1:5
    k = find(strcmpi(glut12_cluster_by_iegscore{end-i},clusteruni));
    axes(ha(i))
    imagesc(iegfram_mat_fold{k},[0,3])
    if i==1
        set(gca,'xtick',[1:5],'xticklabel',[0,2,8,24,26],'ytick',[1:length(ieggenes)],'yticklabel',ieggenes);
    else
        set(gca,'xtick',[1:5],'xticklabel',[0,2,8,24,26],'ytick',[1:length(ieggenes)],'yticklabel',cell(size(ieggenes)));
    end
    colormap('summer')
    title(regexprep(clusteruni{k,2},'_','-'))
end

%%
semagenes = {'Sema4c' ,'Sema6d' ,'Sema4a' ,'Sema6c' ,'Sema3d' ,'Sema3a' ,'Sema3e' ,'Sema3c' ...
    ,'Sema4f' ,'Sema4b' ,'Sema3g' ,'Sema7a' ,'Sema3b' ,'Sema3f' ,'Sema4d' ,'Sema5a' ,'Sema5b' ...
    ,'Sema6b' ,'Sema6a' ,'Sema4g'};
memgenes = {'Creb1' ,'Creb3l1' ,'Creb3l4' ,'Creb3' ,'Creb3l2' ,'Creb5' ...
    ,'Crebl2' ,'Crebzf' ,'Creb3l3' ,'Crebbp' ,'Crebrf' ,'Nrxn3' ,'Nrxn1' ,'Nrxn2' ,'Nlgn3' ,'Nlgn1' ...
    ,'Nlgn2' ,'Nxph2' ,'Nxpe5' ,'Nxph1' ,'Nxph4' ,'Nxpe2' ,'Nxpe4' ,'Nxph3' ,'Nxpe3','Camk2a','Camk2b',...
    'Cttn','Dbn1','Pclo','Bsn','Syn1','Syn2','Syn3','Syp','Ntrk2','Dlg1','Dlg2','Dlg3','Dlg4'};
glutrec = {'Gria1','Gria2','Gria3'...
    ,'Gria4','Grin1','Grin2a','Grin2b','Grin2d','Grin3a','Grin2c','Grid1','Grid2','Grik1','Grik2','Grik3','Grik4',....
    'Grm1','Grm2','Grm3','Grm4','Grm5','Nsmf'};
gabarec = {'Gabra1','Gabra2','Gabra3','Gabra4','Gabra5','Gabra6','Gabrb1','Gabrb2'...
    ,'Gabrb3','Gabrb4','Gabrb5','Gabra6','Gabbr1','Gabbr2','Gabbr3','Gabrg1','Gabrg2','Gabrd','Gabre'};
potassium_channels = loadCellFile('potassium_channel_list.txt');
potassium_channels = lower(potassium_channels);
potassium_channels = cellfun(@(x) [upper(x(1)),x(2:end)], potassium_channels,'uniformoutput',0);
calcium_channels = loadCellFile('calcium_channel_list.txt');
calcium_channels = lower(calcium_channels);
calcium_channels = cellfun(@(x) [upper(x(1)),x(2:end)], calcium_channels,'uniformoutput',0);
% ieggenes = [{'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'},list1];
% ieggenes = [{'Bdnf','Arc','Npas4','Nlgn2' ,'Nrxn2','Creb3'}];
ieg = [{'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'},{'Fosb','Junb','Egr1','Egr2','Egr4','Nr4a3','Dusp1','Dusp6','Cdkn1a','Gadd45b','Inhba'}];
% tmp = cell(length(ieggenes));
% for i=1:length(ieggenes)
%     for j=i+1:length(ieggenes)
%         tmp{i,j} = [ieggenes{i},'-',ieggenes{j}];
%     end
% end
% yt = tmp(~cellfun(@isempty,tmp));
ieggenes = [ieg,semagenes,memgenes,gabarec,glutrec,potassium_channels',calcium_channels'];
[~,loc] = ismember(upper(ieggenes),upper(geneid));
ieggenes(loc==0) = [];
loc(loc==0) = [];
pcth = 95;
d = data(loc,:);
loc(sum(d>0,2)<200) = [];
ieggenes(sum(d>0,2)<200) = []; 
d = data(loc,:);

% rmat = zeros(length(T_cells_tmp_uni),length(ieggenes)*(length(ieggenes)-1)/2);
% hf = figure('position',[100,100,800,800],'color','w');
% ha = tight_subplot(7,8,0.02,0.02,0.02);
mean_p_r = zeros(length(T_cells_tmp_uni),6);
mean_p_r_ieg = zeros(length(T_cells_tmp_uni),6);
mean_p_r_mem = zeros(length(T_cells_tmp_uni),6);
mean_p_r_gabarec = zeros(length(T_cells_tmp_uni),6);
mean_p_r_glutrec = zeros(length(T_cells_tmp_uni),6);
mean_p_r_pc = zeros(length(T_cells_tmp_uni),6);
mean_p_r_cc = zeros(length(T_cells_tmp_uni),6);
mean_p_r_sema = zeros(length(T_cells_tmp_uni),6);
[~,locieg] = ismember(ieg,ieggenes);locieg(locieg==0) = [];
[~,locmem] = ismember(memgenes,ieggenes);locmem(locmem==0) = [];
[~,locsema] = ismember(semagenes,ieggenes);locsema(locsema==0) = [];
[~,locgabarec] = ismember(gabarec,ieggenes);locgabarec(locgabarec==0) = [];
[~,locglutrec] = ismember(glutrec,ieggenes);locglutrec(locglutrec==0) = [];
[~,locpc] = ismember(potassium_channels,ieggenes);locpc(locpc==0) = [];
[~,loccc] = ismember(calcium_channels,ieggenes);loccc(loccc==0) = [];

colorvec = distinguishable_colors(7);
validct = false(length(T_cells_tmp_uni),1);
bdnf_corrgenes = ones(length(geneid),length(T_cells_tmp_uni)*6);
nlgn2_corrgenes = ones(length(geneid),length(T_cells_tmp_uni)*6);
i_bdnf = find(strcmpi('Bdnf',ieggenes));
i_nlgn1 = find(strcmpi('Nlgn1',ieggenes));
for k=1:length(T_cells_tmp_uni)
    gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time>=0);
    bdnfcorr = false(length(ieggenes),6);
    nlgn2corr = false(length(ieggenes),6);
    if min(grsz_batch23(k,:))>20
        validct(k) = 1;
        r = corr_mat(d(:,gr1)');
        r(isnan(r)) = 0;
        t = r.*sqrt((length(gr1)-2)./(1-r.^2));
        p = 2*(1-tcdf(abs(t),length(gr1)-2));
        oo = 1:length(ieggenes);
        for vv = 1:10; [sorted_ind, ~, ~] = neighborhood_sort(r(oo,oo), 1,1); oo = oo(sorted_ind); end
        sorted_ind = oo;
        hf = figure('position',[80,200,1800,350],'color','w');
        ha = tight_subplot(1,6,0.04,0.06,0.06);
        %                 hf2 = figure('position',[80,200,1800,800],'color','w');
        %                 hb = tight_subplot(1,6,0.04,0.06,0.06);
        for jj=1:6
            gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time2==fc_time_uni2(jj));
            r = corr_mat(d(:,gr1)');
            r(isnan(r)) = 0;
            t = r.*sqrt((length(gr1)-2)./(1-r.^2));
            p = (1-tcdf((t),length(gr1)-2));
            inexp = find(sum(data(:,gr1)>0,2)>5);
            [~,prbdnf] = corr(d(i_bdnf,gr1)',data(inexp,gr1)');
            prbdnf(isnan(prbdnf)) = 1;
            [~,xibdnf] = sort(prbdnf);
            bdnf_corrgenes(inexp(xibdnf(1:sum(prbdnf<1e-2))),(k-1)*6+jj) = prbdnf(xibdnf(1:sum(prbdnf<1e-2)))';%1;
            [~,prnlgn2] = corr(d(i_nlgn1,gr1)',data(inexp,gr1)');
            prnlgn2(isnan(prnlgn2)) = 1;
            [~,xinlgn2] = sort(prnlgn2);          
            nlgn2_corrgenes(inexp(xinlgn2(1:sum(prnlgn2<1e-2))),(k-1)*6+jj) = prnlgn2(xinlgn2(1:sum(prnlgn2<1e-2)))';%1;

            %             xytmp = tsne_vdM(d(:,gr1), [], 2, 30, 50, [500],2,1000);
            % %             D = 1-r;
            % %             P = d2p(D, perplexity, 1e-5);
            % %             P = p<0.01;
            % %             xytmp = tsne_p(P, [], 2,500,1000);
            %
            %             axes(hb(jj))
            %             plot(xytmp(locieg,1),xytmp(locieg,2),'.','color',colorvec(1,:));hold on;
            %             text(xytmp(locieg,1),xytmp(locieg,2),ieggenes(locieg),'fontsize',6);
            %             plot(xytmp(locmem,1),xytmp(locmem,2),'.','color',colorvec(2,:));hold on;
            %             text(xytmp(locmem,1),xytmp(locmem,2),ieggenes(locmem),'fontsize',6);
            %             plot(xytmp(locgabarec,1),xytmp(locgabarec,2),'.','color',colorvec(3,:));hold on;
            %             text(xytmp(locgabarec,1),xytmp(locgabarec,2),ieggenes(locgabarec),'fontsize',6);
            %             plot(xytmp(locglutrec,1),xytmp(locglutrec,2),'.','color',colorvec(4,:));hold on;
            %             text(xytmp(locglutrec,1),xytmp(locglutrec,2),ieggenes(locglutrec),'fontsize',6);
            %             plot(xytmp(locpc,1),xytmp(locpc,2),'.','color',colorvec(5,:));hold on;
            %             text(xytmp(locpc,1),xytmp(locpc,2),ieggenes(locpc),'fontsize',6);
            %             plot(xytmp(loccc,1),xytmp(loccc,2),'.','color',colorvec(6,:));hold on;
            %             text(xytmp(loccc,1),xytmp(loccc,2),ieggenes(loccc),'fontsize',6);
            %             plot(xytmp(locsema,1),xytmp(locsema,2),'.','color',colorvec(7,:));hold on;
            %             text(xytmp(locsema,1),xytmp(locsema,2),ieggenes(locsema),'fontsize',6);
            bdnfcorr(:,jj) = p(:,i_bdnf)<0.01;
            nlgn2corr(:,jj) = p(i_nlgn1,:)<0.01;

            axes(ha(jj))
            %     hf = figure('position',[100,100,800,800],'color','w');
            imagesc(-log10(p(sorted_ind,sorted_ind)+1e-5),[0,5]);title(clusteruni{k,2},'Interpreter','none')
            set(gca,'xtick',[1:length(ieggenes)],'xticklabel',ieggenes(sorted_ind),'ytick',[1:length(ieggenes)],'yticklabel',ieggenes(sorted_ind))
            colormap("summer")
%             colormap(cm)
            mean_p_r(k,jj) = mean(-log10(p(triu(ones(size(p)),1)>0)+1e-10));

            ptmp = p(locieg,locieg);
            mean_p_r_ieg(k,jj) = mean(-log10(ptmp(triu(ones(size(ptmp)),1)>0)+1e-10));
            ptmp = p(locmem,locmem);
            mean_p_r_mem(k,jj) = mean(-log10(ptmp(triu(ones(size(ptmp)),1)>0)+1e-10));
            ptmp = p(locsema,locsema);
            mean_p_r_sema(k,jj) = mean(-log10(ptmp(triu(ones(size(ptmp)),1)>0)+1e-10));
            ptmp = p(locgabarec,locgabarec);
            mean_p_r_gabarec(k,jj) = mean(-log10(ptmp(triu(ones(size(ptmp)),1)>0)+1e-10));
            ptmp = p(locglutrec,locglutrec);
            mean_p_r_glutrec(k,jj) = mean(-log10(ptmp(triu(ones(size(ptmp)),1)>0)+1e-10));
            ptmp = p(locpc,locpc);
            mean_p_r_pc(k,jj) = mean(-log10(ptmp(triu(ones(size(ptmp)),1)>0)+1e-10));
            ptmp = p(loccc,loccc);
            mean_p_r_cc(k,jj) = mean(-log10(ptmp(triu(ones(size(ptmp)),1)>0)+1e-10));
        end
        bdnfcorr_cell{k} = bdnfcorr;
        nlgn2corr_cell{k} = nlgn2corr;
        colorbar('position',[0.97,0.02,0.01,0.9])
        linkaxes(ha,'xy')
%         eval(['export_fig Glut12_IEGmemgenes_correlation_220712/',clusteruni{k,2},'_corrP_IEGmemG_p3_',date,'.pdf']);
        
%         hf2 = figure('position',[80,200,800,800],'color','w');
%         [tmp,xi] = sortrows(bdnfcorr,'descend');
%         intmp = find(sum(tmp,2)>0);
%         subplot(1,2,1);
%         imagesc(tmp(intmp,:));
%         colormap('gray')
%         set(gca,'ytick',[1:length(intmp)],'yticklabel',ieggenes(xi(intmp)),'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'})
%         title('Bdnf')
%         [tmp,xi] = sortrows(nlgn2corr,'descend');
%         intmp = find(sum(tmp,2)>0);
%         subplot(1,2,2);
%         imagesc(tmp(intmp,:))
%         set(gca,'ytick',[1:length(intmp)],'yticklabel',ieggenes(xi(intmp)),'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'})
%         title('Nlgn2')
%         eval(['export_fig Glut12_IEGmemgenes_correlation_220704/',clusteruni{k,2},'_Bdnf_Nlgn2_IEGmemG_',date,'.pdf']);
%         pause
        %     eval(['export_fig Glut12_IEGmemgenes_correlation_220606/',clusteruni{k,2},'_corrP_IEGmemG_',date,'.pdf']);
        close(hf)
%         close(hf2)

    %     axes(ha(k))
    % %     imagesc(r,[0,0.5]);title(clusteruni{k,2},'Interpreter','none')
    %     imagesc(-log10(p+1e-3),[0,3]);title(clusteruni{k,2},'Interpreter','none')
    %     axis off
    %     set(gca,'xtick',[1:length(ieggenes)],'xticklabel',ieggenes,'ytick',[1:length(ieggenes)],'yticklabel',ieggenes)
    %     pause
    %     close(hf)
    %     r = triu(r,1);
    %     r = r(r~=0);
    %     rmat(k,:) = r';
    %     hf = figure('Position',[200,200,800,800],'Color','w');
    %     imagesc(r);
    %     set(gca,'YTick',[1:length(ieggenes)],'YTickLabel',ieggenes,'XTick',[1:length(ieggenes)],'XTickLabel',ieggenes)
    %     title(clusteruni{k,2},'Interpreter','none')
    %     pause
    %     close (hf)
    end
end


% x = cell2mat(bdnfcorr_cell(validct));
% tmp = x;
% intmp = find(sum(tmp,2)>4);
% tmp = tmp(intmp,[1:6:end,2:6:end,3:6:end,4:6:end,5:6:end,6:6:end]);
% r = corr_mat(tmp');
% oo = 1:length(intmp);
% for vv = 1:10; [sorted_ind, ~, ~] = neighborhood_sort(r(oo,oo), 1,1); oo = oo(sorted_ind); end
% figure
% imagesc(tmp(sorted_ind,:));
% colormap('gray')
% set(gca,'ytick',[1:length(intmp)],'yticklabel',ieggenes(intmp(sorted_ind)),'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'})
% title('Bdnf')
% 
% x = cell2mat(nlgn2corr_cell(validct));
% tmp = x;
% intmp = find(sum(tmp,2)>4);
% tmp = tmp(intmp,[1:6:end,2:6:end,3:6:end,4:6:end,5:6:end,6:6:end]);
% r = corr_mat(tmp');
% oo = 1:length(intmp);
% for vv = 1:10; [sorted_ind, ~, ~] = neighborhood_sort(r(oo,oo), 1,1); oo = oo(sorted_ind); end
% figure
% imagesc(tmp(sorted_ind,:));
% colormap('gray')
% set(gca,'ytick',[1:length(intmp)],'yticklabel',ieggenes(intmp(sorted_ind)),'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'})
% title('nlgn2')

score_all = max(mean_p_r(:,3:end),[],2)-max(mean_p_r(:,1:2),[],2);
score_ieg = max(mean_p_r_ieg(:,3:end),[],2)-max(mean_p_r_ieg(:,1:2),[],2);
score_mem = max(mean_p_r_mem(:,3:end),[],2)-max(mean_p_r_mem(:,1:2),[],2);
score_gabarec = max(mean_p_r_gabarec(:,3:end),[],2)-max(mean_p_r_gabarec(:,1:2),[],2);
score_glutrec = max(mean_p_r_glutrec(:,3:end),[],2)-max(mean_p_r_glutrec(:,1:2),[],2);
score_pc = max(mean_p_r_pc(:,3:end),[],2)-max(mean_p_r_pc(:,1:2),[],2);
score_cc = max(mean_p_r_cc(:,3:end),[],2)-max(mean_p_r_cc(:,1:2),[],2);
score_sema = max(mean_p_r_sema(:,3:end),[],2)-max(mean_p_r_sema(:,1:2),[],2);

table1 = [ [{'score_all','score_ieg','score_mem','score_gabarec','score_glutrec','score_pc','score_cc','score_sema'}];
    m2c([score_all,score_ieg,score_mem,score_gabarec,score_glutrec,score_pc,score_cc,score_sema])];
table1 = [ [{''};clusteruni(:,2)],table1 ];
% saveCellFile(table1,['corrscore_Glut12_',date,'.txt'])

validct_in = find(validct);
[~,xiscoreall] = sort(score_all(validct_in),'descend');
validct_in_sort = validct_in(xiscoreall);
cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);
figure('position',[400,100,1050,350],'color','w');
[ha, pos] = tight_subplot(1, 8, [0.01,0.01], [0.03,0.08], [0.18,0.01]);
axes(ha(1))
imagesc(mean_p_r(validct_in_sort,:))
colormap(cm)
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',clusteruni(validct_in_sort,2),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('All ');
axes(ha(2))
imagesc(mean_p_r_ieg(validct_in_sort,:))
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',cell(sum(validct_in_sort),1),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('IEG ');
axes(ha(3))
imagesc(mean_p_r_mem(validct_in_sort,:))
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',cell(sum(validct_in_sort),1),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('Mem ');
axes(ha(4))
imagesc(mean_p_r_gabarec(validct_in_sort,:))
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',cell(sum(validct_in_sort),1),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('GABArec ');
axes(ha(5))
imagesc(mean_p_r_glutrec(validct_in_sort,:))
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',cell(sum(validct_in_sort),1),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('GLUTrec ');
axes(ha(6))
imagesc(mean_p_r_pc(validct_in_sort,:))
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',cell(sum(validct_in_sort),1),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('K channels ');
axes(ha(7))
imagesc(mean_p_r_cc(validct_in_sort,:))
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',cell(sum(validct_in_sort),1),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('Ca channels ');
axes(ha(8))
imagesc(mean_p_r_sema(validct_in_sort,:))
set(gca,'ytick',[1:length(T_cells_tmp_uni(validct_in_sort))],'yticklabel',cell(sum(validct_in_sort),1),'ticklabelinterpreter','none',...
    'xtick',[1:6],'xticklabel',{'0','0','2','8','24','recall'});
title('SEMA ');
% eval(['export_fig Corrscore_Glut12_',date,'.pdf']);

% % 
[~,xiscoreall] = sort(score_all(validct),'descend');

bdnf_corrgenes_bin = bdnf_corrgenes<1e-2;
intmp = find(sum(bdnf_corrgenes_bin(:,[3:6:end,4:6:end,5:6:end,6:6:end]),2)>=15);
validct_in = find(validct);
intime = [(validct_in(xiscoreall)-1)*6+repmat([1:6],sum(validct),1)];
intime = intime(:);
% [~,xi] = sortrows(bdnf_corrgenes(intmp,intime([1:6:end,2:6:end,3:6:end,4:6:end,5:6:end,6:6:end])));
% [~,xi] = sortrows(bdnf_corrgenes(intmp,intime));
[~,xi] = sort(sum(bdnf_corrgenes_bin(intmp,[3:6:end,4:6:end,5:6:end,6:6:end]),2),'descend');
figure('position',[80,200,800,800],'color','w');
ax1 = axes('position',[0.1,0.1,0.8,0.8]);
imagesc(-log10(bdnf_corrgenes(intmp(xi),intime)),[0,3]); hold on;
colormap('summer')
% for i=1:sum(validct)
%     plot(i*6*[1,1]+0.5,[0,length(intmp)],'-r')
% end
for i=1:5
    plot(i*sum(validct)*[1,1]+0.5,[0,length(intmp)],'-r')
end
% set(gca,'ytick',[1:length(intmp)],'yticklabel',geneid(intmp(xi)),....
%     'xtick',[1:length(intime)],'xticklabel',repmat({'0','0','2','8','24','recall'},1,sum(validct)))
set(gca,'ytick',[1:length(intmp)],'yticklabel',geneid(intmp(xi)),....
    'xtick',[1:length(intime)],'xticklabel',cell(length(intime),1))
ax2 = axes('position',[0.9,0.1,0.1,0.8]);
barh(sum(bdnf_corrgenes(intmp(xi),intime)<1e-2,2))
axis tight
box off
set(gca,'ydir','reverse','ytick',[1:length(intmp)],'yticklabel',cell(length(intmp),1),'xlim',[0,sum(bdnf_corrgenes(intmp(xi(2)),:)<1e-2)])
ax3 = axes('position',[0.1,0.9,0.8,0.1]);
bar(sum(bdnf_corrgenes(intmp(xi),intime)<1e-2))
axis tight
set(gca,'xtick',[1:length(intime)],'xticklabel',cell(length(intime),1))
box off
% eval(['export_fig Bdnf_15_Corr_Glut12_',date,'.pdf']);
sort(geneid(sum(bdnf_corrgenes(:,[3:6:end,4:6:end,5:6:end,6:6:end]),2)>5))
intmp = find(sum(bdnf_corrgenes_bin(:,[3:6:end,4:6:end,5:6:end,6:6:end]),2)>=5);
[z,xi] = sort(sum(bdnf_corrgenes_bin(intmp,[3:6:end,4:6:end,5:6:end,6:6:end]),2),'descend');
saveCellFile([geneid(intmp(xi)),m2c(z)],['genes_corr_with_Bdnf_Glut12_',date,'.txt'])

ngln2_corrgenes_bin = nlgn2_corrgenes<1e-2;
intmp = find(sum(ngln2_corrgenes_bin(:,[3:6:end,4:6:end,5:6:end,6:6:end]),2)>=10);
validct_in = find(validct);
intime = [(validct_in(xiscoreall)-1)*6+repmat([1:6],sum(validct),1)];
intime = intime(:);
% [~,xi] = sortrows(ngln2_corrgenes(intmp,intime([1:6:end,2:6:end,3:6:end,4:6:end,5:6:end,6:6:end])));
% [~,xi] = sortrows(ngln2_corrgenes(intmp,intime));
[~,xi] = sort(sum(ngln2_corrgenes_bin(intmp,[3:6:end,4:6:end,5:6:end,6:6:end]),2),'descend');
figure('position',[80,200,800,800],'color','w');
ax1 = axes('position',[0.1,0.1,0.8,0.8]);
imagesc(-log10(nlgn2_corrgenes(intmp(xi),intime)),[0,3]); hold on;
colormap('summer')
% for i=1:sum(validct)
%     plot(i*6*[1,1]+0.5,[0,length(intmp)],'-r')
% end
for i=1:5
    plot(i*sum(validct)*[1,1]+0.5,[0,length(intmp)],'-r')
end
% set(gca,'ytick',[1:length(intmp)],'yticklabel',geneid(intmp(xi)),....
%     'xtick',[1:length(intime)],'xticklabel',repmat({'0','0','2','8','24','recall'},1,sum(validct)))
set(gca,'ytick',[1:length(intmp)],'yticklabel',geneid(intmp(xi)),....
    'xtick',[1:length(intime)],'xticklabel',cell(length(intime),1))
ax2 = axes('position',[0.9,0.1,0.1,0.8]);
barh(sum(nlgn2_corrgenes(intmp(xi),intime)<1e-2,2))
axis tight
box off
set(gca,'ydir','reverse','ytick',[1:length(intmp)],'yticklabel',cell(length(intmp),1),'xlim',[0,sum(nlgn2_corrgenes(intmp(xi(2)),:)<1e-2)])
ax3 = axes('position',[0.1,0.9,0.8,0.1]);
bar(sum(nlgn2_corrgenes(intmp(xi),intime)<1e-2))
axis tight
box off
set(gca,'xtick',[1:length(intime)],'xticklabel',cell(length(intime),1))
% eval(['export_fig Nlgn1_10_Corr_Glut12_',date,'.pdf']);
sort(geneid(sum(nlgn2_corrgenes(:,[3:6:end,4:6:end,5:6:end,6:6:end]),2)>5))
intmp = find(sum(ngln2_corrgenes_bin(:,[3:6:end,4:6:end,5:6:end,6:6:end]),2)>=2);
[z,xi] = sort(sum(ngln2_corrgenes_bin(intmp,[3:6:end,4:6:end,5:6:end,6:6:end]),2),'descend');
saveCellFile([geneid(intmp(xi)),m2c(z)],['genes_corr_with_Nlgn1_Glut12_',date,'.txt'])
%%
for k=1:length(T_cells_tmp_uni)
    genelist = {'Snap25','Actg1','Actb','Cfl1'};
    c = k;
    hf = figure('color','w','position',[20,20,877,180],'Name',regexprep(clusteruni{c,2},'_','-'));
    [ha, pos] = tight_subplot(1, length(genelist), [0.05,0.05], [0.08,0.08], [0.05,0.05]);
    b = 1;
    logflag = 0;
    fc_time2 = fc_time;
    fc_time2(fc_time==0 & batch_flag==2)= 1;
    fc_time_uni2 = unique(fc_time2);
    p99 = [3];%prctile(data(g,:),99.99);
    colorvec = [0,0,0;44,34,96;63,162,201]/255;
    for kkk=1:length(genelist)
        p99 = [3];
        axes(ha(kkk))
        gn = genelist{kkk};
        g = find(strcmpi(gn,geneid));

        t_ed = zeros(5,1);
        t_av = zeros(5,1);
        t_75 = zeros(5,1);
        p90 = prctile(data(g,T_cells_tmp(:,b)==c & batch_flag>=2),90);
        if sum(T_cells_tmp(:,b)==c & batch_flag>=2)>100
            for i=[1:6]
                gr2 = find(T_cells_tmp(:,b)==c & fc_time2==fc_time_uni2(i) & batch_flag>=2);%find(fc_time_sorted==fc_time_uni(i));%
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
                    for jj=2:3
                        v = bfgr2==jj;
                        plot(x1(v),y1(v) ,'.','color',colorvec(jj,:)); hold on;
                    end
                else
                    plot(i+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
                end

                %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
                t_ed(i) = median(data(g,gr2));
                t_av(i) = mean(data(g,gr2));
                t_75(i) = prctile(data(g,gr2),90);
                f90(i) = sum(data(g,gr2)>p90);
                pp90(i) = binocdf(f90(i),length(gr2),0.1,'upper');
                p99 = max([p99,prctile(data(g,gr2),97)]);
            end
            plot([0,6],p90*[1,1],'--k');
            plot(t_ed,'sk'); hold on;
            %     plot(t_av,'-or'); hold on;
            plot(t_75,'dm'); hold on;
            axis tight
            yl = get(gca,'ylim');
            set(gca,'xtick',[1:6],'XTickLabel',[{'0','0','2','8','24','recall'}],'yscale','linear','ylim',[-0.5,p99+1],'YScale',yscale,'ytick',round(linspace(0,p99,4)),...
                'yTickLabel',cellfun(@num2str,m2c(round(linspace(0,p99,4))),'UniformOutput',0))

            for ii=1:5
                if -log10(pp90(ii))>1
                    fs = 10;
                else
                    fs = 6;
                end
                text(ii,p99*0.7,cellfun(@(x) num2str(x,2), m2c(-log10(pp90(ii))),'UniformOutput',0),'fontsize',fs);
            end
            title(gn);%[gn,',c=',num2str(c)])
        end
    end
    eval(['export_fig Glut12_Cfl1_Actb_220606/',clusteruni{k,2},'_Cfl1_Actb_',date,'.pdf']);
    close(hf)

end
%%
list1 = {'Sema4c' ,'Sema6d' ,'Sema4a' ,'Sema6c' ,'Sema3d' ,'Sema3a' ,'Sema3e' ,'Sema3c' ...
    ,'Sema4f' ,'Sema4b' ,'Sema3g' ,'Sema7a' ,'Sema3b' ,'Sema3f' ,'Sema4d' ,'Sema5a' ,'Sema5b' ...
    ,'Sema6b' ,'Sema6a' ,'Sema4g' ,'Creb1' ,'Creb3l1' ,'Creb3l4' ,'Creb3' ,'Creb3l2' ,'Creb5' ...
    ,'Crebl2' ,'Crebzf' ,'Creb3l3' ,'Crebbp' ,'Crebrf' ,'Nrxn3' ,'Nrxn1' ,'Nrxn2' ,'Nlgn3' ,'Nlgn1' ...
    ,'Nlgn2' ,'Nxph2' ,'Nxpe5' ,'Nxph1' ,'Nxph4' ,'Nxpe2' ,'Nxpe4' ,'Nxph3' ,'Nxpe3'};
ieggenes = [{'Bdnf','Arc' ,'Npas4','Nlgn2','Nrxn2','Creb3'}];
[~,loc] = ismember(upper(ieggenes),upper(geneid));
ieggenes(loc==0) = [];
loc(loc==0) = [];
d = data(loc,:);
bdnf = data(strcmpi('bdnf',geneid),:);
tmp = cell(length(ieggenes));
for i=1:length(ieggenes)
    for j=i+1:length(ieggenes)
        tmp{j,i} = [ieggenes{i},'-',ieggenes{j}];
    end
end
yt = tmp(~cellfun(@isempty,tmp));

% rmat = zeros(length(T_cells_tmp_uni),length(ieggenes)*(length(ieggenes)-1)/2);
% hf = figure('position',[100,100,800,800],'color','w');
% ha = tight_subplot(7,8,0.02,0.02,0.02);
colorvec = distinguishable_colors(length(yt));
for k=1:length(T_cells_tmp_uni)
    topbdnf = [];
    gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time>=0);
    if length(gr1)>100
        r = corr_mat(d(:,gr1)');
        r(isnan(r)) = 0;
        t = r.*sqrt((length(gr1)-2)./(1-r.^2));
        p = (1-tcdf((t),length(gr1)-2));
        oo = 1:length(ieggenes);
        %     for vv = 1:10; [sorted_ind, ~, ~] = neighborhood_sort(r(oo,oo), 1,1); oo = oo(sorted_ind); end
        sorted_ind = oo;
        %     hf = figure('position',[80,200,1800,410],'color','w');
        %     ha = tight_subplot(1,5,0.04,0.04,0.04);
        pmat = zeros(6,length(yt));
        for jj=1:6
%             if jj==1
%                 gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag==2 & fc_time2==fc_time_uni2(jj));
%             else
                gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time2==fc_time_uni2(jj));
%             end
            if length(gr1)>10
                in = find(sum(data(:,gr1),2)>5);
             rbdnf = corr(bdnf(gr1)',data(in,gr1)');
             [~,xibd] = sort(rbdnf,'descend');  
             t = rbdnf.*sqrt((length(gr1)-2)./(1-rbdnf.^2));
             pbdnf = (1-tcdf((t),length(gr1)-2));
             topbdnf = [topbdnf,[geneid(in(xibd(1:10))),m2c(-log10(1e-20+pbdnf(xibd(1:10)))')]];
            r = corr_mat(d(:,gr1)');
            r(isnan(r)) = 0;
            t = r.*sqrt((length(gr1)-2)./(1-r.^2));
            p = 2*(1-tcdf(abs(t),length(gr1)-2));
            p = -log10(squareform(tril(p,-1),'tovector')+1e-5);
%             p = triu(p,1);
%             p = -log10(p(p~=0));
            pmat(jj,:) = p';
            %         axes(ha(jj))
            %     hf = figure('position',[100,100,800,800],'color','w');
            %         imagesc(-log10(p(sorted_ind,sorted_ind)+1e-3),[0,3]);title(clusteruni{k,2},'Interpreter','none')
            %         set(gca,'xtick',[1:length(ieggenes)],'xticklabel',ieggenes(sorted_ind),'ytick',[1:length(ieggenes)],'yticklabel',ieggenes(sorted_ind))
            %         colormap("summer")
            end
        end
        topbdnf
        hf = figure('position',[100,100,800,800],'color','w');
        for jj=1:length(yt)
            if max(pmat(:,jj))>2
                plot(pmat(:,jj),'-o','color',colorvec(jj,:)); hold on;
            end
        end
        set(gca,'xtick',[1:6],'xticklabel',{'0b3','0b2','2','8','24','recall'}),
        title(clusteruni{k,2},'interpreter','none')
        legend(yt(max(pmat)>2))
        pause
        close(hf)
    end
end

%%

ieggenes = [{'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'}];
% tmp = cell(length(ieggenes));
% for i=1:length(ieggenes)
%     for j=i+1:length(ieggenes)
%         tmp{i,j} = [ieggenes{i},'-',ieggenes{j}];
%     end
% end
% yt = tmp(~cellfun(@isempty,tmp));
[~,loc] = ismember(upper(ieggenes),upper(geneid));
ieggenes(loc==0) = [];
loc(loc==0) = [];
pcth = 95;
d = data(loc,:);
% rmat = zeros(length(T_cells_tmp_uni),length(ieggenes)*(length(ieggenes)-1)/2);
hf = figure('position',[100,100,800,800],'color','w');
ha = tight_subplot(8,9,0.02,0.02,0.02);
for k=1:length(T_cells_tmp_uni)
    gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time==2);
    if length(gr1)>10
    r = corr_mat(d(:,gr1)');
    r(isnan(r)) = 0;
    t = r.*sqrt((length(gr1)-2)./(1-r.^2));
    p = 2*(1-tcdf(abs(t),length(gr1)-2));
%     oo = 1:length(ieggenes);
% %     for vv = 1:10; [sorted_ind, ~, ~] = neighborhood_sort(r(oo,oo), 1,1); oo = oo(sorted_ind); end
%     sorted_ind = oo;
%     hf = figure('position',[80,200,1800,410],'color','w');
%     ha = tight_subplot(1,6,0.04,0.04,0.04);
%     for jj=1:6
%         gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time2==fc_time_uni2(jj));
%         r = corr_mat(d(:,gr1)');
%         r(isnan(r)) = 0;
%         t = r.*sqrt((length(gr1)-2)./(1-r.^2));
%         p = (1-tcdf((t),length(gr1)-2));
%         axes(ha(jj))
%         %     hf = figure('position',[100,100,800,800],'color','w');
%         imagesc(-log10(p(sorted_ind,sorted_ind)+1e-3),[0,3]);title(clusteruni{k,2},'Interpreter','none')
%         set(gca,'xtick',[1:length(ieggenes)],'xticklabel',ieggenes(sorted_ind),'ytick',[1:length(ieggenes)],'yticklabel',ieggenes(sorted_ind))
%         colormap("summer")
%     end
%     colorbar('position',[0.97,0.02,0.01,0.9])
%     pause
%     eval(['export_fig Glut12_IEGmemgenes_correlation_220606/',clusteruni{k,2},'_corrP_IEGmemG_',date,'.pdf']);
%     close(hf)
        mean_p_r(k) = mean(-log10(p(triu(ones(size(p)),1)>0)+1e-10));
        axes(ha(k))
    %     imagesc(r,[0,0.5]);title(clusteruni{k,2},'Interpreter','none')
        imagesc(-log10(p+1e-10),[0,10]);title(clusteruni{k,2},'Interpreter','none')
        axis off
        set(gca,'xtick',[1:length(ieggenes)],'xticklabel',ieggenes,'ytick',[1:length(ieggenes)],'yticklabel',ieggenes)
%         pause
%         close(hf)
%         r = triu(r,1);
%         r = r(r~=0);
%         rmat(k,:) = r';
%         hf = figure('Position',[200,200,800,800],'Color','w');
%         imagesc(r);
        set(gca,'YTick',[1:length(ieggenes)],'YTickLabel',ieggenes,'XTick',[1:length(ieggenes)],'XTickLabel',ieggenes)
        title(clusteruni{k,2},'Interpreter','none')
        colormap("summer")
%         pause
%         close (hf)
    end
end

%%
ieg = [{'Bdnf','Homer1','Fos','Arc','Btg2','Fosl2','Npas4','Nr4a1'},{'Fosb','Junb','Egr1','Egr2','Egr4','Nr4a3','Dusp1','Dusp6','Cdkn1a','Gadd45b','Inhba'}];
ieggenes = [ieg];
[~,loc] = ismember(upper(ieggenes),upper(geneid));
ieggenes(loc==0) = [];
loc(loc==0) = [];
pcth = 95;
d = data(loc,:);
loc(sum(d>0,2)<200) = [];
ieggenes(sum(d>0,2)<200) = []; 
d = data(loc,:);

colorvec = distinguishable_colors(7);
validct = false(length(T_cells_tmp_uni),1);

for k=1:length(T_cells_tmp_uni)
    gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time>=0);
    bdnfcorr = false(length(ieggenes),6);
    nlgn2corr = false(length(ieggenes),6);
    if min(grsz_batch23(k,:))>20
        validct(k) = 1;
        r = corr_mat(d(:,gr1)');
        r(isnan(r)) = 0;
        t = r.*sqrt((length(gr1)-2)./(1-r.^2));
        p = 2*(1-tcdf(abs(t),length(gr1)-2));
        oo = 1:length(ieggenes);
        for vv = 1:10; [sorted_ind, ~, ~] = neighborhood_sort(r(oo,oo), 1,1); oo = oo(sorted_ind); end
        sorted_ind = oo;
        hf = figure('position',[80,200,1800,350],'color','w');
        ha = tight_subplot(1,6,0.04,0.06,0.06);
        %                 hf2 = figure('position',[80,200,1800,800],'color','w');
        %                 hb = tight_subplot(1,6,0.04,0.06,0.06);
        for jj=1:6
            gr1 = find(T_cells_tmp == T_cells_tmp_uni(k) & batch_flag>=2 & fc_time2==fc_time_uni2(jj));
            r = corr_mat(d(:,gr1)');
            r(isnan(r)) = 0;
            t = r.*sqrt((length(gr1)-2)./(1-r.^2));
            p = (1-tcdf((t),length(gr1)-2));

            axes(ha(jj))
            %     hf = figure('position',[100,100,800,800],'color','w');
            imagesc(-log10(p(sorted_ind,sorted_ind)+1e-5),[0,5]);title(clusteruni{k,2},'Interpreter','none')
            set(gca,'xtick',[1:length(ieggenes)],'xticklabel',ieggenes(sorted_ind),'ytick',[1:length(ieggenes)],'yticklabel',ieggenes(sorted_ind))
            colormap("summer")           
        end
        
        colorbar('position',[0.97,0.02,0.01,0.9])
        linkaxes(ha,'xy')
        eval(['export_fig Glut12_IEG_correlation_220712/',clusteruni{k,2},'_corrP_IEG_',date,'.pdf']);        
        close(hf)
    end
end