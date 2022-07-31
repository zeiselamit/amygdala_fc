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

cellidclsuter_gaba = loadCellFile('cellidCluster_step3_GABA_FC_24-Oct-2021.txt');
% cellidclsuter_gaba = [cellidclsuter_gaba, cellfun(@(x,y) ['gaba_',num2str(x),'_',num2str(y)],cellidclsuter_gaba(:,2),cellidclsuter_gaba(:,3),'UniformOutput',0)]; 

gaba_clustername = loadCellFile('ClusterOrder_step3_GABA_FC_24-Oct-2021.txt');
cluster_info = loadCellFile('Clusters_info_hannah.txt');
cluster_info = cluster_info(2:57,:);

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
    T_cells_tmp(ind) = i;
    cluster_loc_flag(ind) = repmat(cluster_location(i),length(ind),1);
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
        tmp3 = (tmp2+1)./(tmp+1);
%         fc50 = (max(tmp(in,:),[],2)+1)./(min(tmp(in,:),[],2)+1);
        fc50 = tmp(in,3)-tmp(in,1);
%         fc90 = tmp2(in,3)-tmp2(in,1);%(max(tmp2(in,:),[],2)+1)./(min(tmp2(in,:),[],2)+1);%
        fc90 = tmp2(in,3)-tmp2(in,1);
        [~,xi] = sort(fc90,'descend');
        topgene_b2 = [topgene_b2; [geneid(in(xi(1:top_g))),m2c([fc90(xi(1:top_g)),fc50(xi(1:top_g)),repmat(k,top_g,1)])]];
    end
end

top_g = 10;
topgene_b23 = [];
fc90_mat = zeros(length(geneid),length(T_cells_tmp_uni));
for k=1:length(T_cells_tmp_uni)
    k
    if min(grsz_batch23(k,1))>30
        tmp = p50_batch23_cluster_tp{k};
%         in = find(max(tmp,[],2)>0);
        tmp2 = p90_batch23_cluster_tp{k};
        in = find(max(tmp2,[],2)>4);
        tmp3 = (tmp2+1)./(tmp+1);
        tmp3 = tmp3./repmat(tmp3(:,1),1,5);
%         fc50 = (max(tmp(in,:),[],2)+1)./(min(tmp(in,:),[],2)+1);
        fc50 = tmp(in,3)-tmp(in,1);
%         fc90 = tmp2(in,3)-tmp2(in,1);%(max(tmp2(in,:),[],2)+1)./(min(tmp2(in,:),[],2)+1);%
        fc90 = max(tmp3(in,:),[],2);%tmp2(in,3)./tmp2(in,1);
        fc90_mat(in,k) = fc90;
        [~,xi] = sort(fc90,'descend');
        topgene_b23 = [topgene_b23; [geneid(in(xi(1:top_g))),m2c([fc90(xi(1:top_g)),fc50(xi(1:top_g)),repmat(k,top_g,1)])]];
    end
end

g = unique(topgene_b23(:,1));
[~,gidx] = ismember(g,geneid);
tmp = cellfun(@(x) x(gidx,:), p90_batch23_cluster_tp,'UniformOutput',0);
tmp = cellfun(@(x) x./repmat(x(:,1),1,5), tmp,'UniformOutput',0);
tmp = cell2mat(tmp);
tmp(isnan(tmp)) = 1;
tmp(isinf(tmp)) = 1;
[leaforder] = get1Dorderlinkage(tmp,1);
n = length(tmp(1,:));
tmp = tmp(:,[1:5:n,2:5:n,3:5:n,4:5:n,5:5:n]);
figure;
tmpcn = tmp;%cent_norm(tmp);%
imagesc(tmpcn(leaforder,:),[prctile(tmpcn(:),2),prctile(tmpcn(:),98)]); hold on;
for i=56:56:n
    plot(i+0.5*[1,1],[0,length(g)],'k');
end
set(gca,'ytick',[1:length(g)],'yticklabel',g(leaforder))

figure;
tmpcn = tmp;%cent_norm(tmp);%
imagesc(tmpcn(leaforder,:),[prctile(tmpcn(:),2),prctile(tmpcn(:),98)]); hold on;
for i=5:5:n
    plot(i+0.5*[1,1],[0,length(g)],'k');
end
set(gca,'ytick',[1:length(g)],'yticklabel',g(leaforder))

%%
gn = 'dgkb';
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
% p99 = prctile(data(g,:),99);
for k=1:length(T_cells_tmp_uni)
    if min(grsz_batch23(k,1))>30
        axes(ha(k))
        b23_50 = p50_batch23_cluster_tp{k}(g,:)+1;
        b23_90 = p90_batch23_cluster_tp{k}(g,:)+1;       
        
        plot([1,2,3,4,5],b23_50,'--sb'); hold on;
        plot([1,2,3,4,5],b23_90,'-sc'); hold on;      
        
        axis tight
        yl = get(gca,'ylim');
        if k>=(length(T_cells_tmp_uni)-ceil(length(unique(T_cells_tmp(:,b)))/8))
            set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear','YScale',yscale)
        else
            set(gca,'xtick',[1:5],'XTickLabel',cell(1,5),'yscale','linear','YScale',yscale)
        end
        
        title(regexprep(clusteruni{k},'_','-'));%[gn,',c=',num2str(c)])
    end
end
linkaxes(ha,'y');
%%
gn = 'syngr1';
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
p99 = [];%prctile(data(g,:),99.99);
colorvec = distinguishable_colors(3);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    c=k;
    t_ed = zeros(5,1);
    t_av = zeros(5,1);
    t_75 = zeros(5,1);
    p90 = prctile(data(g,T_cells_tmp(:,b)==c & batch_flag>=2),90);
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
        f90(i) = sum(data(g,gr2)>p90);
        pp90(i) = binocdf(f90(i),length(gr2),0.1,'upper');
        p99 = max([p99,prctile(data(g,gr2),99)]);
    end
    plot([0,6],p90*[1,1],'--k');
    %     plot(t_ed,'-sk'); hold on;
    %     plot(t_av,'-or'); hold on;
    %     plot(t_75,'-dg'); hold on;
    axis tight
    yl = get(gca,'ylim');
    if c>(length(T_cells_tmp_uni)-7)
        set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear','ylim',[-0.5,p99+1],'YScale',yscale,'ytick',round(linspace(0,p99,4)),...
            'yTickLabel',cellfun(@num2str,m2c(round(linspace(0,p99,4))),'UniformOutput',0))
    else
        set(gca,'xtick',[1:5],'XTickLabel',cell(1,5),'yscale','linear'...
            ,'ylim',[-0.5,p99+1],'YScale',yscale,'ytick',round(linspace(0,p99,4)),'yTickLabel',cellfun(@num2str,m2c(round(linspace(0,p99,4))),'UniformOutput',0))
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
linkaxes(ha,'y');
if 1==100
    %     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
    eval(['export_fig ',gn,'_violin_FC_GABA_',date,'.pdf']);    
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
saveCellFile([hitgenes],['hitgenes_gaba_',date,'.txt'])

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
save(['pbino_gaba_',date],'pbino_mat','geneid','clusteruni')
save(['p90mat_gaba_',date],'p90_mat','geneid','clusteruni')
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

saveCellFile(g,'top300_mostdiff_FC_GABA.txt')
tp_order = reshape([1:5*length(T_cells_tmp_uni)],5,length(T_cells_tmp_uni))';
tp_order = tp_order(:);
msum = reshape(sum(m(:,tp_order)>=1.5),length(T_cells_tmp_uni),5);
figure('position',[200,200,1200,1000],'color','w');
[ha, pos] = tight_subplot(1, 5, [0.05,0.05], [0.02,0.02], [0.05,0.05]);
for i=1:length(fc_time_uni)
    s = msum(:,i);
    [ssort,sxi] = sort(s,'descend');
    axes(ha(i));
    scatter(ssort,[1:length(T_cells_tmp_uni)],clustersize(sxi));
    set(ha(i),'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',regexprep(clusteruni(sxi),'_','-'),'fontsize',6,'xlim',[0,80])
    title(['t=',num2str(fc_time_uni(i))])
end

msum = reshape(sum(m_ieg(:,tp_order)>1),length(T_cells_tmp_uni),5);
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
    title(regexprep(clusteruni(k),'_','-'));
    set(gca,'ylim',[0,100])
end
    

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
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',g(leaforder),'fontsize',10,'xtick'...
    ,[length(T_cells_tmp_uni)/2:length(T_cells_tmp_uni):5*length(T_cells_tmp_uni)],'XTickLabel',{'0','2','8','24','recall'})
if savefig_flag==100
    %     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
    eval(['export_fig top100_mostdiff_FC_GABA_',date,'.pdf']);    
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
%%
% % % % % % % % % % % % % % % % % % % % % % % % 
in = find(sum(data>0,2)>100 & sum(data>0,2)<length(data(1,:))*0.9);
ave_flag = 1;
b = 1;
c_top = cell(length(unique(T_cells_tmp(:,b))),1);
c_top2 = cell(length(unique(T_cells_tmp(:,b))),1);
c_top3 = cell(length(unique(T_cells_tmp(:,b))),1);
c_top_rnd = cell(length(unique(T_cells_tmp(:,b))),1);
fc_time_rnd = fc_time_sorted;
for c1=1:length(unique(T_cells_tmp(:,b)))
    c1
    fc_time_uni = [0,2,8,24,28];
    top_g = 20;
    t = fc_time_rnd(T_cells_tmp(:,b)==c1);
    fc_time_rnd(T_cells_tmp(:,b)==c1) = t(randperm(length(t)));
    gr1 = find(T_cells_tmp(:,b)==c1 & fc_time_sorted==fc_time_uni(1));
    gr1_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(1));
    x1 = mean(log2(data(:,gr1)+1),2);
    x1_rnd = mean(log2(data(:,gr1_rnd)+1),2);
    p90gene = prctile(data(:,T_cells_tmp(:,b)==c1),90,2);
    k = 0;
    qmat = ones(length(geneid_all),4);
    qmat2 = ones(length(geneid_all),4);
    qmat3 = ones(length(geneid_all),4);
    qmat_rnd = ones(length(geneid_all),4);
    for i=[2,3,4,5]
        k = k+1;
        gr2 = find(T_cells_tmp(:,b)==c1 & fc_time_sorted==fc_time_uni(i));
        gr2_not = find(T_cells_tmp(:,b)==c1 & fc_time_sorted~=fc_time_uni(i));
        gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(i));
        x2 = mean(log2(data(:,gr2)+1),2);
        x2_rnd = mean(log2(data(:,gr2_rnd)+1),2);
        d = x1-x2 ;
        x90 = sum(data(:,gr2)>repmat(p90gene,1,length(gr2)),2);
        x90(p90gene==0) = nan;
        pbino = binocdf(x90,length(gr2),0.1,'upper');
        pbino(isnan(pbino)) = 1;
        %         [~,ptt] = ttest2(data(in,gr1)',data(in,gr2)','tail','both');
        ptt = ones(length(in),1);
        ptt2 = ones(length(in),1);
        ptt_rnd = ones(length(in),1);
        if length(gr1)>10 & length(gr2)>10
            for s=1:length(in)
                %             s
                %             ptt(s) = ranksum(data(in(s),gr1),data(in(s),gr2),'tail','left');
                [~,ptt(s)] = kstest2(data(in(s),gr1),data(in(s),gr2),'tail','larger');
                [~,ptt2(s)] = kstest2(data(in(s),gr2_not),data(in(s),gr2),'tail','larger');
                [~,ptt_rnd(s)] = kstest2(data(in(s),gr1_rnd),data(in(s),gr2_rnd),'tail','larger');
            end
        end
        ptt(isnan(ptt)) = 1;
        ptt2(isnan(ptt2)) = 1;
        ptt_rnd(isnan(ptt_rnd)) = 1;
        %         in2 = (in(ptt<1));
        %         q = qval_from_pval(ptt(ptt<1));
        %         q(q<1e-10) = 1e-10;
        qmat(in,k) = ptt;
        qmat2(in,k) = ptt2;
        qmat3(in,k) = pbino(in);
        qmat_rnd(in,k) = ptt_rnd;
    end
    c_av_tp = zeros(length(geneid_all),length(fc_time_uni));
    c_av_tp_rnd = zeros(length(geneid_all),length(fc_time_uni));
    rank_tp = zeros(length(geneid_all),length(fc_time_uni));
    [~,xii] = sort(data(:,T_cells_tmp(:,b)==c1)');
    tp_c1 = fc_time_sorted(T_cells_tmp(:,b)==c1);
    xii = xii';
    grsize = zeros(5,1);
    for i=[1,2,3,4,5]
        gr2 = find(T_cells_tmp(:,b)==c1 & fc_time_sorted==fc_time_uni(i));
        gr2_rnd = find(T_cells_tmp(:,b)==c1 & fc_time_rnd==fc_time_uni(i));
        grsize(i) = length(gr2);
        if ave_flag==1
            c_av_tp(:,i) = median(data(:,gr2),2);
            c_av_tp_rnd(:,i) = median(data(:,gr2_rnd),2);
        else
            c_av_tp(:,i) = mean(data(:,gr2)>0,2);
        end
        rank_tp(:,i) = mean(xii(:,tp_c1==fc_time_uni(i)),2)/length(tp_c1);
    end
    
    [qsort,xi] = sort(qmat);
    [qsort2,xi2] = sort(qmat2);
    [qsort3,xi3] = sort(qmat3);
    [qsort_rnd,xi_rnd] = sort(qmat_rnd);
    meddiff = c_av_tp(:,2:5)-repmat(c_av_tp(:,1),1,4);
    meddiff_rnd = c_av_tp_rnd(:,2:5)-repmat(c_av_tp_rnd(:,1),1,4);
    c_top{c1} = [xi(1:200,:),qsort(1:200,:),[meddiff(xi(1:200,1),1),meddiff(xi(1:200,2),2),meddiff(xi(1:200,3),3),meddiff(xi(1:200,4),4)]];
    c_top2{c1} = [xi2(1:200,:),qsort2(1:200,:),[meddiff(xi(1:200,1),1),meddiff(xi(1:200,2),2),meddiff(xi(1:200,3),3),meddiff(xi(1:200,4),4)]];
    c_top3{c1} = [xi3(1:200,:),qsort3(1:200,:),[meddiff(xi(1:200,1),1),meddiff(xi(1:200,2),2),meddiff(xi(1:200,3),3),meddiff(xi(1:200,4),4)]];
    c_top_rnd{c1} = [xi_rnd(1:200,:),qsort_rnd(1:200,:),[meddiff_rnd(xi(1:200,1),1),meddiff_rnd(xi(1:200,2),2),meddiff_rnd(xi(1:200,3),3),meddiff_rnd(xi(1:200,4),4)]];
  
end


%%
% [ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp(:,1),data_sorted_all,1);
list =  {'Ptk2b','Wfs1','Six3','Satb1','Tshz2','Ppp1r1b','Pax6','Foxp2','Prkcd','Adora2a','Drd1','Pdyn','Tac1','Tac2','Gbx1','Nts','Vdr','Lhx6','Prlr'};
% list = {'Meis2','Maf','Zeb2','Prlr','Lhx6','Npas1','Igfbp4','Sst','Cryab','Penk','Moxd1','Prkcd','Unc5b','Ucn3','Pthlh','Vip','Chat','Htr3a','Vdr'};
% list = {'Fth1','Ftl1','Scrib','Tfr2','Tfrc','Slc11a2','Slc40a1','Aco1','Ireb2'};
figure;
set(gcf,'color','w','position',[20,20,1000,1000]);
[ha, pos] = tight_subplot(length(list), 1, [0.01,0.01], [0.02,0.01], [0.09,0.01]);
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
    eval(['export_fig strAmy_markers_violin_final_step3_FC_GABA_',date,'.pdf']);
    
end


%%
% [ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp(:,1),data_sorted_all,1);
list =  {'Maf','Nfib','Zeb2','Htr3a','Sncg','Vip','Chat','Reln','Ndnf','Lamp5','Pvalb','Sst','Chodl'};%geneid(ind_gr_tmp_mark(1:end));
% list = {'Meis2','Maf','Zeb2','Prlr','Lhx6','Npas1','Igfbp4','Sst','Cryab','Penk','Moxd1','Prkcd','Unc5b','Ucn3','Pthlh','Vip','Chat','Htr3a','Vdr'};
% list = {'Fth1','Ftl1','Scrib','Tfr2','Tfrc','Slc11a2','Slc40a1','Aco1','Ireb2'};
figure;
set(gcf,'color','w','position',[20,20,1000,1000]);
[ha, pos] = tight_subplot(length(list), 1, [0.01,0.01], [0.02,0.01], [0.09,0.01]);
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
    eval(['export_fig INT_markers_violin_final_step3_FC_GABA_',date,'.pdf']);
end


%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
markers = {'Sema4c' ,'Sema6d' ,'Sema4a' ,'Sema6c' ,'Sema3d' ,'Sema3a' ,'Sema3e' ,'Sema3c' ...
    ,'Sema4f' ,'Sema4b' ,'Sema3g' ,'Sema7a' ,'Sema3b' ,'Sema3f' ,'Sema4d' ,'Sema5a' ,'Sema5b' ...
    ,'Sema6b' ,'Sema6a' ,'Sema4g' ,'Creb1' ,'Creb3l1' ,'Creb3l4' ,'Creb3' ,'Creb3l2' ,'Creb5' ...
    ,'Crebl2' ,'Crebzf' ,'Creb3l3' ,'Crebbp' ,'Crebrf' ,'Nrxn3' ,'Nrxn1' ,'Nrxn2' ,'Nlgn3' ,'Nlgn1' ...
    ,'Nlgn2' ,'Nxph2' ,'Nxpe5' ,'Nxph1' ,'Nxph4' ,'Nxpe2' ,'Nxpe4' ,'Nxph3' ,'Nxpe3'};
% markers = IEG_markers;
t1 = [];
t2 = [];
for ii=1:length(markers)
    tmpmark = markers{ii};
    marker = data(strcmpi(geneid,tmpmark),:);
    marker_percent = zeros(1, length(T_cells_tmp_uni));
    p_marker = zeros(1, length(T_cells_tmp_uni));
    M = length(marker);
    K = sum(marker>0);
    for j=1:length(T_cells_tmp_uni)
        c1 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j));
        c2 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j) & marker'>0);
        marker_percent(j) = 100*c2/c1;
        p_marker(j) = hygecdf(c2,M,K,c1,'upper');
    end
    marker_percent(isnan(marker_percent)) = 0;
    eval([markers{ii},'_percent = marker_percent;'])
    eval([markers{ii},'_p = p_marker;']);
    eval(['t1 = [t1;',markers{ii},'_percent];'])
    eval(['t2 = [t2;',markers{ii},'_p];'])
end
figure;
set(gcf,'color','w','position',[20,20,900,1200])
[ha, pos] = tight_subplot(1, 2, [0.05,0.00], [0.05,0.05], [0.05,0.02]);
axes(ha(1));
imagesc(t1');
colormap('summer');
set(gca,'ytick',[1:length(T_cells_tmp_uni)],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
axes(ha(2));
imagesc(-log10(t2)',[0,10]);
colormap('summer');
set(gca,'ytick',[],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
linkaxes([ha(1),ha(2)],'y');


%%


markers = IEG_markers;
t1 = [];
t2 = [];
for ii=1:length(markers)
    tmpmark = markers{ii};
    marker = data(strcmpi(geneid,tmpmark),:);
    marker_percent = zeros(1, length(T_cells_tmp_uni));
    p_marker = zeros(1, length(T_cells_tmp_uni));
    M = length(marker);
    K = sum(marker>0);
    for j=1:length(T_cells_tmp_uni)
        c1 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j));
        c2 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j) & marker'>0);
        marker_percent(j) = 100*c2/c1;
        p_marker(j) = hygecdf(c2,M,K,c1,'upper');
    end
    marker_percent(isnan(marker_percent)) = 0;
    eval([markers{ii},'_percent = marker_percent;'])
    eval([markers{ii},'_p = p_marker;']);
    eval(['t1 = [t1;',markers{ii},'_percent];'])
    eval(['t2 = [t2;',markers{ii},'_p];'])
end
figure;
set(gcf,'color','w','position',[20,20,900,1200])
[ha, pos] = tight_subplot(1, 2, [0.05,0.00], [0.05,0.05], [0.05,0.02]);
axes(ha(1));
imagesc(t1');
colormap('summer');
set(gca,'ytick',[1:length(T_cells_tmp_uni)],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
axes(ha(2));
imagesc(-log10(t2)',[0,10]);
colormap('summer');
set(gca,'ytick',[],'xtick',[1:length(markers)],'XTickLabel',markers,'XTickLabelRotation',45,'FontSize',8)
linkaxes([ha(1),ha(2)],'y');
eval(['export_fig IEG_heatmap_per_cluster_',date,'.pdf']);
%%
markers = {'fos','arc','npas4','rrad'};
t1 = [];
t2 = [];
t3 = [];
for ii=1:length(markers)
    tmpmark = markers{ii};
    marker = data(strcmpi(geneid,tmpmark),:);
    marker_percent = zeros(1, length(T_cells_tmp_uni));
    marker_ave = zeros(1, length(T_cells_tmp_uni));
    p_marker = zeros(1, length(T_cells_tmp_uni));
    M = length(marker);
    K = sum(marker>0);
    for j=1:length(T_cells_tmp_uni)
        c1 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j));
        c2 = sum( T_cells_tmp(:,1)==T_cells_tmp_uni(j) & marker'>0);
        marker_percent(j) = 100*c2/c1;
        marker_ave(j) = mean(marker(T_cells_tmp(:,1)==T_cells_tmp_uni(j)));
        p_marker(j) = hygecdf(c2,M,K,c1,'upper');
    end
    marker_percent(isnan(marker_percent)) = 0;
    marker_ave(isnan(marker_ave)) = 0;
    eval([markers{ii},'_percent = marker_percent;'])
    eval([markers{ii},'_ave = marker_ave;'])
    eval([markers{ii},'_p = p_marker;']);
    eval(['t1 = [t1;',markers{ii},'_percent];'])
    eval(['t2 = [t2;',markers{ii},'_ave];'])
    eval(['t3 = [t3;',markers{ii},'_p];'])
end

figure;
set(gcf,'color','w','position',[20,20,900,1800])
[ha, pos] = tight_subplot(1, 3, [0.02,0.00], [0.02,0.02], [0.05,0.02]);
axes(ha(1));
barh(t1');
set(gca,'ytick',[1:length(T_cells_tmp_uni)])
title('%>0')

axes(ha(2));
barh(t2');
set(gca,'ytick',[1:length(T_cells_tmp_uni)],'YTickLabel',cell(length(T_cells_tmp_uni),1),'XLim',[0,16])
title('Average')
axes(ha(3));
barh(-log10(t3)');
set(gca,'ytick',[1:length(T_cells_tmp_uni)],'YTickLabel',cell(length(T_cells_tmp_uni),1),'XLim',[0,10])
legend(markers);
title('-log10(p)')
eval(['export_fig Fos_arc_npas4_per_cluster_',date,'.pdf']);
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
[ha, pos] = tight_subplot(1, 3, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
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
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
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
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight
% eval(['export_fig scatter_volcano_',gr1name,'_',gr2name,'_',clusteruni((j)),'_',date,'.pdf']);

%%

% % % % % % % % % % % % %
list = {'Sst','Vip','Htr3a','Pvalb','Cnr1','Foxp2','Prlr','Lypd1','Meis2','Rrad','Npas4','Arc','Jun','Fos'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_orig_all(strcmpi(geneid_all,genePlot),:));
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
if savefig_flag==1
    % eval(['export_fig tsne_markers2_AmyPiri_FC_',date,'.pdf']);
    savefig(gcf,['tsne_FC_GABA_markers1_',date,'.fig'])
end
if savefig_pdf==1
    eval(['export_fig tsne_FC_GABA_markers1_',date,'.pdf']);
end
%%


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
%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
top_g = 50;
gr1name = '0';
gr2name = '2';
sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d','Erdr1','Mid1',};%,'Snca','Gm9843','Gm12027','Gm12031'
gr1 = find((T_cells_tmp==9 & fc_time==0 & batch_flag==3));
gr2 = find((T_cells_tmp==9 & fc_time==2 & batch_flag==3));
% gr1 = find((strcmpi(cluster_loc_flag,'cge-int') & fc_time==0 & batch_flag==2));
% gr2 = find((strcmpi(cluster_loc_flag,'cge-int') & fc_time==24 & batch_flag==2));
% in = find(mean(data_orig_all_sorted(:,gr1)>0,2)>0.03 | mean(data_orig_all_sorted(:,gr2)>0,2)<0.7);
in = find(mean(data>0,2)>0.05 & mean(data>0,2)<0.5 & ~ismember(geneid(:,1),sex_genes));
ptt = zeros(length(in),1);
for s=1:length(in)
    if mod(s,100)==0
        s
    end
%     ptt(s) = ranksum(data(in(s),gr1),data(in(s),gr2),'tail','both');
    [~,ptt(s)] = kstest2(data(in(s),gr1),data(in(s),gr2));
end
ptt(isnan(ptt)) = 1;
ptt(ptt<1) = qval_from_pval(ptt(ptt<1));
ptt(ptt<1e-300) = 1e-300;
x1 = median(log2(data(:,gr1)+1),2);
x2 = median(log2(data(:,gr2)+1),2);
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

top_g = 100;
x1 = mean(data(:,gr1)>0,2);
x2 = mean(data(:,gr2)>0,2);
d = x1-x2 ;
d(max([x1,x2],[],2)<0.2) = 0;
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

top_g = 150;
x1 = log2(prctile(data(:,gr1),90,2)+1);
x2 = log2(prctile(data(:,gr2),90,2)+1);
d = x1-x2 ;
d(max([x1,x2],[],2)<3) = 0;
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
top_g = 200;
x1 = log2(prctile(data(:,gr1),90,2)+1);
x2 = log2(prctile(data(:,gr2),90,2)+1);
% x1 = (prctile(data(:,gr1),90,2));
% x2 = (prctile(data(:,gr2),90,2));
d = x1-x2 ;
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
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight
% eval(['export_fig scatter_volcano_',gr1name,'_',gr2name,'_',clusteruni((j)),'_',date,'.pdf']);
%%
for c=1:length(clusteruni)
    c
    % c = 3;
    list = {'Arc','Fos','Fosb','Fosl2','Npas4','Btg2','Jun','Junb','Egr1','Egr2','Egr4','Nr4a1','Nr4a3','Dusp1'};
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
            set(gca,'xtick',[1:5],'XTickLabel',[{'0','2','8','24','recall'}],'yscale','linear','YScale',yscale,'xlim',[1,5],'ylim',[0,9])
            
            %[gn,',c=',num2str(c)])
            if jjj==1
                ylabel(clusteruni{c},'Interpreter','none','fontsize',7)
                title(gn);
%                 title([clusteruni{c},',',gn],'Interpreter','none');
            else
                title(gn);
            end
            if jjj==length(list)
                text([1,2,3,4,5],8-[linspace(0,2,5)],cellfun(@num2str, m2c(grsz_batch23(c,:)),'UniformOutput',0),'fontsize',6)
            end
        end
    end
    linkaxes(ha,'xy');
    eval(['export_fig GABA_IEG_FC_220411/',clusteruni{c},'_IEG_FC',date,'.pdf']);
    close(hf)
end
%%
list = {'Arc','Fos','Fosb','Fosl2','Npas4','Btg2','Jun','Junb','Egr1','Egr2','Egr4','Nr4a1','Nr4a3','Dusp1'};
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
    set(ha(i),'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',cell(1,length(T_cells_tmp_uni)),'xlim',[1,10],'TickLabelInterpreter','none','fontsize',6)
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
set(gca,'ytick',[1:length(T_cells_tmp_uni)],'yticklabel',cell(1,length(T_cells_tmp_uni))...
    ,'fontsize',6,'ydir','reverse','ylim',[0.5,length(T_cells_tmp_uni)+0.5])




