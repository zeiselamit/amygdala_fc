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
fc_time_uni = unique(fc_time);
% % % % % % % % % % % % % % % % % % % 
%%


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
gn = 'gsk1';
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
        
        title(regexprep(clusteruni{k,2},'_','-'));%[gn,',c=',num2str(c)])
    end
end
linkaxes(ha,'y');
%%
%%
% Cdkn1a Nrgn Fth1 Sgk1 Gstp1 Mt1 Mt2 Ntm Negr1 Hes1 Hes5 Ccl3 Ccl4 Atf3
% Nfkbia C1qc Ctsd Fos Btg2 Jun Egr1 Cenpf Cenpa Cdc20 Ccnb1 Ccnd2
gn = 'cx3cr1';
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
p99 = [3];%prctile(data(g,:),99.99);
colorvec = distinguishable_colors(3);
fc_time2 = fc_time;
fc_time2(fc_time==0 & batch_flag==2)= 1;
fc_time_uni2 = unique(fc_time2);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    c=k;
    t_ed = zeros(5,1);
    t_av = zeros(5,1);
    t_75 = zeros(5,1);
    p90 = prctile(data(g,T_cells_tmp(:,b)==c & batch_flag==2),90);
    for i=[1:6]
        gr2 = find(T_cells_tmp(:,b)==c & fc_time2==fc_time_uni2(i) & batch_flag==2);%find(fc_time_sorted==fc_time_uni(i));%
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
            for jj=1:3
                v = bfgr2==jj;
                if sum(v)>10
            [f,xi] = ksdensity(y(v));
            fi = interp1(xi,f,y(v));
            fi = fi/max(fi);
%             if i==1 & jj==2
%             x1 = i-1 + fi'.*(0.9*rand(length(gr2(v)),1)-0.45);
%             y1 = 0.5*rand(length(gr2(v)),1)-0.1+y(v)';
%             else
                x1 = i + fi'.*(0.9*rand(length(gr2(v)),1)-0.45);
            y1 = 0.5*rand(length(gr2(v)),1)-0.1+y(v)';
%             end
            
            plot(x1,y1 ,'.','color',colorvec(jj,:)); hold on;
                end
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
        p99 = max([p99,prctile(data(g,gr2),99)]);
    end
    plot([0,6],p90*[1,1],'--k');
        plot(t_ed,'sk'); hold on;
    %     plot(t_av,'-or'); hold on;
        plot(t_75,'dg'); hold on;
    axis tight
    yl = get(gca,'ylim');
    if c>(length(T_cells_tmp_uni)-7)
        set(gca,'xtick',[1:6],'XTickLabel',[{'0','0','2','8','24','recall'}],'yscale','linear','ylim',[-0.5,p99+1],'YScale',yscale,'ytick',round(linspace(0,p99,4)),...
            'yTickLabel',cellfun(@num2str,m2c(round(linspace(0,p99,4))),'UniformOutput',0))
    else
        set(gca,'xtick',[1:6],'XTickLabel',cell(1,5),'yscale','linear'...
            ,'ylim',[-0.5,p99+1],'YScale',yscale,'ytick',round(linspace(0,p99,4)),'yTickLabel',cellfun(@num2str,m2c(round(linspace(0,p99,4))),'UniformOutput',0))
    end
    for ii=1:6
        if -log10(pp90(ii))>1
            fs = 10;
        else
            fs = 6;
        end
        text(ii,p99*0.7,cellfun(@(x) num2str(x,2), m2c(-log10(pp90(ii))),'UniformOutput',0),'fontsize',fs);
    end
    title(regexprep(clusteruni{c,1},'_','-'));%[gn,',c=',num2str(c)])
end
linkaxes(ha,'y');
if 1==100
    %     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
%     eval(['export_fig ',gn,'_violin_FC_GABA_',date,'.pdf']);    
end
%%
%%  scatters per time point median

top_g = 150;
gr1name = 'ctrl';
gr2name = '';
for k=1:length(T_cells_tmp_uni)
    k
    sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d','Erdr1','Mid1',};%,'Snca','Gm9843','Gm12027','Gm12031'
    if min(grsz_batch23(k,:),[],2)>20
        hf = figure('position',[200,200,1400,1000],'color','w','name',clusteruni{k});
        [ha, pos] = tight_subplot(4, 3, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
        s = 1;
        for jj=2:5
            if jj==2 | jj==5
                gr1 = find((T_cells_tmp==k & fc_time==0 & batch_flag==3));
                gr2 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag==3));
            else
                gr1 = find((T_cells_tmp==k & fc_time==0 & batch_flag==2));
                gr2 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag==2));
            end
            if min([length(gr1),length(gr2)])>5
                in = find(mean(data>0,2)>0.05 & mean(data>0,2)<0.7 & ~ismember(geneid(:,1),sex_genes));
                ptt = zeros(length(in),1);
                for z=1:length(in)
                    ptt(z) = ranksum(data(in(z),gr1),data(in(z),gr2),'tail','both');
                end
                ptt(isnan(ptt)) = 1;
                ptt(ptt<1) = qval_from_pval(ptt(ptt<1));
                ptt(ptt<1e-300) = 1e-300;
                ptmp = ones(length(geneid),1);
                ptmp(in) = ptt;
                %                 p_top_pertype(:,k) = -log10(ptmp(indcorrg));
                x1 = mean(log2(data(:,gr1)+1),2);
                x2 = mean(log2(data(:,gr2)+1),2);
                d = x2-x1 ;

                axes(ha(s)); s=s+1;
                plot(d(in),-log10(ptt),'.'); hold on;
                [~,xi] = sort(ptt);
                plot(d(in(xi(1:100))),-log10(ptt(xi(1:100))),'or');
                text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid(in(xi(1:100))))
                xlabel([num2str(fc_time_uni(jj)),'-',gr1name]);
                ylabel(['-log10(q)']);
                grid on;
                title(['comparing ',gr1name,' to ',num2str(fc_time_uni(jj)),', gr1=',num2str(length(gr1))....
                    ,', gr2=',num2str(length(gr2))])
                [~,xi] = sort(d);
                axes(ha(s));s=s+1;

                %
                d(max([x1,x2],[],2)<1) = 0;
                [~,xi] = sort(d);
                %                 d_top_pertype(:,k) = d(indcorrg);
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
                ylabel('gr2-gr1')
                axis tight
                yl = get(gca,'ylim');
                yl(2) = max(yl(2),1);yl(1) = min(yl(1),-1);
                set(gca,'ylim',yl);

                top_g = 100;
                x1 = mean(data(:,gr1)>0,2);
                x2 = mean(data(:,gr2)>0,2);
                d = x2-x1 ;
                d(max([x1,x2],[],2)<0.2) = 0;
                [~,xi] = sort(d);
                axes(ha(s));s=s+1;
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
            else
                s = s+3;
            end
            
        end
        linkaxes(ha([3,6,9,12]),'y')
        linkaxes(ha([2,5,8,11]),'y')
%             eval(['export_fig nonNeurons_median_vs_ctrl_220510/',clusteruni{k},'_scatter_median_FC',date,'.pdf']);
        pause
        close(hf)
    end

end


%%  scatters per time point 90th percentime

top_g = 150;
gr1name = 'ctrl';
gr2name = '';
for k=1:length(T_cells_tmp_uni)
    k
    sex_genes = {'Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d','Erdr1','Mid1',};%,'Snca','Gm9843','Gm12027','Gm12031'
    if min(grsz_batch23(k,:),[],2)>20
        hf = figure('position',[200,200,1400,600],'color','w','name',clusteruni{k});
        [ha, pos] = tight_subplot(1, 4, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
        s = 1;
        for jj=2:5
            if jj==2 | jj==5
                gr1 = find((T_cells_tmp==k & fc_time==0 & batch_flag==3));
                gr2 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag==3));
            else
                gr1 = find((T_cells_tmp==k & fc_time==0 & batch_flag==2));
                gr2 = find((T_cells_tmp==k & fc_time==fc_time_uni(jj) & batch_flag==2));
            end
            top_g = 100;
                x1 = prctile(log2(data(:,gr1)+1),90,2);
                x2 = prctile(log2(data(:,gr2)+1),90,2);
                d = x2-x1;
                d(max([x1,x2],[],2)<3) = 0;
                [~,xi] = sort(d);
                axes(ha(s));s=s+1;
                %
                plot((x1 + x2)/2, x2-x1, '.');hold on;
                xmax = max(x1);
                plot([0,xmax],[0,0],'-k'); grid on
                plot([0,xmax],[0,0]+0.5,'--k'); grid on
                plot([0,xmax],[0,0]-0.5,'--k'); grid on
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
                yl = [min(d((x1+x2)>4)),max(d((x1+x2)>4))];
                set(gca,'ylim',yl);
                title(['t=',num2str(fc_time_uni(jj)),', gr1=',num2str(length(gr1))....
                    ,', gr2=',num2str(length(gr2))])
        end
        linkaxes(ha,'y')
%         eval(['export_fig nonNeurons_prc90_vs_ctrl_220510/',clusteruni{k},'_scatter_median_FC',date,'.pdf']);
        pause
        close(hf)
    end
end

%%
k = 7;
% for k=1:length(T_cells_tmp_uni)
%     genelist = {'Cx3cr1','P2ry12','C1qb','Ctsd','Ccl3'};
%     genelist = {'Mt1','Ntm','Neat1','Id3','Apoe'};
    genelist = {'Mt1','Sgk1','Fth1','Gstp1','Apoe'};
    c = k;
    hf = figure('color','w','position',[20,20,877,180],'Name',regexprep(clusteruni{c,1},'_','-'));
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
%     eval(['export_fig GABA_Cfl1_Actb_220606/',clusteruni{k,2},'_Cfl1_Actb_',date,'.pdf']);
%     close(hf)

% end
