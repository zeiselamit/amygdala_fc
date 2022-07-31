tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;


load afterloading_QC3000_FC_06-Oct-2021
IEG_markers = loadCellFile('IEG_cluster.txt');
[~,in_ieg] = ismember(IEG_markers,geneid);

tot_mol = sum(data);
tot_mol(tot_mol>3e4) = 3e4;
tot_genes = sum(data>0);

cellidclsuter_gaba = loadCellFile('cellidCluster_step3_Glut2_FC_24-Oct-2021.txt');
% cellidclsuter_gaba = [cellidclsuter_gaba, cellfun(@(x,y) ['gaba_',num2str(x),'_',num2str(y)],cellidclsuter_gaba(:,2),cellidclsuter_gaba(:,3),'UniformOutput',0)]; 

gaba_clustername = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');

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
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,4),clusteruni{i}));
    T_cells_tmp(ind) = i;
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp);
fc_time_uni = unique(fc_time);
% % % % % % % % % % % % % % % % % % % 
%%
gn = 'Prkaca';
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
p99 = prctile(data(g,:),99.99);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    c=k;
    t_ed = zeros(5,1);
    t_av = zeros(5,1);
    t_75 = zeros(5,1);
    p90 = prctile(data(g,T_cells_tmp(:,b)==c),90);%max(ceil(prctile(data(g,T_cells_tmp(:,b)==c),90)),1);
    for i=[1:5]
        gr2 = find(T_cells_tmp(:,b)==c & fc_time==fc_time_uni(i));%find(fc_time_sorted==fc_time_uni(i));%
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
            plot(i + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.'); hold on;
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
    text([1:5],p99*0.7*ones(1,5),cellfun(@(x) num2str(x,2), m2c(-log10(pp90)),'UniformOutput',0),'fontsize',6);
    title(regexprep(clusteruni{c},'_','-'));%[gn,',c=',num2str(c)])
end
if 1==1
    %     savefig(gcf,['clustermarkers_violin_final_step3_FC_GABA_',date,'.fig'])
    eval(['export_fig ',gn,'_violin_FC_Glut2_',date,'.pdf']);    
end
%%
in = find(sum(data>0,2)>100 & sum(data>0,2)<length(data(1,:))*0.9);
k = 0;
pbino_mat = ones(length(geneid),5*length(T_cells_tmp_uni));
for c1=1:length(unique(T_cells_tmp(:,b)))
    c1
    fc_time_uni = [0,2,8,24,28];
    p90gene = ((prctile(data(:,T_cells_tmp(:,b)==c1),90,2)));
    for i=[1,2,3,4,5]
        k = k+1
        gr2 = find(T_cells_tmp(:,b)==c1 & fc_time==fc_time_uni(i));
        x90 = sum(data(:,gr2)>repmat(p90gene,1,length(gr2)),2);
        x90(p90gene==0) = nan;
        pbino = binocdf(x90,length(gr2),0.1,'upper');
        pbino(isnan(pbino)) = 1;
        pbino(pbino<1e-300) = 1e-300;
        %         [~,ptt] = ttest2(data(in,gr1)',data(in,gr2)','tail','both');
        pbino_mat(:,k) = pbino;
    end
end
x = sum(pbino_mat(in,:)<1e-2,2);
[~,xi] = sort(x,'descend');
g = [geneid(in(xi(1:100))),m2c(x(xi(1:100)))]
m = -log10(pbino_mat(in(xi(1:100)),:));
m_ieg = -log10(pbino_mat(in_ieg,:));

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
set(gcf,'color','w','position',[20,20,1200,1000],'Name',gn);
[ha, pos] = tight_subplot(8, ceil(length(unique(T_cells_tmp(:,1)))/8), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for k=1:length(T_cells_tmp_uni)
    axes(ha(k))
    bar(msum(k,:));
    title(regexprep(clusteruni(k),'_','-'))
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
    eval(['export_fig top100_mostdiff_FC_Glut2_',date,'.pdf']);    
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


% for i=1:length(T_cells_tmp_uni)-1
%     plot(i*5*[1,1]+0.5,yl,'-k');
% end
% set(gca,'ytick',[1:length(leaforder)],'YTickLabel',g(leaforder),'fontsize',10,'xtick'...
%     ,[3:6:5*length(T_cells_tmp_uni)],'XTickLabel',regexprep(clusteruni,'_','-'),'XTickLabelRotation',45)

%%

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
    eval(['export_fig markers1_violin_final_step3_FC_Glut2_',date,'.pdf']);    
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
    eval(['export_fig markers2_violin_final_step3_FC_Glut2_',date,'.pdf']);
end


%%







