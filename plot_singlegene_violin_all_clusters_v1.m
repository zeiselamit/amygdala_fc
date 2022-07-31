tic
clear all
close all
addpath('/data/matlab_functions/')


load geneid
load(['amyFC_NONneurons_sc_celltypes_muhammad_25-Nov-2021.mat'])%%%,'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time')
data_nn = data;
cellid_nn = cellid;
c_nn = c;

load(['amyFC_neurons_sc_celltypes_muhammad_25-Nov-2021.mat'])%%%,'data','cellid','sample','cellid_clusters','c','amy_flag','batch_flag','fc_time')
data = [data,data_nn];
data = normalize(data,'norm');
data = round(data./repmat(sum(data),length(data(:,1)),1)*20e3);
cellid = [cellid;cellid_nn];
c = [c;c_nn];


gabahex = loadCellFile_turbo('GABA-hexa.txt',1);
gabahex = gabahex(2:end,:);
glut1hex = loadCellFile_turbo('Glut1-hexa.txt',1);
glut1hex = glut1hex(2:end,:);
glut2hex = loadCellFile_turbo('Glut2-hexa.txt',1);
glut2hex = glut2hex(2:end,:);
nnhex = loadCellFile_turbo('non-neuronal-hexa.txt',1);
nnhex = nnhex(2:end,:);

gabaorder = loadCellFile('ClusterOrder_step3_GABA_FC_04-Nov-2021.txt');
glut1order = loadCellFile('ClusterOrder_step3_Glut1_FC_12-Dec-2021.txt');
glut2order = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');

clusteruni = [gabaorder(:,2);glut1order(:,2);glut2order(:,2)];

cluster_info = loadCellFile('Clusters_info_hannah.txt');
cluster_info(1,:) = [];
clusteruni = [clusteruni,cellfun(@(x,y) [x,'-',y] , cluster_info(:,4),cluster_info(:,5),'UniformOutput',0)];
for i=1:length(clusteruni)
    in = find(strcmpi(clusteruni{i,1},c));
    c(in) = repmat(clusteruni(i,2),length(in),1);
end

cuni = [clusteruni(:,2);....
    {'Astro';'Astro_agt';'Astro_SC';'Epend';'OPC_cycling';'OPC';'COP';'OL';'EC';'Peri';'VSM';'VLMC';'microglia';'pvm'}];
r1 = cuni(1:56);
r2 = cuni(57:88);
r3 = cuni(89:130);
r4 = {'Astro';'Astro_agt';'Astro_SC';'Epend';'OPC_cycling';'OPC';'COP';'OL';'EC';'Peri';'VSM';'VLMC';'Microglia';'PVM'};
nn_colors = distinguishable_colors(length(r4));
in = find(sum(data>0,2)>50 );
logflag = 0;
%%
list = geneid(in);
% list = {'Rspo2','Ptk2b','Meis2','Wfs1','Maf','Zeb2','Satb1','Tshz1','Drd1','Pdyn','Drd2','Adora2a','Foxp2','Nfib','Lhx6','Resp18','Ppp1r1b'};
% [~,in] = ismember(list,geneid);
% c = cuni;
fname = dir('/data/Technion_analysis/Amygdala/Amy_FC_final/singlegene_violoin_all_clusters_220717/*png');
fsize = [fname.bytes];
fname = {fname(:).name}';
fname = regexprep(fname,'_alltypes_violinSingle.png','');
list = setdiff(list,fname);
[~,in] = ismember(list,geneid);

poolobj = parpool(8);
parfor jjj=1:length(list)
    jjj
    hf1 = figure;
    set(gcf,'color','w','position',[20,20,1000,900],'Visible','off');    
    g = in(jjj);%find(strcmpi(gn,geneid));
    gn = list{jjj}
    p99 = prctile(data(g,data(g,:)>0),98);
    a1 = axes('Position',[0.05,0.15,0.9,0.1]);
%     p90 = max(ceil(prctile(data(g,:),90)),1);
%     t_ed = zeros(length(T_cells_tmp_uni),1);
    
%     t_75 = zeros(length(T_cells_tmp_uni),1);
t_av = zeros(length(r1),1);
    for k=1:length(r1)
        k;
        gr2 = find(strcmpi(r1{k},c));
        y = (data(g,gr2));
        if length(y)>10
            [f,xi] = ksdensity(y);
%             [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
            fi = interp1(xi,f,y);
            fi = fi/max(fi);            
            plot(k + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.',...
                'color',['#',gabahex{k,3}]); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+y'...
                ,'.','color',['#',gabahex{k,3}]); hold on;
        end
%         t_ed(k) = median(data_orig_all_sorted(g,gr2));
        t_av(k) = mean(y);
        p99 = max([p99,prctile(y,98)]);
%         t_75(k) = prctile(data_orig_all_sorted(g,gr2),75);
    end
    plot([1:length(t_av)],t_av,'sk','MarkerEdgeColor','none','MarkerFaceColor','k','markersize',2)
    axis tight
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
    xt = regexprep(r1,'GABA-','');
    xt = regexprep(xt,'_','-');
    ylabel('UMI');
    set(a1,'xtick',[1:length(r1)],'xticklabel',xt,'xticklabelrotation',45....
        ,'ylim',[-0.5,p99+1],'fontsize',6);
    box off
    text(length(r1)/2,p99,'GABA')
    
    a2 = axes('Position',[0.05,0.35,0.9,0.1]);
    t_av = zeros(length(r2),1);
    for k=1:length(r2)
        k;
        gr2 = find(strcmpi(r2{k},c));
        y = (data(g,gr2));
        if length(y)>10
%             [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
            [f,xi] = ksdensity(y);
            fi = interp1(xi,f,y);
            fi = fi/max(fi);
            plot(k + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.',...
                'color',['#',glut1hex{k,3}]); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+y'...
                ,'.','color',['#',glut1hex{k,3}]); hold on;
        end
%         t_ed(k) = median(data_orig_all_sorted(g,gr2));
        t_av(k) = mean(y);
        p99 = max([p99,prctile(y,98)]);
%         t_75(k) = prctile(data_orig_all_sorted(g,gr2),75);
    end
    plot([1:length(t_av)],t_av,'sk','MarkerEdgeColor','none','MarkerFaceColor','k','markersize',2)
    axis tight
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
    xt = regexprep(r2,'VGLUT1-','');
    xt = regexprep(xt,'_','-');
    ylabel('UMI');
    set(a2,'xtick',[1:length(r2)],'xticklabel',xt,'xticklabelrotation',45....
        ,'ylim',[-0.5,p99+1],'fontsize',6);
    box off
    text(length(r2)/2,p99,'vGlut1')
    
    a3 = axes('Position',[0.05,0.55,0.9,0.1]);
    t_av = zeros(length(r3),1);
    for k=1:length(r3)
        k;
        gr2 = find(strcmpi(r3{k},c));
        y = (data(g,gr2));
        if length(y)>10
%             [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
            [f,xi] = ksdensity(y);
            fi = interp1(xi,f,y);
            fi = fi/max(fi);
            plot(k + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.',...
                'color',['#',glut2hex{k,3}]); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+y'...
                ,'.','color',['#',glut2hex{k,3}]); hold on;
        end
%         t_ed(k) = median(data_orig_all_sorted(g,gr2));
        t_av(k) = mean(y);
        p99 = max([p99,prctile(y,98)]);
%         t_75(k) = prctile(data_orig_all_sorted(g,gr2),75);
    end
    plot([1:length(t_av)],t_av,'sk','MarkerEdgeColor','none','MarkerFaceColor','k','markersize',2)
    axis tight
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
    xt = regexprep(r3,'VGLUT2-','');
    xt = regexprep(xt,'_','-');
    ylabel('UMI');
    set(a3,'xtick',[1:length(r3)],'xticklabel',xt,'xticklabelrotation',45....
        ,'ylim',[-0.5,p99+1],'fontsize',6);
    box off
    text(length(r3)/2,p99,'vGlut2')
    
    a4 = axes('Position',[0.05,0.75,0.9,0.1]);    
    t_av = zeros(length(r4),1);
    for k=1:length(r4)
        k;
        gr2 = find(strcmpi(r4{k},c));
        y = (data(g,gr2));
        if length(y)>10
%             [f,xi] = ksdensity(y,'support',[0-0.1,0.1+max(y)]);
            [f,xi] = ksdensity(y);
            fi = interp1(xi,f,y);
            fi = fi/max(fi);
            plot(k + fi'.*(0.9*rand(length(gr2),1)-0.45), 0.5*rand(length(gr2),1)-0.1+y','.',...
                'color',['#',nnhex{k,3}]); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+y'...
                ,'.','color',['#',nnhex{k,3}]); hold on;
        end
%         t_ed(k) = median(data_orig_all_sorted(g,gr2));
        t_av(k) = mean(y);
        p99 = max([p99,prctile(y,98)]);
%         t_75(k) = prctile(data_orig_all_sorted(g,gr2),75);
    end
    plot([1:length(t_av)],t_av,'sk','MarkerEdgeColor','none','MarkerFaceColor','k','markersize',2)
    axis tight
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
    xt = regexprep(r4,'_','-');
    ylabel('UMI');
    set(a4,'xtick',[1:length(r4)],'xticklabel',xt,'xticklabelrotation',45....
        ,'ylim',[-0.5,p99+1],'fontsize',6);
    box off
    text(length(r4)/2,p99,'Non-Neurons')
    linkaxes([a1,a2,a3,a4],'y');
    title(gn,'fontsize',16)
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png(['singlegene_violoin_all_clusters_220717/',gn,'_alltypes_violinSingle'],gcf,700)
    close(hf1)
end
delete(poolobj)
%%