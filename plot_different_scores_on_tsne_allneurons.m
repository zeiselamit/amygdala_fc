%%
tsneall = loadCellFile('celltype_tsne_all_xy_size_colors_30-Jun-2022.txt');%'oldname','newname','hexcol','size','x','y'};
% tsneall = tsneall(2:end,:);
xy = cell2mat(tsneall(:,5:6));

clusterIEGscore = loadCellFile('GABA_IEGfrac_score_p90_grsz20.txt');
clusterIEGscore = [clusterIEGscore;loadCellFile('Glut12_IEGfrac_score_p90_grsz20.txt')];


[~,loc] = ismember(tsneall(:,2),clusterIEGscore(:,1));
iegscore = cell2mat(clusterIEGscore(loc,2));
mingrsz = cell2mat(clusterIEGscore(loc,4));
% iegscore(mingrsz<20) = -1;
% iegscore(mingrsz<10) = -2;
markergene = iegscore';
markergene(markergene>1) = 1;
cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);

% markergene_color = [interp1(linspace(min(markergene),max(markergene),length(cm)),cm(:,1),markergene'),...
%         interp1(linspace(min(markergene),max(markergene),length(cm)),cm(:,2),markergene')...
%         ,interp1(linspace(min(markergene),max(markergene),length(cm)),cm(:,3),markergene')];
markergene_color = [interp1(linspace(min(markergene),1,length(cm)),cm(:,1),markergene'),...
        interp1(linspace(min(markergene),1,length(cm)),cm(:,2),markergene')...
        ,interp1(linspace(min(markergene),1,length(cm)),cm(:,3),markergene')];

figure;
set(gcf,'color','w','position',[20,20,900,800])
scatter(xy(:,1),xy(:,2),cell2mat(tsneall(:,4)),markergene_color,'filled','o','MarkerEdgeColor','k'); hold on;
axis tight;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
rr = 0.05*diff(xl);
for i=1:length(tsneall)
    if i>=0
        in = i;
        th = atan(xy(in,2)/xy(in,1)); %(rand*360)*pi/180;
        xo = rr*cos(th); yo = rr*sin(th);
        plot([xy(in,1),xy(in,1)+xo] , ...
            [xy(in,2),xy(in,2)+yo],'-','color',0.7*[1,1,1]);
        ht = text(xy(in,1)+xo,xy(in,2)+yo,regexprep(tsneall{i,2},'_','-'));
        set(ht,'fontsize',8,'FontName','arial');%,'BackgroundColor',0.8*[1,1,1]
    end
end

axis tight;
axis equal
axis off
eval(['export_fig tsne_allneurons_IEGscore_',date,'.pdf']);

%%
tsneall = loadCellFile('celltype_tsne_all_xy_size_colors_30-Jun-2022.txt');%'oldname','newname','hexcol','size','x','y'};
% tsneall = tsneall(2:end,:);
xy = cell2mat(tsneall(:,5:6));

gabacorrscore = loadCellFile('corrscore_GABA_12-Jul-2022.txt');
glutcorrscore = [loadCellFile('corrscore_Glut12_12-Jul-2022.txt')];
clusterIEGscore = [gabacorrscore(2:end,:);glutcorrscore(2:end,:)];

cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);
figure('position',[400,100,1050,350],'color','w');
[ha, pos] = tight_subplot(1, 8, [0.01,0.01], [0.03,0.08], [0.05,0.01]);
dotsize = 0.2*cell2mat(tsneall(:,4));
[~,xidotsize] = sort(dotsize,'descend');
for i=1:8
    iegscore = cell2mat(clusterIEGscore(:,i+1));
    markergene = iegscore';
    maxscore = prctile(markergene,99);
    markergene(markergene>maxscore) = maxscore;
    markergene_color = [interp1(linspace(min(markergene),maxscore,length(cm)),cm(:,1),markergene'),...
        interp1(linspace(min(markergene),maxscore,length(cm)),cm(:,2),markergene')...
        ,interp1(linspace(min(markergene),maxscore,length(cm)),cm(:,3),markergene')];
    markergene_color(iegscore==0,:) = repmat([0.2,0.2,0.2],sum(iegscore==0),1);

    axes(ha(i))
    
    scatter(xy(xidotsize,1),xy(xidotsize,2),dotsize(xidotsize),markergene_color(xidotsize,:),'filled','o','MarkerEdgeColor','k'); hold on;
    axis tight;
    axis off
    axis equal
    title(gabacorrscore{1,i+1},'interpreter','none')
%     pause
end

eval(['export_fig tsne_allneurons_CorrScores_',date,'.pdf']);