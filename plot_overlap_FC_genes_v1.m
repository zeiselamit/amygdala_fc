tic
clear all
close all
addpath('/data/matlab_functions/')

%%
load pbino_glut1_14-Nov-2021.mat %'pbino_mat','geneid','clusteruni'
pbino_mat_g1 =  pbino_mat;
clusteruni_g1 = clusteruni;
nc_g1 = length(clusteruni_g1);
hgenes_g1 = loadCellFile('hitgenes_glut1_14-Nov-2021.txt');
load p90mat_glut1_16-Nov-2021.mat %'
p90mat_glut1 = p90_mat;
load pbino_glut2_14-Nov-2021.mat %'pbino_mat','geneid','clusteruni'
pbino_mat_g2 =  pbino_mat;
clusteruni_g2 = clusteruni;
nc_g2 = length(clusteruni_g2);
load p90mat_glut2_16-Nov-2021.mat %'
p90mat_glut2 = p90_mat;
hgenes_g2 = loadCellFile('hitgenes_glut2_14-Nov-2021.txt');
load pbino_gaba_14-Nov-2021.mat %'pbino_mat','geneid','clusteruni'
pbino_mat_gaba =  pbino_mat;
clusteruni_gaba = clusteruni;
nc_gaba = length(clusteruni_gaba);
load p90mat_gaba_16-Nov-2021.mat %'
p90mat_gaba = p90_mat;
hgenes_gaba = loadCellFile('hitgenes_gaba_14-Nov-2021.txt');



c1 = length(hgenes_g1(:,1));
c2 = length(hgenes_g2(:,1));
c3 = length(hgenes_gaba(:,1));
c12 = length(intersect(hgenes_g1(:,1),hgenes_g2(:,1)));
c13 = length(intersect(hgenes_g1(:,1),hgenes_gaba(:,1)));
c23 = length(intersect(hgenes_g2(:,1),hgenes_gaba(:,1)));
c123 = length(intersect(hgenes_gaba(:,1),intersect(hgenes_g1(:,1),hgenes_g2(:,1))));

figure;
venn([c1,c2,c3],[c12,c13,c23,c123])
axis tight
axis equal
axis off
legend({'Glut1','Glut2','GABA'})
%%
glist = unique([hgenes_g1(:,1);hgenes_g2(:,1);hgenes_gaba(:,1)]);

[~,inhit] = ismember(glist,geneid);
m = -log10([pbino_mat_g1(inhit,:),pbino_mat_g2(inhit,:),pbino_mat_gaba(inhit,:)]);
z = linkage(m,'ward','correlation');
d = corr_mat(m');
leaforder = optimalleaforder(z,squareform(1-d,'tovector'));
tp_order = reshape([1:5*(nc_g1+nc_g2+nc_gaba)],5,(nc_g1+nc_g2+nc_gaba))';
tp_order = tp_order(:);
figure('position',[200,200,1200,1000],'color','w');
set(gca,'Position',[0.1,0.03,0.85,0.9])
imagesc((m(leaforder,tp_order)),[2,3]); hold on;
colormap('summer');

freezeColors(gca);
yl = get(gca,'ylim');
for i=1:5
    plot(i*(nc_g1+nc_g2+nc_gaba)*[1,1]+0.5,yl,'-k');
    plot(i*(nc_g1+nc_g2+nc_gaba)*[1,1]+0.5-nc_gaba,yl,'--k');
    plot(i*(nc_g1+nc_g2+nc_gaba)*[1,1]+0.5-nc_gaba-nc_g2,yl,'--k');
end
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',glist(leaforder,1),'fontsize',10,'xtick'...
    ,[(nc_g1+nc_g2+nc_gaba)/2:(nc_g1+nc_g2+nc_gaba):5*(nc_g1+nc_g2+nc_gaba)],'XTickLabel',{'0','2','8','24','recall'})
n = (nc_g1+nc_g2+nc_gaba);
x = [mean(m(:,tp_order(1:n)),2),mean(m(:,tp_order(n+1:2*n)),2),....
    mean(m(:,tp_order(2*n+1:3*n)),2),mean(m(:,tp_order(3*n+1:4*n)),2),mean(m(:,tp_order(4*n+1:5*n)),2)];
[~,imax] = max(x,[],2);
figure; 
hist(imax,5)
%%

glist = unique([hgenes_g1(:,1);hgenes_g2(:,1);hgenes_gaba(:,1)]);

[~,inhit] = ismember(glist,geneid);
m = [p90mat_glut1(inhit,:),p90mat_glut2(inhit,:),p90mat_gaba(inhit,:)];
for i=1:length(m(1,:))/5
    ind = (i-1)*5+1:i*5;
    m(:,ind) = m(:,ind) - repmat(m(:,ind(1)),1,5);
end
z = linkage(m,'ward','correlation');
d = corr_mat(m');
leaforder = optimalleaforder(z,squareform(1-d,'tovector'));
tp_order = reshape([1:5*(nc_g1+nc_g2+nc_gaba)],5,(nc_g1+nc_g2+nc_gaba))';
tp_order = tp_order(:);
figure('position',[200,200,1200,1000],'color','w');
set(gca,'Position',[0.1,0.03,0.85,0.9])
mtmp = cent_norm(m(leaforder,tp_order));
imagesc(mtmp,[prctile(mtmp(:),1),prctile(mtmp(:),99)]); hold on;
colormap('summer');

freezeColors(gca);
yl = get(gca,'ylim');
for i=1:5
    plot(i*(nc_g1+nc_g2+nc_gaba)*[1,1]+0.5,yl,'-k');
    plot(i*(nc_g1+nc_g2+nc_gaba)*[1,1]+0.5-nc_gaba,yl,'--k');
    plot(i*(nc_g1+nc_g2+nc_gaba)*[1,1]+0.5-nc_gaba-nc_g2,yl,'--k');
end
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',glist(leaforder,1),'fontsize',10,'xtick'...
    ,[(nc_g1+nc_g2+nc_gaba)/2:(nc_g1+nc_g2+nc_gaba):5*(nc_g1+nc_g2+nc_gaba)],'XTickLabel',{'0','2','8','24','recall'})










