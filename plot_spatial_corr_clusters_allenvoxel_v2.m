tic
clear all
close all
addpath('/data/matlab_functions/')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load all_neurons_cluster_ave.mat %%%    'meangr_mat_all','clusteruni','geneid_all'  geneid
gabaclusters = loadCellFile('ClusterOrder_step3_GABA_FC_24-Oct-2021.txt');
glut1clusters = loadCellFile('ClusterOrder_step3_Glut1_FC_24-Oct-2021.txt');
glut2clusters = loadCellFile('ClusterOrder_step3_Glut2_FC_24-Oct-2021.txt');
tmp = [gabaclusters(:,2);glut1clusters(:,2);glut2clusters(:,2)];
[~,loc] = ismember(tmp,clusteruni);
meangr_mat_all = meangr_mat_all(:,loc);
clusteruni = clusteruni(loc);


data_m_all = meangr_mat_all;

[~,ind_gr_tmp_mark] = ismember(geneid,geneid_all);

ind_gene_selection = ind_gr_tmp_mark;
data_m_all_sorted = data_m_all(ind_gene_selection,:);
data_m_all_sorted_cn = cent_norm(log2(data_m_all+1));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 25 micron volume size
sizegridd25um = [528 320 456];
sizeGrid = [67 41 58];
% ANO = 3-D matrix of annotation labels
fid = fopen('/data/Technion_analysis/allen_atlas/annotation.raw', 'r', 'l' );
ANOGD = fread( fid, prod(sizegridd25um), 'uint32' );
fclose( fid );
ANOGD = reshape(ANOGD,sizegridd25um);
ANOGD = ANOGD(round(linspace(1,sizegridd25um(1),sizeGrid(1))),:,:);
ANOGD = ANOGD(1:sizeGrid(1),:,1:sizegridd25um(3)/2);

annot2RGB = cell2mat(loadCellFile('/data/Technion_analysis/allen_atlas/Allen_ISH_annotation_to_RGB.txt'));
ANOGD_r = uint8(zeros(size(ANOGD)));
ANOGD_g = ANOGD_r;
ANOGD_b = ANOGD_r;
region_ID_uni = unique(ANOGD(:));
annot2RGB = [annot2RGB;[0,round(1*[255,255,255])]];
for i=1:length(region_ID_uni)
    rgbtmp = uint8(annot2RGB(annot2RGB(:,1)==region_ID_uni(i),2:end));
    if ~isempty(rgbtmp)
        idx = ANOGD==region_ID_uni(i);
        ANOGD_r(idx) = rgbtmp(1);
        ANOGD_g(idx) = rgbtmp(2);
        ANOGD_b(idx) = rgbtmp(3);
    end
end

annotation_key = loadCellFile('/data/Technion_analysis/allen_atlas/AllenBrainReference_forAmit.txt');
annotation_key = annotation_key(2:end,:);
regionid = cell2mat(annotation_key(:,1));
regionname = annotation_key(:,2);
% ANOGD_lev6 = uint32(zeros(size(ANOGD)));
ANOGD_lev6 = uint32(ANOGD);
for i=1:length(region_ID_uni)
    i
    if ~isempty(find(cell2mat(annotation_key(:,1))==region_ID_uni(i)))
        levstring = annotation_key{find(cell2mat(annotation_key(:,1))==region_ID_uni(i)),3};
        lev = length( regexp(levstring,'/'))-1;
        if lev>10
            sublev = strsplit(levstring,'/');
            sublev = sublev(2:end-1);
            sublev = uint32(str2double(sublev(10)))
            ANOGD_lev6(ANOGD==region_ID_uni(i)) = sublev;
        end
    end
end
region_ID_uni = unique(ANOGD_lev6(:));
coronal_contours = cell(sizeGrid(1),length(region_ID_uni));
coronal_contours_flip = cell(sizeGrid(1),length(region_ID_uni));
for i=2:length(region_ID_uni)
    i
    for j=1:sizeGrid(1)
        idx = squeeze(ANOGD_lev6(j,:,:))==region_ID_uni(i);
%         lev = cellfun(@length, regexp(annotation_key(find(cell2mat(annotation_key(:,1))==region_ID_uni(i)),3),'/'))-1;
        if sum(idx(:))>0
%             j
            bw = false(sizegridd25um(2),sizegridd25um(3)/2);
            bw(idx) = true;
%             bw = imdilate(bw,strel('disk',3));
            bw1 = fliplr(bw);
            B = bwboundaries(bw);
            B1 = bwboundaries(bw1);
            for jj=1:length(B)
                B{jj} = [smoothn(B{jj}(:,1),5), smoothn(B{jj}(:,2),5)];
%                 B{jj} = [B{jj};B{jj}(1,:)];
                B1{jj} = [smoothn(B1{jj}(:,1),5), smoothn(B1{jj}(:,2),5)];
%                 B1{jj} = [B1{jj};B1{jj}(1,:)];
            end
            coronal_contours{j,i} = B;
            coronal_contours_flip{j,i} = B1;
        end
    end       
end

% sizeGrid = [67 41 58];
% ANOGD = 3-D matrix of grid-level annotation labels 200um volxel
fid = fopen( '/data/Technion_analysis/allen_atlas/gridAnnotation.raw', 'r', 'l' );
ANOGD_200 = fread( fid, prod(sizeGrid), 'uint32' );
fclose( fid );
ANOGD_200 = reshape(ANOGD_200,sizeGrid);


load /data/Technion_analysis/allen_atlas/voxel_200um_data.mat %,'voxel_energy','gene_folder','ind_coro','ind_sagi','sizeGrid'

lefthalf = cat(3, true(67,41,29),false(67,41,29));
lefthalf_lin = find(lefthalf(:));
righthalf_lin = find(~lefthalf(:));
% voxel_energy(ind_coro,lefthalf_lin) = (voxel_energy(ind_coro,lefthalf_lin) + fliplr(voxel_energy(ind_coro,righthalf_lin)) )/2;
voxel_energy = voxel_energy(:,lefthalf_lin);
ANOGD_200 = ANOGD_200(1:67,1:41,1:29);
outvoxels = find(ANOGD_200(:)==0);
invoxels = find(ANOGD_200(:)>0);
voxel_energy(:,outvoxels) = -1;

uniANO = unique(ANOGD_200(:));
uniANOfreq = zeros(size(uniANO));
ANOGDlike = zeros(size(ANOGD_200));
regionind_200 = cell(length(uniANO),1);
for i=1:length(uniANO)
    ind = ANOGD_200(:)==uniANO(i);
    uniANOfreq(i) = sum(ind);
    ANOGDlike(ind) = i-1;
    regionind_200{i} = find(ind);
end

regioniduni = cell2mat(annotation_key(:,1));
region_alllevels_mat = zeros(length(regioniduni),10);
for i=1:length(regioniduni)
    t = annotation_key{i,3};
    t = strsplit(t,'/');
    t = cell2mat(cellfun(@str2double,t(2:end-1),'uniformoutput',0));
    region_alllevels_mat(i,1:length(t)) = t;
end
%%
ANOGD_200_vec = ANOGD_200(:);
region_ANOGD_mat = false(length(ANOGD_200_vec),length(regioniduni));
for i=2:length(uniANO)
    ind = find(ANOGD_200_vec==uniANO(i));
    [I,J] = ind2sub(size(region_alllevels_mat),find(region_alllevels_mat==uniANO(i)));    
    t = region_alllevels_mat(find(sum(region_alllevels_mat==uniANO(i),2)),1:unique(J));
    t(t==0) = [];
    [~,loc] = ismember(t,regioniduni);
%     k = sub2ind(size(region_ANOGD_mat),ind,loc);
    region_ANOGD_mat(ind,loc) = 1;
end



%%  
qualgenes = find(sum(voxel_energy(:,invoxels)>5,2)>30 & mean(voxel_energy(:,invoxels),2)>0.2);

% geneid_coro = gene_folder(ind_coro);
% rcoro = corr_mat(voxel_energy(ind_coro,:)'); 
geneidvalid = gene_folder(qualgenes);

voxel_energy_cn = zeros(size(voxel_energy));
for i=1:length(gene_folder)
    in = voxel_energy(i,:)>0;
    m = mean(log2(voxel_energy(i,in)+1));
    s = std(log2(voxel_energy(i,in)+1));
    voxel_energy_cn(i,in) = (log2(voxel_energy(i,in)+1) - m)/s;
end
   
[valid_genes, ia, ib] = intersect(geneid,geneidvalid);

% r = corr(log2(data_m_sorted(ia,:)+1), voxel_energy_cn(qualgenes(ib),:));
r = corr(data_m_all_sorted_cn(ind_gr_tmp_mark(ia),:), voxel_energy_cn(qualgenes(ib),:));
% r = corr(log2(data_m_sorted(ia,:)+1), voxel_energy(qualgenes(ib),:));
r(:,outvoxels) = 0;

rint = zeros(length(r(:,1)),67*sizegridd25um(2)*sizegridd25um(3)/2);
for i=1:length(r(:,1))
    i
    V = reshape(r(i,:),[67 41 29]);
    [x1,x2,x3] = meshgrid([1:67],linspace(1,41,sizegridd25um(2)),linspace(1,29,sizegridd25um(3)/2));
    Vq = interpn(V,x1,x2,x3,'linear');
    Vq = permute(Vq,[2,1,3]);
    rint(i,:) = Vq(:)';
end


clusterid_all2 = clusteruni;
% for i=1:length(clusterid_all)
%     k = strfind(clusterid_all{i},'-');
%     k2 = strfind(clusterid_all{i},'(');
%     clusterid_all2{i} = [clusterid_all{i}(k(1)+1:k2-1),'-',clusterid_all{i}(1:k(1)-1),clusterid_all{i}(k2:end)];
% end
table_clusters_region = cell(100,2*length(r(:,1)));
zmax = 0;
for i=1:length(r(:,1))
    sizeGrid = [67 41 29];
    x = reshape(r(i,:),sizeGrid);
    [~,xi] = sort(r(i,:),'descend');
    toppixel = sort(ANOGD_200(xi(1:500)));
    [toppixeluni, ia] = unique(toppixel,'first');
    [~, ib] = unique(toppixel,'last');
    [~,loc] = ismember( toppixeluni, uniANO);
    reffreq = uniANOfreq(loc);
    foldenrich = (ib-ia+1)./reffreq;
    foldenrich = foldenrich/sum(foldenrich);
    
    [~,loc] = ismember( toppixeluni, regionid);
    [~,xi] = sort(foldenrich,'descend');
    z = [m2c(foldenrich(xi)), regionname(loc(xi))];
    table_clusters_region(1:length(z(:,1)),(i*2-1):i*2) = z;
    zmax = max([zmax,length(z(:,1))]);    
end
table_clusters_region = table_clusters_region(1:zmax,:);
h = [clusterid_all2';clusterid_all2'];
h = h(:)';
table_clusters_region = [h;repmat({'foldEnrich','regionName'},1,length(r(:,1)));table_clusters_region];
saveCellFile(table_clusters_region,'table_clusters_region_test.txt')
   
%% 

celltype_region_corrsum = zeros(length(r(:,1)),length(regioniduni));
[~,loc] = ismember( uniANO, regionid);
regionnamemat = annotation_key(:,2);%[{'root'};regionname(loc(2:end))];
for i=1:length(r(:,1))
    sizeGrid = [67 41 29];
    x = reshape(r(i,:),sizeGrid);
    for j=1:length(region_ANOGD_mat(1,:))
        celltype_region_corrsum(i,j) = mean(x(region_ANOGD_mat(:,j)));
    end    
end
% celltype_region_corrsum(:,uniANOfreq<5) = 0;
hannah_regions = loadCellFile('ABregions_for_score_13-Feb-2022.txt');
[~,inamyregions] = ismember(hannah_regions,regionnamemat);
% [~,inamyregions] = ismember(hannah_regions,regionname);
% inamyregions = find(contains(lower(regionnamemat),'amy') | contains(lower(regionnamemat),'lateral olfactory'));
tmp = celltype_region_corrsum(:,inamyregions)';
z = linkage(tmp,'average','correlation');
d = squareform(pdist(tmp,'correlation'),'tomatrix');
leaforder1 = optimalleaforder(z,d);%[1:length(inamyregions)];% 
tmp = celltype_region_corrsum(:,inamyregions);
z = linkage(tmp,'average','correlation');
d = squareform(pdist(tmp,'correlation'),'tomatrix');
leaforder2 = [1:length(clusteruni)];%optimalleaforder(z,d);% 

figure('position',[100,100,1500,800],'color','w');
ax1 = axes('Position',[0.13,0.05,0.82,0.6]);
imagesc(celltype_region_corrsum(leaforder2,inamyregions(leaforder1))')
set(gca,'ytick',[1:length(inamyregions)],'yticklabel',regionnamemat(inamyregions(leaforder1))....
    ,'xtick',[1:length(leaforder2)],'xticklabel',(leaforder2),'XTickLabelRotation',45)
cm = [[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']];
colormap((cm))
ax2 = axes('Position',[0.13,0.73,0.82,0.25]);
hden = dendrogram(z,length(leaforder2),'Reorder',leaforder2,'Orientation','top');
eval(['export_fig regions_vs_celltypes_NOTclustered_heatmap_',date,'.pdf']);

table1 = [ [{'clustername';'regionvol'};clusterid_all2], [regionnamemat';m2c([sum(region_ANOGD_mat);celltype_region_corrsum]) ] ];
saveCellFile(table1',['celltype_by_regions_score_',date,'.txt'])
%% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% clusterid_all2 = clusterid_all;
% for i=1:length(clusterid_all)
%     k = strfind(clusterid_all{i},'-');
%     k2 = strfind(clusterid_all{i},'(');
%     clusterid_all2{i} = [clusterid_all{i}(k(1)+1:k2-1),'-',clusterid_all{i}(1:k(1)-1),clusterid_all{i}(k2:end)];
% end


poolobj = parpool(10);
% sizeGrid = [67 41 29];
sizeGrid = [67 sizegridd25um(2) sizegridd25um(3)/2];
parfor i = 1:length(r(:,1))
    i
    clusterid_all2{i}
    
    x = reshape(rint(i,:),sizeGrid);
    x_borders = [0,prctile(x(:),99.9)];
    hf1 = figure('visible','off','position',[0,36,1920,1100],'color','w','name',clusterid_all2{i});
    [ha, pos] = tight_subplot(2, 10, [0.00,0.00], [0.00,0.00], [0.00,0.00]);
    k=0;    
    
    for j=31:40;% j=round(linspace(8,63,25))
        k =k+2;
        %         axes(ha(k));
        set(gcf,'CurrentAxes',ha(k))
        
        tmp = squeeze(x(j,:,:));
        imagesc(fliplr(tmp),x_borders); hold on; %imagesc(fliplr(tmp),[0.05,0.35]); hold on;
        
        cmap = colormap('pink');
        cmap = flipud(cmap);
        %         cmap = [ones(64,1), linspace(1,0,64)', linspace(1,0,64)'];
        set(gcf,'colormap',cmap); freezeColors;
        btmp = coronal_contours_flip(j,:);
        for kkk=1:length(btmp)
            btmp2 = btmp{kkk};
            if ~isempty(btmp2)
                for jjj=1:length(btmp2)
                    plot(btmp2{jjj}(:,2), btmp2{jjj}(:,1), 'w','linewidth',0.5,'color',[0.3*[1,1,1],0.2]);
                end
            end
        end
        axis equal;
        axis off
        axis tight
        set(ha(k),'xlim',[41,190],'ylim',[170,320])
        %         tmp = squeeze(ANOGDlike(j,:,:));
        %         bw = tmp>0;
        %         bw = fliplr(bw);
        %         B = bwboundaries(bw);
        %         for jj=1:length(B)
        %             plot(B{jj}(:,2),B{jj}(:,1),'-r');
        %         end
        %         axes(ha(k-1));
        set(gcf,'CurrentAxes',ha(k-1))
        tmp = cat(3,squeeze(ANOGD_r(j,:,:)), squeeze(ANOGD_g(j,:,:)), squeeze(ANOGD_b(j,:,:)) );
        imagesc(tmp); hold on;
        colormap(colorcube);freezeColors;
        axis equal;
        axis off
        
        btmp = coronal_contours(j,:);
        for kkk=1:length(btmp)
            btmp2 = btmp{kkk};
            if ~isempty(btmp2)
                for jjj=1:length(btmp2)
                    plot(btmp2{jjj}(:,2), btmp2{jjj}(:,1), 'w','linewidth',1,'color',[0.4*[1,1,1],0.4]);
                end
            end
        end
        axis equal;
        axis off
        axis tight
        set(ha(k-1),'xlim',[41,190],'ylim',[170,320])
        text(10,20,[num2str(j)],'color','r','fontsize',12)
        %         if k==6
        %             text(-1,1,clusterid_all2{leaforder(i)},'color','k','fontsize',20,'HorizontalAlignment','center');
        %         end
    end
    axes('position',[0.5,0.9,0.01,0.01])
    title([num2str(i),'-',regexprep(clusterid_all2{(i)},'_','-')],'fontsize',14)
    axis off
    %     eval(['export_fig -depsc /mnt/sanger-data2/10X_adol_mice/L2_3_4_test_build_20171205/spatial_correlation_Dec6_2017/'.....
    %         ,clusterid_all2{leaforder(i)},'_',date,'.pdf']);
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png(['/data/Technion_analysis/Amygdala/Amy_FC_final/spatial_correlation_211010/',num2str(i),'_',....
        clusterid_all2{(i)},'_',date,''],gcf,500)
    close(hf1)
    
end
delete(poolobj)

%%


poolobj = parpool(10);
% sizeGrid = [67 41 29];
sizeGrid = [67 sizegridd25um(2) sizegridd25um(3)/2];
parfor i = 1:length(r(:,1))
    i
    clusterid_all2{i}
    
    x = reshape(rint(i,:),sizeGrid);
    x_borders = [0,prctile(x(:),99.9)];
    hf1 = figure('visible','on','position',[0,36,1920,500],'color','w','name',clusterid_all2{i});
    [ha, pos] = tight_subplot(1, 10, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
    k=0;    
    
    for j=31:40;% j=round(linspace(8,63,25))
        k =k+1;
        %         axes(ha(k));
        set(gcf,'CurrentAxes',ha(k))
        
        tmp = squeeze(x(j,:,:));
        imagesc(fliplr(tmp),x_borders); hold on; %imagesc(fliplr(tmp),[0.05,0.35]); hold on;
        
        cmap = colormap('pink');
        cmap = flipud(cmap);
        %         cmap = [ones(64,1), linspace(1,0,64)', linspace(1,0,64)'];
        set(gcf,'colormap',cmap); freezeColors;
        btmp = coronal_contours_flip(j,:);
        for kkk=1:length(btmp)
            btmp2 = btmp{kkk};
            if ~isempty(btmp2)
                for jjj=1:length(btmp2)
                    plot(btmp2{jjj}(:,2), btmp2{jjj}(:,1), 'w','linewidth',0.5,'color',[0.3*[1,1,1],0.2]);
                end
            end
        end
        axis equal;
        axis off
        axis tight
        set(ha(k),'xlim',[41,190],'ylim',[170,320])
%         %         tmp = squeeze(ANOGDlike(j,:,:));
%         %         bw = tmp>0;
%         %         bw = fliplr(bw);
%         %         B = bwboundaries(bw);
%         %         for jj=1:length(B)
%         %             plot(B{jj}(:,2),B{jj}(:,1),'-r');
%         %         end
%         %         axes(ha(k-1));
%         set(gcf,'CurrentAxes',ha(k-1))
%         tmp = cat(3,squeeze(ANOGD_r(j,:,:)), squeeze(ANOGD_g(j,:,:)), squeeze(ANOGD_b(j,:,:)) );
%         imagesc(tmp); hold on;
%         colormap(colorcube);freezeColors;
%         axis equal;
%         axis off
%         
%         btmp = coronal_contours(j,:);
%         for kkk=1:length(btmp)
%             btmp2 = btmp{kkk};
%             if ~isempty(btmp2)
%                 for jjj=1:length(btmp2)
%                     plot(btmp2{jjj}(:,2), btmp2{jjj}(:,1), 'w','linewidth',1,'color',[0.4*[1,1,1],0.4]);
%                 end
%             end
%         end
%         axis equal;
%         axis off
%         axis tight
%         set(gca,'xdir','reverse')
%         set(ha(k),'xlim',[41,190],'ylim',[170,320])
% %         text(100,-20,[num2str(j)],'color','r','fontsize',12)
%         %         if k==6
%         %             text(-1,1,clusterid_all2{leaforder(i)},'color','k','fontsize',20,'HorizontalAlignment','center');
%         %         end
%     end
%     axes('position',[0.5,0.9,0.01,0.01])
%     title([num2str(i),'-',regexprep(clusterid_all2{(i)},'_','-')],'fontsize',14)
%     axis off
    %     eval(['export_fig -depsc /mnt/sanger-data2/10X_adol_mice/L2_3_4_test_build_20171205/spatial_correlation_Dec6_2017/'.....
    %         ,clusterid_all2{leaforder(i)},'_',date,'.pdf']);
    set(gcf,'Renderer','OpenGL');
    set(gcf,'PaperPositionMode','auto');
    save2png(['/data/Technion_analysis/Amygdala/Amy_FC_final/spatial_correlation_211107/',num2str(i),'_',....
        clusterid_all2{(i)},'_',date,''],gcf,500)
    close(hf1)
    
    end
end
delete(poolobj)
%%