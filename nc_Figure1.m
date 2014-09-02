function nc_Figure1
% Code for Figure 1 in Yeatman et al. 2014
%
% nc_Figure1
%
% This code will reproduce the top and middle panels of Figure 1.
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications

% Order the fiber groups based on the amount of developmental change
[fgnumsr1, fgnumsmd, fgnumsfa] = nc_SortByGrowth;

% Change to directory with example data
cd(fullfile(nc_Path,'data','exampleSubject'))

% Load segmentation image. This image marks left hemisphere white matter as
% 3 and right hemisphere as 4. It is registered to the subjects dMRI data.
seg = readFileNifti('t1_class_2DTI.nii.gz');

% Convert the segmentation into a binary image where all white matter is
% denoted with a 1 and everything else a zero
seg.data = seg.data==3 | seg.data==4;

% Create a mesh of the cortical surface. The 'boxfilter' argument smooths
% the image before creating a mesh. This leads to a more pleasant
% appearance.
msh = AFQ_meshCreate(seg ,'boxfilter', 5);

% Cut away some of the cortical surface so we can visualize the white
% matter beneath.
plane1 = [-100 0 0;-20 0 0];
plane2 = [0 0 0; 0 0 0];
msh = AFQ_meshCut(msh,plane1);
msh = AFQ_meshCut(msh,plane2);

% Load the fiber groups for this subject. All the cleaned fiber groups
% defined by AFQ were saved in a variable containing the fiber group
% structure and the fiber group names
load exampleFibers.mat

% Render the fibers. Many variables can be set to change the appearance of
% the rendering. The ones here were used in figure 1 of the manuscript
numfibers=250; colors = AFQ_colormap('bgr',24);
lightH = AFQ_RenderFibers(fg(fgnumsr1(1)),'color',colors(1,:),'radius',.5,'numfibers',numfibers);
for ii = 2:length(fgnumsr1)
    AFQ_RenderFibers(fg(fgnumsr1(ii)),'color',colors(ii,:),'radius',.5,'numfibers',numfibers,'newfig',0);
end
axis image
axis off
set(gcf,'color',[0 0 0],'inverthardcopy','off')

% Add the cortical surface to the rendering
p = patch(msh.tr); shading interp

% Rotate the surface to a nice view
camorbit(45,20);

% Set the lighting
camlight(lightH,40,20)

% Save the image
print('-dpng','-r300','Figure1_wmandCortex.png')

% Now save out figures showing different views
camorbit(-45,-20);
camlight(lightH,'right')
delete(p)
print('-dpng','-r300','Figure1_wm1.png')
camorbit(90,0)
print('-dpng','-r300','Figure1_wm2.png')
camorbit(90,0)
camlight(lightH,'right')
print('-dpng','-r300','Figure1_wm3.png')
camorbit(90,0)
print('-dpng','-r300','Figure1_wm4.png')
camorbit(180,90)
camlight(lightH,'infinite')
print('-dpng','-r300','Figure1_wm5.png')

%% Make a legend
figure; hold on
n=0;
for ii = 1:length(fgnumsr1)
    n = n+1;
    if ismember(ii,1:2:100)
        plot(1,n./2,'ks','markerfacecolor',colors(ii,:),'markersize',10);
        text(1.1,n./2,fgNames{fgnumsr1(ii)},'fontname','times');
    else
        plot(2.1,(n-1)./2,'ks','markerfacecolor',colors(ii,:),'markersize',10);
        text(2.2,(n-1)./2,fgNames{fgnumsr1(ii)},'fontname','times');
    end
end
axis([.5 3 -5 20]);
axis off
set(gcf,'color',[1 1 1],'inverthardcopy','off','paperpositionmode','auto')
print('-depsc','-r300','Figure1_legend.eps')