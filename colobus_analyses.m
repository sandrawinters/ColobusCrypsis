% analyzes colobus images using QCPA in MICA 
%
% Sandra Winters <sandra.winters@helsinki.fi>
% updated 22 Nov 2021

%note: requires cone-catch models to be pre-generated in MICA (i.e., 'Generate Cone Mapping Model from Spectral Sensitivities')

%% set up files
cd '~/Dropbox/Color_scripts' % SET TO WORKING DIRECTORY

baseFolder = '~/Dropbox/Colobus'; % SET TO PROJECT FOLDER

if ispc
    foldDelim = '\';
else
    foldDelim = '/';
end

if isfolder(fullfile(baseFolder,'mica_results'))==0
    mkdir(fullfile(baseFolder,'mica_results'))
end
if isfolder(fullfile(baseFolder,'mica_results','rescale_values'))==0
    mkdir(fullfile(baseFolder,'mica_results','rescale_values'))
end

micaDir = fullfile(baseFolder,'mica_files');
roiDir = fullfile(micaDir,'ROIs');

visSystems = {'Cat','Chimp','Raptor'}; 
illuminants = {'D65','CivilTwilight'}; 
distances = [20,100,200,400]; %unit = monkey length, ~0.5 m

roiNames = {'background','body','head','tail','whole'};

%% start MIJ
javaaddpath([matlabroot foldDelim 'java' foldDelim 'mij.jar'])
javaaddpath([matlabroot foldDelim 'java' foldDelim 'ij.jar'])

MIJ.start([matlabroot foldDelim 'java' foldDelim])
IJ = ij.IJ();

%% process images
ims = dir([micaDir foldDelim '*.mspec']);
    
%loop through distances (making this first loop to break things up more for start/restart...)
for dist = distances %2^[1:8] %2, 4, 8, 16, 32, 64, 128, 256 viewing distance (unit = monkey length, ~0.5 m)
    
    %loop through images
    for i = 1:length(files) 
        
        %get image info
        mspecName = strrep(ims(i).name,'.mspec','');

        if isempty(dir([micaDir foldDelim mspecName '*.CR2']))==0
            imName = dir([micaDir foldDelim mspecName '*.CR2']);
        elseif isempty(dir([micaDir foldDelim mspecName '*.NEF']))==0
            imName = dir([micaDir foldDelim mspecName '*.NEF']);
        else
            error(['image not found: ' ims(i).name])
        end
        imName = imName(1).name(6:end);

        if strcmp(imName(1:4),'2008') || strcmp(imName(1:5),'Ghana') %Fernando Canon 40D
            camera = 'Canon 40D'; 
        elseif strcmp(imName(1:4),'2017') %Sandra Canon 5D mkIII
            camera = 'Canon 5D MKIII Tamron 150 to 600mm'; 
        elseif strcmp(imName(1:7),'MONKCBW') %Suzi Canon 1D mkII or 1D X
            camera = 'Canon 1D MKIII';  %similar model
        elseif strcmp(imName(1:10),'EvaWikberg') %Eva Nikon D80
            camera = 'Nikon D80 Nikkor AF 60mm'; 
        else
            error(['unknown camera: ' ims(i).name])
        end

        %loop through illuminants
        for j = 1:length(illuminants)
            illum = illuminants{j};

            %loop through visual systems
            for k = 1:length(visSystems)
                vs = visSystems{k};

                %load mspec image & rois
                disp([mspecName ' ' vs ' ' illum ' ' num2str(dist) ':'])
                MIJ.run(' Load Multispectral Image',['select=' micaDir foldDelim ims(i).name ' image=[Non-linear Colour Image]']);
                IJ.open([roiDir foldDelim mspecName '.zip']);

                %convert to cone catch
                MIJ.run('Convert to Cone Catch',['model=[' camera ' D65 to ' vs ' ' illum '] desaturate desaturation=0.010 remove replace=0.001']);

                if strcmp(vs,'Cat')
                    MIJ.run('Create Luminance Channel','lw'); %cat luminance based on lw cone catch
                elseif strcmp(vs,'Chimp')
                    MIJ.run('Create Luminance Channel','lw mw'); %chimp luminance based on mean of lw & mw cone catches
                else %raptor
                    MIJ.run("Next Slice [>]");
                    MIJ.run("Next Slice [>]");
                    MIJ.run("Next Slice [>]");
                    MIJ.run("Set Label...", "label=lum"); %raptor luminance based on double cone (already in cone catch image; just renaming for consistency)
                end

                %measure ROIs
                MIJ.run("Measure ROIs");
                IJ.saveAs('results',[baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_roiResults.csv']);
                MIJ.run('Clear Results')
                IJ.selectWindow("Results"); 
                IJ.run("Close");

                %get cone catch image
                im_vm = MIJ.getCurrentImage();
                saveImage(im2uint16(im_vm), ... 
                          vs, ... 
                          [baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_1-VM'])

                %set visual system data (for weber fractions, using v = 0.05 for all; order lw, mw, sw)
                if strcmp(vs,'Cat')
                    MRA = 1/10; %Wassle 1971 doi.org/10.1016/0042-6989(71)90219-7
                    weberStr = 'weber_1=0.05 weber_2=0.1225'; %based on receptor density (s:l) 1:6; Loop et al. 1987 doi.org/10.1113/jphysiol.1987.sp016383, Linberg et al. 2001 doi.org/10.1002/1096-9861(20010212)430:3%3C343::AID-CNE1035%3E3.0.CO;2-U
                    weberLum = '0.07'; %Chiao et al. 2000 doi.org/10.1016/S0042-6989(00)00156-5
                elseif strcmp(vs,'Chimp')
                    MRA = 1/64; %Spence 1934 doi.org/10.1037/h0075291
                    weberStr = 'weber_1=0.05 weber_2=0.05 weber_3=0.2'; %based on receptor density (s:m:l) 1:16:16; Knoblauch et al. 2006 doi.org/10.1017/S0952523806233157
                    weberLum = '0.08'; %Osorio & Vorobyev 1996; Osorio et al. 2004
                else %raptor
                    MRA = 1/140; %Reymond 1985 doi.org/10.1016/0042-6989(85)90226-3 (for eagle Aquila audax)
                    weberStr = 'weber_1=0.05 weber_2=0.0707 weber_3=0.0707'; %based on receptor density (uv:s:m:l) 1:2:2:4 Lind et al. 2013 doi.org/10.1242/jeb.082834 --> (s:m:l) 1:2:2 
                    weberLum = '0.1'; %Lind et al. 2013 doi.org/10.1242/jeb.082834 (as in Potier et al. 2018 doi.org/10.1098/rspb.2018.1036)
                end

                %get scale bar length
                roiNames = unzip([micaDir foldDelim 'ROIs' foldDelim mspecName '.zip']);
                scale = roiNames{cellfun(@isempty,regexp(roiNames,'Scale*'))==0};
                scale = split(scale,':');
                scale = str2double(scale{2});

                %generate acuity image
                pxMRA = 3;
                w = size(im_vm,2)/scale; %real width of image in monkey units
                
                [im_av,rescale] = GenerateAcuityImage(im_vm,w,dist,MRA,pxMRA); 
                
                save([baseFolder foldDelim 'mica_results' foldDelim 'rescale_values' foldDelim mspecName vs '_' illum '_' num2str(dist) '_rescale.mat'], ... 
                      'rescale')
                saveImage(im2uint16(im_av), ... 
                          vs, ... 
                          [baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_2-AV'])
                      
                %
                if rescale>1
                    pxMRA = floor(pxMRA/rescale);
                    if pxMRA<1
                        pxMRA = 1;
                    end
                    
                    [im_av,rescale] = GenerateAcuityImage(im_vm,w,dist,MRA,pxMRA); 
                end
                
                %rescale cone catch image (to get ROIs to match acuity image)
                MIJ.run("Multispectral Image Scaler No Scale Bar", ['scaling=' num2str(rescale)]);

                %send acuity image to ImageJ
                MIJ.createImage([mspecName 'acuity'],im_av,1);
                MIJ.run("32-bit")
                if strcmp(vs,'Cat')
                    MIJ.run("Set Label...", "label=lw")
                    MIJ.run("Next Slice [>]")
                    MIJ.run("Set Label...", "label=sw")
                    MIJ.run("Next Slice [>]")
                    MIJ.run("Set Label...", "label=lum") 
                else
                    MIJ.run("Set Label...", "label=lw")
                    MIJ.run("Next Slice [>]")
                    MIJ.run("Set Label...", "label=mw")
                    MIJ.run("Next Slice [>]")
                    MIJ.run("Set Label...", "label=sw")
                    MIJ.run("Next Slice [>]")
                    MIJ.run("Set Label...", "label=lum") 
                end

                %generate RNL ranked filter image
                MIJ.run('RNL Ranked Filter', ... 
                        [weberStr ' lum=' weberLum ' iterations=5 radius=' num2str(pxMRA) ' falloff=3']); 
                im_rnl = MIJ.getCurrentImage();
                saveImage(im2uint16(im_rnl), ... 
                          vs, .....
                          [baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_3-RNL'])

                %RNL clustering
                MIJ.run('RNL Clustering', ... 
                        ['colour_jnd_threshold=1.000 luminance_jnd_threshold=1.000 loops=20 radius=2 minimum=' num2str(pxMRA) ' compare=6 stop=1 record=20 image=' imName ' ' weberStr ' luminance_weber_fraction=' weberLum]);
                im_clust = MIJ.getImage([imName '_Clustered']);
                saveImage(im_clust, ... 
                          vs, ...
                          [baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_4-clustered'])

                im_clustID = MIJ.getImage([imName '_Cluster_IDs']);
                save([baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_clusterIDs.mat'], ... 
                     'im_clustID')
                imwrite(im_clustID./(max(im_clustID(:))), ... 
                          [baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_clusterIDs.tiff'])

                MIJ.selectWindow('Cluster Results');
                IJ.saveAs('results',[baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_clusterResults.csv']);

                %save log
                MIJ.selectWindow('Log');
                IJ.saveAs('text',[baseFolder foldDelim 'mica_results' foldDelim mspecName vs '_' illum '_' num2str(dist) '_log.txt']);

                %close windows (save first to prevent popup... SIGH)
                windows = MIJ.getListImages();
                for w = 1:length(windows)
                    MIJ.selectWindow(windows(w))
                    IJ.saveAs('tiff',[baseFolder foldDelim 'mica_results' foldDelim 'tmp.tiff']); 
                end
                MIJ.run('Clear Results')
                MIJ.closeAllWindows;

            end %visual systems loop
        end %illuminants loop
    end %images loop
end %distances loop

%% exit MIJ
MIJ.run('Clear Results')
MIJ.exit

%% compile rescale values
rescaleValues = CompileFiles([baseFolder foldDelim 'mica_results' foldDelim 'rescale_values'],'rescale');
rescaleValues = array2table(rescaleValues,'VariableNames',{'image','rescale'});
rescaleValues.image = strrep(rescaleValues.image,'_rescale','');
writetable(rescaleValues,[baseFolder foldDelim 'mica_results' foldDelim 'rescale_values.csv'])

%% calculate image sizes
ims = dir([baseFolder foldDelim 'mica_results' foldDelim '*2-AV_lum.tiff']);
imSize = [array2table({ims.name}','VariableNames',{'image'}), ... 
             array2table(zeros(length(ims),3),'VariableNames',{'nRow','nCol','nPx'})];
for i = 1:length(ims)
    im = imread(fullfile(ims(i).folder,ims(i).name));
    imSize.nRow(i) = size(im,1);
    imSize.nCol(i) = size(im,2);
    imSize.nPx(i) = numel(im);
    clear im
end
clear ims i

writetable(imSize,[baseFolder foldDelim 'mica_results' foldDelim 'image_size.csv'])
         
%basic image info
files = dir([micaDir foldDelim '*.mspec']);
files = {files.name};
files = strrep(files,'_.mspec','');
imInfo = table('Size',[length(files) 4], ...
             'VariableTypes',{'string', ...
                              'double', ...
                              'double', ...
                              'double'}, ...
             'VariableNames',{'image', ...
                              'scaleBar', ...
                              'imWidth', ...
                              'imHeight'});
for i = 1:length(files)
    imMatch = dir(fullfile(micaDir,'tiffs',[files{i} '*.tiff']));
    im = imread(fullfile(micaDir,'tiffs',imMatch.name));

    imInfo.image{i} = files{i};
    imInfo.imHeight(i) = size(im,1);
    imInfo.imWidth(i) = size(im,2);
    
    sbMatch = dir(fullfile(baseFolder,'ROIs',[files{i} '_Scale Bar*.roi']));
    imInfo.scaleBar(i) = str2double(strrep(strrep(sbMatch.name,[files{i} '_Scale Bar:'],''),':1.roi',''));

    clear imMatch im sbMatch
end

imInfo.nSquareMonkeys = (imInfo.imWidth./imInfo.scaleBar) .* (imInfo.imHeight./imInfo.scaleBar);

%% functions
function [] = saveImage(im,type,loc)
    %save original image as .mat file
    imwrite(im,[loc '.tiff'])
    
    %duplicate lw channel in dochromat images so chromatic image will be lw,lw,sw
    if strcmp(type,'Cat') 
        im = cat(3,im(:,:,1),im);
    end
    
    %save color & luminance images (w/ square root transformation)
    imwrite(real(sqrt(im2double(im(:,:,1:3)))),[loc '_col.tiff'])
    imwrite(real(sqrt(im2double(im(:,:,4)))),[loc '_lum.tiff'])
end

function [] = showImages(imDir,match,saveDir)
    ims = dir(fullfile(imDir,['*' match '*']));
    maxSz = 512;
    imMat = uint8(zeros(maxSz*maxSz*3,length(ims)));
    for i = 1:length(ims)
        %get image
        im = imread(fullfile(ims(i).folder,ims(i).name));
        
        %reduce resolution of large images
        if max(size(im))>maxSz 
            im = imresize(im,maxSz/max(size(im)));
        end
        
        %give greyscale images rgb values
        if size(im,3)==1
            im = repmat(im,[1 1 3]);
        end
        
        %add image to matrix, padding with 0s as necessary
        imAdd = uint8(zeros(maxSz,maxSz,3));
        offR = floor((maxSz-size(im,1))/2)+1;
        offC = floor((maxSz-size(im,2))/2)+1;
        imAdd(offR:offR+size(im,1)-1,offC:offC+size(im,2)-1,:) = im;
        imMat(:,i) = imAdd(:);
    end
    
    figure; montage(reshape(imMat,maxSz,maxSz,3,length(ims)))
    saveas(gcf,fullfile(saveDir,[match '.tiff']));
    close gcf
end

