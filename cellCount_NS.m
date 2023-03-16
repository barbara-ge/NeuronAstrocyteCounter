%%cell count
clear all
close all

for i = 1 :2
    cd Example_images/

    filename1 = ['220311_NS_1_cs2_10X_' num2str(i) '_DAPI.tif'];
    filename = extractBefore(filename1,'DAPI');
    dapi = imread(filename1);

    [px, py] = size(dapi);
    % figure
    % subplot(1,2,1)
    % imshow(dapi)
    % subplot(1,2,2)
    % imhist(dapi)

    %%
    dapi = adapthisteq(dapi);

    % figure
    % imhist(dapi)

    %%
    dapiBorder = imclearborder(dapi);
    % figure
    % imshow(dapiBorder)

    %%
    dapiNoise = wiener2(dapiBorder, [5 5]);
    f1=figure(1);
    imshow(dapiNoise)
    saveas(f1,[filename, 'dapi_border.png'])

    %%
    dapiBW = im2bw(dapiNoise, graythresh(dapiNoise));
    f2=figure(2);
    imshow(dapiBW)
    saveas(f2,[filename, 'dapi_BW.png'])

    %%
    bw = imopen(dapiBW, strel('disk', 2));
    bw2 = bwareaopen(bw, 100);
    % figure
    % imshow(bw2)
    %%
    bw2_perim = bwperim(bw2);
    overlay1 = imoverlay(dapiBW, bw2_perim, [1 .3 .3]);

    % figure
    % imshow(overlay1)
    %%
    maxs = imextendedmax(dapiNoise, 5);
    maxs = imclose(maxs, strel('disk', 3));
    maxs = imfill(maxs, 'holes');
    maxs= bwareaopen(maxs,1);
    %maxs= bwareaopen(maxs,2);

    overlay2 = imoverlay(dapiBW, bw2_perim | maxs, [1 .3 .3]);

    % figure
    % imshow(overlay2)
    %%
    dapiJC = imcomplement(dapiNoise);
    dapi_mod = imimposemin(dapiJC, ~bw2 | maxs);
    f3=figure(3);
    imshow(dapi_mod)
    saveas(f3,[filename, 'dapi_maxima.png'])

    %%
    L = watershed(dapi_mod);
    labeledImage = label2rgb(L);

    f4=figure(4);
    imshow(labeledImage)
    saveas(f4, [filename,'dapi_labels.png'])

    %%
    [L, num] = bwlabel(L);

    mask = im2bw(L,1);
    overlay3 = imoverlay(dapiNoise, mask, [1 .3 .3]);

    f5=figure(5);
    imshow(overlay3)
    saveas(f5,[filename, 'dapi_detection_overlay.png'])

    %remove debris
    props = regionprops(L, 'PixelList');
    % L2 = L;
    % for cell = 2: num
    %     area = length(props(cell).PixelList);
    %     if area < 50
    %         for xx = 1: area
    %             c=props(cell).PixelList(xx,1);
    %             r=props(cell).PixelList(xx,2);
    %             L2(r,c)= 1;
    %         end
    %     end
    % end

    % [L2, num_clean] = bwlabel(L2);
    %
    % mask_debris = im2bw(L2,1);
    % overlay3 = imoverlay(dapiNoise, mask_debris, [1 .3 .3]);
    %%
    %load MAP2
    map2 = imread([filename, '647.tif']);
    [px, py] = size(map2);
    map2 = adapthisteq(map2);
    map2border = imclearborder(map2);
    %map2border = map2border.*2; %increase contrast
    % figure
    % imshow(map2border)
    map2Noise = wiener2(map2, [5 5]);
    %map2Noise = map2Noise.*2;
    map2BW = im2bw(map2Noise, graythresh(map2Noise));
    figure
    imshow(map2BW)

    black = zeros(size(map2), 'uint8');
    map2GR = cat(3, black, map2, black);
    figure
    imshow(map2GR)

    dapiAndMap2 = imfuse(map2, dapiNoise);
    % figure
    % imshow(dapiAndMap2)

    dapi_perim = imdilate(L,ones(3,3)) > imerode(L, ones(3,3));
    % figure
    % imshow(dapi_perim)

    map2GR = map2GR(:, :, 2);
    greenMF = medfilt2(map2GR, [3 3]);
    greenBW = im2bw(greenMF);
    % figure
    % imshow(greenMF)
    %
    % figure
    % imshow(greenBW)

    neuron = 0;
    astro = 0;
    for cell = 2: num
        one = [];
        [r, c] = find(L == cell);
        for row = 1: length(r)
            a = r(row);
            b = c(row);

            if map2BW(a,b)==1
                one(row,1) = 1;
            else
                one(row,1) = 0;
            end
        end

        if length(find(one == 1))>0.1*length(one)
            neuron = neuron+1;
            NeuronLoc(neuron) = cell;
        else
            astro = astro +1;
            AstroLoc(astro) = cell;
        end
    end

    maskAs = zeros(px,py);
    for ll = 1: length(AstroLoc)
        [r,c] = find(L == AstroLoc(ll));
        for rr = 1: length(r)
            a = r(rr);
            b = c(rr);
            maskAs(a,b) = 1;
        end
    end

    maskNeu = zeros(px,py);
    for ll = 1: length(NeuronLoc)
        [r,c] = find(L == NeuronLoc(ll));
        for rr = 1: length(r)
            a = r(rr);
            b = c(rr);
            maskNeu(a,b) = 1;
        end
    end

    LocNeu = imfuse(map2, maskNeu);
    f6=figure(6);
    imshow(LocNeu)
    saveas(f6,[filename, 'neurons_detection.png'])

    LocAs = imfuse(map2, maskAs);
    f7=figure(7);
    imshow(LocAs)
    saveas(f7,[filename, 'astro_detection.png'])

    LocTot = imfuse(LocNeu, maskAs);
    f8=figure(8);
    imshow(LocTot)
    saveas(f8,[filename, 'tot_detection.png'])

    %save('counts_NS_B1__.mat', 'neuron', 'NeuronLoc', 'astro', 'AstroLoc', 'maskNeu', 'maskAs', 'map2GR', 'dapiBW', 'dapiNoise')
    %%
    %GFAP
    gfap = imread([filename, '555.tif']);
    [px, py] = size(gfap);


    for i = 1: px
        for j = 1: py

            if gfap(i,j) < 30
                gfap(i,j) = 0;
            end
        end
    end
    gfap = adapthisteq(gfap);
    gfapborder = imclearborder(gfap);
    gfapNoise = wiener2(gfap, [2 2]);
    gfapBW = im2bw(gfapNoise, graythresh(gfapNoise));
    % figure
    % imshow(gfapBW)

    black = zeros(size(gfap), 'uint8');
    gfapRD = cat(3, gfap, black, black);
    figure
    imshow(gfapRD)

    dapiAndgfap = imfuse(gfapRD, dapiNoise);
    % figure
    % imshow(dapiAndgfap)

    astro2 = 0;
    for cell2 = 2: astro
        cell = AstroLoc(1,cell2);
        one_1 = [];
        [r, c] = find(L == cell);
        for row = 1: length(r)
            a = r(row);
            b = c(row);

            if gfapBW(a,b)==1
                one_1(row,1) = 1;
            else
                one_1(row,1) = 0;
            end
        end

        if length(find(one_1 == 1))>0.71*length(one_1)
            astro2 = astro2+1;
            AstroLoc2(astro2) = cell;
        end
    end

    maskAs2 = zeros(px,py);
    for ll = 1: length(AstroLoc2)
        [r,c] = find(L == AstroLoc2(ll));
        for rr = 1: length(r)
            a = r(rr);
            b = c(rr);
            maskAs2(a,b) = 1;
        end
    end

    LocAs2 = imfuse(gfapRD, maskAs2);
    f9=figure(9);
    imshow(LocAs2)
    saveas(f9,[filename, 'astro_detection_clean.png'])

    LocAs3 = imfuse(gfapRD, maskAs);
    f10=figure(10);
    imshow(LocAs3)
    saveas(f10,[filename, 'astro_detection_not_clean.png'])
    % dapi_perim = imdilate(L,ones(3,3)) > imerode(L, ones(3,3));
    % figure
    % imshow(dapi_perim)

    ratioA = astro2/cell;
    ratioN = neuron/cell;
    ratioA2 = astro2/(astro2+neuron)
    ratioN2 = neuron/(astro2+neuron)

    cd ..
    save([filename, 'counts.mat'],'L', 'neuron', 'NeuronLoc', 'astro', 'astro2', 'AstroLoc2', 'AstroLoc', 'maskNeu', 'maskAs', 'map2GR', 'dapiBW', 'dapiNoise', 'gfapRD', 'gfapBW', 'maskAs2', 'LocAs2')
end


