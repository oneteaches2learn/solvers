
function cMap = colormap_creator2()
% colormap downloaded from: https://www.kennethmoreland.com/color-advice/
% save it as a csv file (which is the default download) then just load it
% note: save the 'float' version of the csv file from the website

    input = readmatrix('fast1024');
    cMap = input(:,2:4);

    %cMap = flipud(grad_rgb);

end


