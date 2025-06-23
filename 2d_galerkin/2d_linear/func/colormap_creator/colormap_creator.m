
function cMap = colormap_creator()
    % values created by https://davidjohnstone.net/lch-lab-colour-gradient-picker
    % note: manually edited to change '[' and ']' to '{' and '}'
    grad_hex = {'#fcbaba', '#fab9ba', '#f7b9ba', '#f5b8ba', '#f2b7bb', '#f0b6bb', '#edb6bb', '#ebb5bb', '#e8b4bb', '#e6b4bb', '#e3b3bb', '#e0b2bc', '#deb2bc', '#dbb1bc', '#d9b0bc', '#d6afbc', '#d4afbc', '#d1aebc', '#ceadbd', '#ccadbd', '#c9acbd', '#c7abbd', '#c4aabd', '#c1aabd', '#bfa9bd', '#bca8bd', '#b9a8be', '#b7a7be', '#b4a6be', '#b1a5be', '#aea5be', '#aca4be', '#a9a3be', '#a6a2be', '#a3a2be', '#a0a1bf', '#9da0bf', '#9aa0bf', '#989fbf', '#959ebf', '#929dbf', '#8e9dbf', '#8b9cbf', '#889bbf', '#859ac0', '#829ac0', '#7f99c0', '#7b98c0', '#7897c0', '#7497c0', '#7196c0', '#6d95c0', '#6a95c0', '#6694c0', '#6293c1', '#5e92c1', '#5992c1', '#5591c1', '#5090c1', '#4c8fc1', '#498ec0', '#488cbe', '#478bbb', '#4689b9', '#4587b7', '#4486b5', '#4384b2', '#4282b0', '#4181ae', '#417fac', '#407da9', '#3f7ca7', '#3e7aa5', '#3d79a3', '#3c77a1', '#3b759e', '#3a749c', '#39729a', '#397098', '#386f96', '#376d93', '#366c91', '#356a8f', '#34688d', '#33678b', '#336589', '#326487', '#316284', '#306082', '#2f5f80', '#2e5d7e', '#2e5c7c', '#2d5a7a', '#2c5978', '#2b5776', '#2a5674', '#295472', '#295270', '#28516d', '#274f6b', '#264e69', '#254c67', '#244b65', '#244963', '#234861', '#22465f', '#21455d', '#20435b', '#204259', '#1f4057', '#1e3f55', '#1d3d53', '#1c3c51', '#1c3b4f', '#1b394e', '#1a384c', '#19364a', '#193548', '#183346', '#173244'};

    % Convert hex to RGB
    for i = 1:length(grad_hex)
        grad_rgb(i,:) = hex2rgb(grad_hex{i});
    end

    cMap = flipud(grad_rgb);

end


