function chanlocsOut = bvef_to_chanlocs(bvefFile)

    doc = xmlread(bvefFile);
    nodes = doc.getElementsByTagName('Electrode');
    n = nodes.getLength();

    labels = cell(n,1);
    Theta = nan(n,1); Phi = nan(n,1); Radius = nan(n,1);
    valid = false(n,1);

    for i = 0:n-1
        node = nodes.item(i);
        try
            labels{i+1} = char(node.getElementsByTagName('Name').item(0).getFirstChild.getData);
            Theta(i+1)  = str2double(node.getElementsByTagName('Theta').item(0).getFirstChild.getData);
            Phi(i+1)    = str2double(node.getElementsByTagName('Phi').item(0).getFirstChild.getData);
            Radius(i+1) = str2double(node.getElementsByTagName('Radius').item(0).getFirstChild.getData);
            valid(i+1)  = true;
        catch
            valid(i+1) = false;  
        end
    end

    labels = labels(valid); Theta = Theta(valid); Phi = Phi(valid); Radius = Radius(valid);

    % --- Brainstorm's verified spherical -> Cartesian conversion ---
    phi_rad   = deg2rad(180 + Phi);
    theta_rad = deg2rad(90 - Theta);
    sx = Radius .* cos(theta_rad) .* cos(phi_rad);
    sy = Radius .* cos(theta_rad) .* sin(phi_rad);
    sz = Radius .* sin(theta_rad);

    X = -sy;        % toward nose
    Y =  sx;        % toward left ear
    Z =  sz + 0.05; % toward vertex

    chanlocsOut = struct('labels', labels(:), ...
        'X', num2cell(X(:)), ...
        'Y', num2cell(Y(:)), ...
        'Z', num2cell(Z(:)));


    if exist('convertlocs', 'file')
    chanlocsOut = convertlocs(chanlocsOut, 'cart2all');
    else
        warning(['bvef_to_chanlocs: convertlocs (EEGLAB) not found on path -- ' ...
            'returning X/Y/Z only, without 2-D theta/radius/sph_* fields. ' ...
            'Run eeglab (or eeglab nogui) first if you need those.']);
    end

    fprintf('bvef_to_chanlocs: parsed %d electrode positions from %s\n', ...
        numel(labels), bvefFile);

    fprintf('bvef_to_chanlocs: parsed %d electrode positions from %s\n', ...
        numel(labels), bvefFile);
end