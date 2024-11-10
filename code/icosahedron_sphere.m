function [vertices, faces] = icosahedron_sphere(refinement, plot_fig)
% Credit: ChatGPT (unknown)
    % Step 1: Create base icosahedron
    t = (1 + sqrt(5)) / 2;
    vertices = [
        -1,  t,  0;
         1,  t,  0;
        -1, -t,  0;
         1, -t,  0;
         0, -1,  t;
         0,  1,  t;
         0, -1, -t;
         0,  1, -t;
         t,  0, -1;
         t,  0,  1;
        -t,  0, -1;
        -t,  0,  1;
    ];
    % Normalize vertices to lie on the unit sphere
    vertices = vertices ./ vecnorm(vertices, 2, 2);
    
    % Define the 20 faces of the base icosahedron
    faces = [
        1, 12,  6;
        1,  6,  2;
        1,  2,  8;
        1,  8, 11;
        1, 11, 12;
        2,  6, 10;
        6, 12,  5;
        12, 11,  3;
        11,  8,  7;
        8,  2,  9;
        4, 10,  5;
        4,  5,  3;
        4,  3,  7;
        4,  7,  9;
        4,  9, 10;
        5, 10,  6;
        3,  5, 12;
        7,  3, 11;
        9,  7,  8;
        10,  9,  2;
    ];
    
    % Step 2: Refine the icosahedron
    for i = 1:refinement
        [vertices, faces] = refine_icosahedron(vertices, faces);
    end

    % Step 3: Plot the refined icosahedron sphere
    if plot_fig
        figure;
        trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'black');
        axis equal;
        title(['Icosahedron Sphere with Refinement Level ', num2str(refinement)]);
    end
end

function [vertices, faces] = refine_icosahedron(vertices, faces)
    % This function subdivides each triangle face into 4 smaller triangles
    edgeMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    newFaces = zeros(size(faces, 1) * 4, 3);
    newVertexIdx = size(vertices, 1) + 1;
    
    for i = 1:size(faces, 1)
        f = faces(i, :);
        
        % Create midpoints and add them if not present
        [a, edgeMap, vertices, newVertexIdx] = get_midpoint(f(1), f(2), edgeMap, vertices, newVertexIdx);
        [b, edgeMap, vertices, newVertexIdx] = get_midpoint(f(2), f(3), edgeMap, vertices, newVertexIdx);
        [c, edgeMap, vertices, newVertexIdx] = get_midpoint(f(3), f(1), edgeMap, vertices, newVertexIdx);
        
        % Create 4 new faces
        newFaces((i-1)*4 + 1, :) = [f(1), a, c];
        newFaces((i-1)*4 + 2, :) = [f(2), b, a];
        newFaces((i-1)*4 + 3, :) = [f(3), c, b];
        newFaces((i-1)*4 + 4, :) = [a, b, c];
    end
    
    faces = newFaces;
    vertices = vertices ./ vecnorm(vertices, 2, 2); % Normalize vertices to unit sphere
end

function [midpointIdx, edgeMap, vertices, newVertexIdx] = get_midpoint(v1, v2, edgeMap, vertices, newVertexIdx)
    % Generate a unique key for the edge
    edgeKey = sprintf('%d-%d', min(v1, v2), max(v1, v2));
    
    if isKey(edgeMap, edgeKey)
        % Edge midpoint already exists
        midpointIdx = edgeMap(edgeKey);
    else
        % Calculate new midpoint
        midpoint = (vertices(v1, :) + vertices(v2, :)) / 2;
        vertices(newVertexIdx, :) = midpoint;
        edgeMap(edgeKey) = newVertexIdx;
        midpointIdx = newVertexIdx;
        newVertexIdx = newVertexIdx + 1;
    end
end