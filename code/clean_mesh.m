function [v,f] = clean_mesh(v,f,arg3)
%% remove unreferenced vertices
unref_v = setdiff(1:length(v),unique(reshape(f,1,3*length(f))));
if ~isempty(unref_v)
    for i = length(unref_v):-1:1
        f(f>=unref_v(i)) = f(f>=unref_v(i))-1;
    end
    v = v(setdiff(1:length(v),unref_v),1:3);
end

%% remove extra faces
[row,~] = find(f >length(v));
f(row,:) = [];

%% remove repeated faces
[~,unirow] = unique(sort(f,2), 'rows');
f = f(unirow,1:3);

%% remove zero area surfaces
e1 = v(f(:,3),:) - v(f(:,2),:);
e2 = v(f(:,1),:) - v(f(:,3),:);
a = cross(e1,e2);
area = ((a(:,1).^2 + a(:,2).^2 + a(:,3).^2).^(1/2))/2;
f = f(area>1e-15,1:3);

if nargin == 3
    %%
    for i = length(v):-1:1
        count = sum(sum((f == i)));
        if count == 1
            %% face with 1 connection 
            [row, ~] = find(f == i);
            f(row,:) = [];
            v(i,:) = [];
            f(f>=i) = f(f>=i)-1;
        end
    end
end