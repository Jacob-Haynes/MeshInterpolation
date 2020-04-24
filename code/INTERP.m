%% This script runs mesh interperlation

%% Load start and target mesh
obj = readmesh('man.obj');
% number of verticies and faces
Vnum = size(obj.v,1);
Fnum = size(obj.f,1);
% NOTE: WOULD ALSO DO THIS WITH THE TARGET MESH
%   HOWEVER HERE A STRUCT FROM PREVIOUS CW IS USED AND IS ALREADY EXPORTED
%   FROM .OBJ TO A MATLAB DATA STRUCT
target = importdata('end_mesh2.mat'); %this file is from last cw

%% linear interpolation of vertices
%paramaters
frames = 10;
%position(t) = (1-t)position(end) + (t)position(start)
pos = cell(frames,1);
pos{1} = obj.v;
for t = 1:frames
    tf = 1-(1/frames)*t;
    pos{t+1} = (1-tf).*target.v + (tf).*obj.v;
end

animate(pos, obj);


%% transformation-based interpolation
% compute A for each triangle
% Ap = q
A = cell(length(obj.f),1);
B = cell(length(obj.f),1);
Ra = cell(length(obj.f),1);
Rb = cell(length(obj.f),1);
Rg = cell(length(obj.f),1);
D = cell(length(obj.f),1);
S = cell(length(obj.f),1);
a = zeros(2,6);
b = zeros(6,1);
for i = 1:length(obj.f)
    %Vals
    At = cell(frames,1);
    V_idx = obj.f(i,:);
    q = target.v(V_idx,:);
    p = obj.v(V_idx,:);
    p(:,3)=1;
    for j = 1:3
        a((2*j-1),1:3)= p(j,:);
        a((2*j),4:6) = p(j,:);
        b(2*j-1) = q(j,1);
        b(2*j) = q(j,2);
    end
    
    %ax=q' : ax-q' minimise
    x = a\b;
    A{i} = [x(1:2)';x(4:5)']; %only care about rotation
    
    %SVD of A to get Rgamma
    [Ra{i},D{i},Rb{i}] = svd(A{i});
    Rg{i} = Ra{i}*Rb{i};
    S{i} = Rb{i}'*D{i}*Rb{i};
    
    %linerarly interpolate to get A_t for each triangle
    for f = 1:frames
        tf = 1-(1/frames)*f;
        At{f+1} = Rg{i}*((1-tf).*eye(2) + (tf).*S{i});
    end
    A{i} = At; %A{T}{f} = A for each triangle T and each frame f
end

% compute V for each triangle
V = cell(frames,1);
V{1} = obj.v;
V{frames} = target.v;
for f = 2:(frames-1)
    uc = [];
    bc = [];
    for i = 1:length(obj.f)
        V_idx = obj.f(i,:);
        c1 = V{f-1}(V_idx(1),:);
        c2 = V{f-1}(V_idx(2),:);
        c3 = V{f-1}(V_idx(3),:);
        c1(3) = 1;
        c2(3) = 1;
        c3(3) = 1;
        M = inv([c1,0,0,0;0,0,0,c1;c2,0,0,0;0,0,0,c2;c3,0,0,0;0,0,0,c3]);
        M = [M(1,:);M(2,:);M(4,:);M(5,:)];
        u = zeros(length(obj.v)*2,4);
        for j = 1:4
            u(V_idx(1)*2-1:V_idx(1)*2,j) = M(j,1:2);
            u(V_idx(2)*2-1:V_idx(2)*2,j) = M(j,3:4);
            u(V_idx(3)*2-1:V_idx(3)*2,j) = M(j,5:6);
        end
        % fix first vertex
        u = u(3:end,:)';
        uc = [uc;u];
        %form bd
        Ad = A{i}{f};
        bd = [Ad(1,1),Ad(1,2),Ad(2,1),Ad(2,2)]';
        bc = [bc;bd];

    end
    V_target = uc\bc;
    V{f}= pos{f}(1,:); %linearly interpolate first vertex from part 1
    for r = 1:2:(length(V_target)-1)
        rw = [V_target(r),V_target(r+1),0];
        V{f} = [V{f};rw];
    end
end

