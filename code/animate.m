function animate(film, obj)
% play the animation at a designated frame rate
% film is a cell array of the vertex coords for each frame
% obj is the mesh - only for the struct and face component

for t = 1:length(film)
    obj.v = film{t};
    cla
    %xlim([-1, 1])
    %ylim([-1, 1.5])
    dispmodel(obj);
    pause(0.5);
end
end
