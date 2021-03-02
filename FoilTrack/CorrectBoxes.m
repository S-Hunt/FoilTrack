files = ls('*boxes.mat')


for x = 1:length(files)
    
    load(files(x,:))
    boxX = round(boxX);
    if boxX(1,2) == boxX(2,1)
        boxX(:,2) = boxX(:,2) - 1;

    end
    boxX
    save(files(x,:), 'boxX', '-append')
end