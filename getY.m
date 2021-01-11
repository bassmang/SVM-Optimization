function y = getY(species,class)
   y = zeros(length(species),1);
    for i = 1:length(y)
        if strcmp(species(i),class)
            y(i) = 1;
        else
            y(i) = -1;
        end
    end
    y = y';
end