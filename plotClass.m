function plotClass(X,species,class,alpha,b,K)
    y = getY(species,class);
    x1 = linspace(0,8);
    x2 = linspace(1,5);
    [X1,X2] = meshgrid(x1,x2);
    Z = zeros(length(x1),length(x2));
    for i = 1:length(x2)
        for j = 1:length(x1)
            Z(i,j) = (value(X,y,alpha,K,b,[x1(j),x2(i)]));
        end
    end
    range = [-inf 0 inf];
    contourf(X1,X2,Z,range); hold on;
    xc1 = X(find(y+1),:);
    xc2 = X(find(y-1),:);
    p1 = plot(xc1(:,1),xc1(:,2),"go",'DisplayName','cos(3x)'); hold on;
    p2 = plot(xc2(:,1),xc2(:,2),"ro",'DisplayName','cos(2x)');
    xlabel('Pedal Length');
    ylabel('Sepal Width');
    legend([p1 p2],{string(class), 'other'});
end