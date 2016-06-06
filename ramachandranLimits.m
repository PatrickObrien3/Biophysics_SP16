function ramachandranLimits()
% RAMACHANDRANLIMITS Plot the hard-sphere limits of the Ramachandran map
% to the current figure.
%
% The current figure is expected to be in radians.
    file_name='tau115data.csv';
    tempMat=csvread(file_name);
    ind=find(tempMat(:,2)>180);
    tempMat(ind,2)=tempMat(ind,2)-360;
    tempMat(3,2)=180;

    region1=tempMat(3:26,:);
    region2=tempMat([27:28,1:2],:);
    region3=tempMat([29:34,29],:);
    region4=tempMat([35:41,35],:);
    region5=tempMat([42:45,42],:);
    region6=tempMat(35:38,:);

    hold all
    plot(region1(:,1), region1(:,2), '-m', 'LineWidth',2)
    plot(region2(:,1), region2(:,2), '-m', 'LineWidth',2)
    plot(region3(:,1), region3(:,2), '-b', 'LineWidth',2)
    plot(region5(:,1), region5(:,2), '-m', 'LineWidth',2)
    plot(region4(:,1), region4(:,2), '-b', 'LineWidth',2)
    plot(region6(:,1), region6(:,2), '-b', 'LineWidth',2)
end