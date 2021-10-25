function yy=mysmooth(data,j)
if j==1 || j==length(data)
    yy=data(j);
elseif j==2 || j==length(data)-1
    yy=(data(j-1)+data(j)+data(j+1))/3;
else
    yy=(data(j-2)+data(j-1)+data(j)+data(j+1)+data(j+2))/5;
end
end