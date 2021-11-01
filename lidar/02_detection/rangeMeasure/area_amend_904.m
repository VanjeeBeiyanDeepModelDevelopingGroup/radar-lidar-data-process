function [area_data,wide_data,max1,x_dif,x_dif2,x_dif3,x_max,wide_data2,area_data2,x_dif4]=area_amend_904(data,max_index)
if data(max_index)<20 || length(data)-max_index<3 || max_index<4
    area_data=NaN;
    wide_data=NaN;
    max1=NaN;
    x_dif=NaN;
    x_dif2=NaN;
    x_dif3=NaN;
    x_dif4=NaN;
    x_max=NaN;
    wide_data2=NaN;
    area_data2=NaN;
    return
end
wide1=20;
for i=max_index:-1:1
    if data(i)<=wide1
        try 
            if max(data(i-1),data(i-2))-data(i)<10 && max(data(i-1),data(i-2))-min(data(i-4),data(i-3))<10  %Ê¶±ðÕñµ´²¨ÐÎ
                break
            else
                continue
            end
        catch
            break
        end
    end
end
if data(i)>wide1 || length(data)-i<=1
    area_data=NaN;
    wide_data=NaN;
    max1=NaN;
    x_dif=NaN;
    x_dif2=NaN;
    x_dif3=NaN;
    x_dif4=NaN;
    x_max=NaN;
    wide_data2=NaN;
    area_data2=NaN;
    return
    
else
    data_i=mysmooth(data,i);
    data_i1=mysmooth(data,i+1);
end

wide_t1=i+(wide1-data_i)/(data_i1-data_i);
for j=max_index:length(data)
    if data(j)<=wide1
        try 
            if max(data(j+1),data(j+2))-data(j)<10 && max(data(j+1),data(j+2))-min(data(j+4),data(j+3))<10
                break
            else
                continue
            end
        catch
            break
        end
    end
end
if data(j)>wide1 || j-1<=1
    area_data=NaN;
    wide_data=NaN;
    max1=NaN;
    x_dif=NaN;
    x_dif2=NaN;
    x_dif3=NaN;
    x_dif4=NaN;
    x_max=NaN;
    wide_data2=NaN;
    area_data2=NaN;
    return
else
data_j=mysmooth(data,j);
data_j1=mysmooth(data,j-1);
end
wide_t2=j-(wide1-data_j)/(data_j1-data_j);
wide_data=wide_t2-wide_t1;
area1=20;
for ii=max_index:-1:1
    if data(ii)<=area1
        try 
            if max(data(ii-1),data(ii-2))-data(ii)<10 && max(data(ii-1),data(ii-2))-min(data(ii-4),data(ii-3))<10
                break
            else
                continue
            end
        catch
            break
        end
    end
end
for jj=max_index:length(data)
    if data(jj)<=area1
        try 
            if max(data(jj+1),data(jj+2))-data(jj)<10 && max(data(jj+1),data(jj+2))-min(data(jj+4),data(jj+3))<10
                break
            else
                continue
            end
        catch
            break
        end
    end
end
if data(ii)>area1 || data(jj)>area1
    area_data=NaN;
    wide_data=NaN;
    max1=NaN;
    x_dif=NaN;
    x_dif2=NaN;
    x_dif3=NaN;
    x_dif4=NaN;
    x_max=NaN;
    wide_data2=NaN;
    area_data2=NaN;
    return
end
area_data=0;
for k=ii+1:jj-1
    if k==ii+1||k==jj-1
        area_data=area_data+0.5*(data(k)-area1);
    else
        area_data=area_data+(data(k)-area1);
    end
end

for mm=max_index:-1:1
    if data(mm)<=10
        try 
            if max(data(mm-1),data(mm-2))-data(mm)<10 && max(data(mm-1),data(mm-2))-min(data(mm-4),data(mm-3))<10
                break
            else
                continue
            end
        catch
            break
        end
    end
end
wide_10=mm+(10-data(mm))/(data(mm+1)-data(mm));
for kk=mm:max_index
    if data(kk+1)-data(kk)<=0
        break
    end
end
max1=data(kk);
x_max=kk;
if data(max_index)>40
    for kkk=mm:max_index
        if data(kkk)>=40
            break
        end
    end
    try
        data_kkk=data(kkk);
        data_kkk1=data(kkk-1);
        wide_40=kkk-(data_kkk-40)/(data_kkk-data_kkk1);
    catch
        area_data=NaN;
        wide_data=NaN;
        max1=NaN;
        x_dif=NaN;
        x_dif2=NaN;
        x_dif3=NaN;
        x_max=NaN;
        wide_data2=NaN;
        area_data2=NaN;
        x_dif4=NaN;
        return
    end
    x_dif=wide_40-wide_10;
    if max1>60
        for m=kkk:max_index
            if data(m)>=60
                break
            end
        end
        data_m=data(m);
        data_m1=data(m-1);
        wide_60=m-(data_m-60)/(data_m-data_m1);
        x_dif3=wide_60-wide_40;
    else 
        x_dif3=NaN;
    end
    if max1>70
        for m=mm:max_index
            if data(m)>=70
                break
            end
        end
        data_m=data(m);
        data_m1=data(m-1);
        wide_70=m-(data_m-70)/(data_m-data_m1);
        x_dif2=wide_70-wide_10;
    else 
        x_dif2=NaN;
    end
else
    max1=NaN;
    x_dif=NaN;
    x_dif2=NaN;
    x_dif3=NaN;
end
if data(max_index)>70
    for i1=max_index:-1:1
        if data(i1)<=70
            break
        end
    end
    wide_t3=i1+(70-data(i1))/(data(i1+1)-data(i1));
    for j1=max_index:length(data)
        if data(j1)<=10
            break
        end
    end
    data_j=mysmooth(data,j1);
    data_j1=mysmooth(data,j1-1);
    wide_t4=j1-(10-data_j)/(data_j1-data_j);
    wide_data2=wide_t4-wide_t3;
    
    area_data2=0.5*(i1+1-wide_t3)*(70+data(i1+1)-20);
    for k=i1+1:j1-1
        if k==i1+1||k==j1-1
            area_data2=area_data2+0.5*(data(k)-10);
        else
            area_data2=area_data2+(data(k)-10);
        end
    end
else
    wide_data2=NaN;
    area_data2=NaN;
end
if data(max_index)>80
    for i1=max_index:-1:1
        if data(i1)<=80
            break
        end
    end
    wide_80=i1+(80-data(i1))/(data(i1+1)-data(i1));
    for j1=max_index:-1:1
        if data(j1)<=60
            break
        end
    end
    wide_60_2=j1+(60-data(j1))/(data(j1+1)-data(j1));
    x_dif4=wide_80-wide_60_2;
else
    x_dif4=NaN;
end
end
