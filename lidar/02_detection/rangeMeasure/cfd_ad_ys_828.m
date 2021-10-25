function [x]=cfd_ad_ys_828(data,max_index,s,y)
%% 到处都是阈值
if data(max_index)<=5
    x=NaN;
    return
end
if data(max_index)<=40
    min_d=8;
else
    min_d=20;
end
data1=zeros(length(data)+y,1);
data2=zeros(length(data)+y,1);
data1(1:length(data))=floor(s*data);
data1(length(data)+1:length(data)+y)=data1(length(data));
data2(1:y,1)=data(1);
data2(y+1:length(data)+y,1)=data(1:length(data));
% data2(1:1024+y,1)=data(1-y:1024);
% data2(1024+y+1:1024,1)=data(1024);
diff_data=data2-data1;
for i=max_index+y:-1:3
    try
        if data2(i)<=min_d
            try 
                if data2(i-1)-data2(i)<=0 && data2(i-2)-data2(i)<=0
                    break
                else
                    continue
                end
            catch
            break
            end
        end
    catch
        x=NaN;
        return
    end
end
for j=i:length(data1)-1
    if diff_data(j)<=0 && diff_data(j+1)>=0
        break
    end
end
%     x_d=1:length(data);
%         figure(100)
%         plot(x_d,data1)
%         hold on
%         plot(x_d,data2)
%         hold on
if ~isempty(j)
    if j>=length(data) || data(j)<5
        x=NaN;
    else
        s_data_1=mysmooth(data,j+1);
        s_data_2=mysmooth(data,j);
        s_data2_1=mysmooth(data2,j+1);
        s_data2_2=mysmooth(data2,j);
        k1=(s_data_1-s_data_2)*s;
        b1=s_data_2*s-k1*j;
        k2=s_data2_1-s_data2_2;
        b2=s_data2_2-k2*j;
        x=-(b1-b2)/(k1-k2);
%         s_data1=smooth(data);
%         s_data2=smooth(data2);
%         k1=(s_data1(j+1)-s_data1(j))*s;
%         b1=s_data1(j)*s-k1*j;
%         k2=s_data2(j+1)-s_data2(j);
%         b2=s_data2(j)-k2*j;
%         x=-(b1-b2)/(k1-k2)-63.4;
    end
else
    x=NaN;
end
end
