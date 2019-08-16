% function depth_ranges = depth_segmentation_conventional(Cut,total_wrp)
% 
% len = numel(Cut);
% 
% % 
% size_of_one_depth_range = round(length(Cut) / total_wrp);
% n=numel(Cut);
% m=fix(n/size_of_one_depth_range);
% p=mod(n,size_of_one_depth_range); 
% depth_ranges =[mat2cell(Cut(1:m*size_of_one_depth_range),size_of_one_depth_range*ones(m,1),1);{Cut(size_of_one_depth_range*m+1:size_of_one_depth_range*m+p)}];
% 
% if(cellfun(@isempty, depth_ranges(end)))
%     depth_ranges(end)=[];
% end
%  
% end
function depth_ranges = depth_segmentation_conventional(Cut,total_wrp)

len = numel(Cut);
% if mod(len,2)==1
%     Cut(end+1)=100;
% end
% 
size_of_one_depth_range = round(length(Cut) / total_wrp);
n=numel(Cut);
m=fix(n/size_of_one_depth_range);
p=mod(n,size_of_one_depth_range); 

if mod(len,2)==1
       last = size_of_one_depth_range*m+p;
else
    last = size_of_one_depth_range*m+p;
end

depth_ranges =[mat2cell(Cut(1:m*size_of_one_depth_range),size_of_one_depth_range*ones(m,1),1);{Cut(size_of_one_depth_range*m+1:last)}];

if(cellfun(@isempty, depth_ranges(end)))
    depth_ranges(end)=[];
end
 
end
