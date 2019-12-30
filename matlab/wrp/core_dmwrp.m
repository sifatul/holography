project_name = 'dmwrp_center';

segmentation_number = 0;

while segmentation_number < total_wrp %first step of segmentation
    if segmentation_number == 0
        [Cut,~,idx] = unique(obj_z);
        n = accumarray(idx(:),1);
        X = [0,~diff(Cut',3),0];
        [B,E] = regexp(char('0'+X),'0+1*');
        C = arrayfun(@(b,e){Cut(b:e)},B,E);
    else
        
        %      B_temp = B;
        %      E_temp = E;
        [s,ddif] = cellfun(@size,C);
        [~,index_max] = max([s,ddif]);
        
        second_order_obj_z = Cut(B(index_max):E(index_max));
        
        A = obj_z( ismember( obj_z, second_order_obj_z ) );
        
        [~,~,idx2] = unique(A);
        n = accumarray(idx2(:),1);
        X = [0,1,~diff(n',2)];
        [B_new,E_new] = regexp(char('0'+X),'0+1*');
        C_new = arrayfun(@(b,e){n(b:e)},B_new,E_new);
        
        B_new = B_new + E(index_max-1);
        E_new = E_new +  E(index_max-1);
        
        %
        % n = [ 3 3 1 1 1 ];
        % X = [0,1,~diff(n,2)];
        
        
        
        B_temp = [B(1:index_max-1), B_new(1:end), B(index_max+1:end)];
        E_temp = [E(1:index_max-1), E_new(1:end), E(index_max+1:end)];
        C_temp = {C{1:index_max-1},C_new{1:end}, C{index_max+length(C_new)-1:end}};
        C(1:end) = [];     B(1:end) = [];     E(1:end) = [];
        C = C_temp;  B = B_temp;  E = E_temp;
        
        if length(C) ~= length(B)
            disp("######");
            return;
        end
        
        
    end
    %%
    
    segmentation_number = segmentation_number+1;
    
end


% for i = 1:length(C)       %iterate each depth range
%     range_range_all_index = B(i):E(i) ;
%     range = Cut(range_range_all_index);
%     ranges{i, 1} =range;
% end
for i = 1:length(C)       %iterate each depth range
    range_range_all_index = B(i):E(i) ;
    range = Cut(range_range_all_index);
    if length(range)<=3
        Lia = ismember(obj_z,range);
        idx = find(Lia);
        obj_z(idx) = median(range);
%         ranges{i, 1} = median(range);
    else
%         ranges{i, 1} = range;
    end
    
end
% obj_z = cell2mat(ranges);
unique(obj_z)

