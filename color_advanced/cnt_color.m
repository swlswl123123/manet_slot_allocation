function [color_set, color_cnt] = cnt_color(src, dst, link)
%cnt_color - Description
%
% Syntax: [color_set, color_cnt] = cnt_color(src, dst, link)
%
% Long description
    color_set = [];
    for i = 1:size(link, 1)
        if link(i,3) ~= 0 && (link(i,1) == src || link(i,2) == src || link(i,1) == dst || link(i,2) == dst)
            color_set = [color_set, link(i,3)];
        end
    end
    color_cnt = length(color_set);
end