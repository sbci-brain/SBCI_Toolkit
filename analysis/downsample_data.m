function new_data = downsample_data(sbci_map,data)

new_data = zeros(sbci_map.shape(1), 1);

for i = 1:length(new_data)
    new_data(i) = mean(data(sbci_map.map(1, sbci_map.map(2,:) == i)));
end

end

