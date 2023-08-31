function new_data = upsample_data(sbci_map,data)

new_data = zeros(sbci_map.shape(4), 1);

for i = 1:length(data)
    new_data(sbci_map.map(1, sbci_map.map(2,:) == i)) = data(i);
end

end

