% Multiplex bits of Nx3 uint16 xyz (each having depth<=16 bits).
function mortonCode = xyzToMorton(xyz,depth)
    mortonCode = zeros(size(xyz,1),1,'uint64');
    for b = 1:depth
        mortonCode = bitset(mortonCode,3*(b-1)+3,bitget(xyz(:,1),b));
        mortonCode = bitset(mortonCode,3*(b-1)+2,bitget(xyz(:,2),b));
        mortonCode = bitset(mortonCode,3*(b-1)+1,bitget(xyz(:,3),b));
    end
end