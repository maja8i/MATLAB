% De-multiplex bits of Nx1 uint64 mortonCode into x, y, z
% (each having depth<=16 bits).
function xyz = mortonToXyz(mortonCode,depth)
    xyz = zeros(size(mortonCode,1),3,'uint16');
    for b = 1:depth
        xyz(:,1) = bitset(xyz(:,1),b,bitget(mortonCode,3*(b-1)+3));
        xyz(:,2) = bitset(xyz(:,2),b,bitget(mortonCode,3*(b-1)+2));
        xyz(:,3) = bitset(xyz(:,3),b,bitget(mortonCode,3*(b-1)+1));
    end
end