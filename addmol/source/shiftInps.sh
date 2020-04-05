foreach i (0 1)
foreach j (0 1 2 3 4 5 6 7 8 9)
zshift ../../BCD/run/inpfiles/watBrO/smWatBrO{$i}{$j}3.equ ../../BCD/run/inpfiles/watBrO/smWatBrO{$i}{$j}4.as 0 2.5 1
writebin ../../BCD/run/inpfiles/watBrO/smWatBrO{$i}{$j}4.as ../../BCD/run/inpfiles/watBrO/smWatBrO{$i}{$j}4.equ 
psxyz ../../BCD/run/inpfiles/watBrO/smWatBrO{$i}{$j}4.equ ../../BCD/run/inpfiles/watBrO/smWatBrO{$i}{$j}4.xyz
end
end
