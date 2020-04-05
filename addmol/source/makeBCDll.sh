foreach j (0 1)
foreach i (0 1 2 3 4 5 6 7 8 9)
foreach z (000 002 004 006 008 010 012 014 -02 -04 -06 -08 -10 -12 -14)
insBCDll ../../BCD/run/inpfiles/watBrO/watBrO{$j}{$i}2.equ BCDfromPDB.min 3944 {$z} ../../BCD/run/inpfiles/BCDll/BCDll{$z}_{$j}{$i}0.as
writebin ../../BCD/run/inpfiles/BCDll/BCDll{$z}_{$j}{$i}0.as ../../BCD/run/inpfiles/BCDll/BCDll{$z}_{$j}{$i}0.min
psxyz ../../BCD/run/inpfiles/BCDll/BCDll{$z}_{$j}{$i}0.min ../../BCD/run/inpfiles/BCDll/BCDll{$z}_{$j}{$i}0.xyz
end
end
end
