#
foreach i (0 1)
foreach j (0 1 2 3 4 5 6 7 8 9)
insBCDll ../../BCD/run/data/watBrO/watBrO{$i}{$j}44.rsp ../../BCD/run/inpfiles/BCDb/BCDb{$i}{$j}4.equ 3944 -2.0 Bll/bigBll{$i}{$j}0.as
writebin Bll/bigBll{$i}{$j}0.as Bll/bigBll{$i}{$j}0.min
psxyz Bll/bigBll{$i}{$j}0.min Bll/bigBll{$i}{$j}0.xyz
cutLL Bll/bigBll{$i}{$j}0.min Bll/equBll{$i}{$j}0.as 50.0 980 1
writebin Bll/equBll{$i}{$j}0.as Bll/equBll{$i}{$j}0.min
psxyz Bll/equBll{$i}{$j}0.min Bll/equBll{$i}{$j}0.xyz
end
end
