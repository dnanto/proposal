#!/usr/bin/env awk -f

# 01 query acc.ver
# 02 subject acc.ver
# 03 % identity
# 04 alignment length
# 05 mismatches
# 06 gap opens
# 07 q. start
# 08 q. end
# 09 s. start
# 10 s. end
# 11 evalue
# 12 bit score
# 13 subject strand
# 14 query length

# plus: sstart - flank1,   send + flank2
# else:   send - flank2, sstart + flank1

/^[^#]/ { 
  f1 = $7 - 1; f2 = $14 - $8;
  if ($13 == "plus") { p1 = $9 - f1; p2 = $10 + f2; }
  else { p1 = $10 - f2; p2 = $9 + f1; }
  p1 = (p1 < 1) ? 1 : p1;
  printf "%s %d-%d %s\n", $2, p1, p2, $13;
}
