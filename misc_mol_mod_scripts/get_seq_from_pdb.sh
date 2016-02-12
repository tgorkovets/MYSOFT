#!/bin/bash
#awk 'BEGIN {pr = -1; pch = "qw"; } { cr = $6; cch = $5; } pch != cch; { printf( "chain %s: \n", cch); pch = cch; } pr != cr; { printf("%s",$4); } {pr = cr;} END {print "";} ' $1;
awk 'BEGIN {pr = -1; pch = 4; } { cr = $6; cch = $5;  if(pch != cch) { printf( "\n chain %s: \n", cch); pch = cch;} if(pr != cr){printf("%s",$4); pr = cr;}  } END {print " ";} ' $1;
