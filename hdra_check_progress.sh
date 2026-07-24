#!/bin/bash
cd /home/yehia/hdra_rdp1
done=0
todo=0
for s in hdra_S_none hdra_S_010 hdra_S_005 hdra_M_none hdra_M_001 hdra_M_003 hdra_M_005 hdra_M_010 hdra_G_none hdra_G_020 hdra_G_010 hdra_G_005 hdra_G_002 hdra_L_none hdra_L_08 hdra_L_05 hdra_L_02 hdra_L_100; do
  if [ -f "${s}_K10.Q" ]; then
    echo "OK  $s"
    done=$((done+1))
  else
    echo "     $s"
    todo=$((todo+1))
  fi
done
echo "Done: $done/18  Remaining: $todo"