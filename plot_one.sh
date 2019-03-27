#!/bin/bash 

plot_one() {
gnuplot<<EOF
set terminal png
set output '$1.png'
plot '$1' matrix with image 
EOF
}

for i in `ls *.dat` ; do
  if [[ "$1" != "noskip" ]]; then
    if [[ -e $i.png ]]; then
      echo "skipping $i"
      continue;
    fi
  fi
  echo "plotting $i"
  plot_one $i  
done
