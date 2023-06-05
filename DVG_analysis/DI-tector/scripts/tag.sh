#!/bin/bash

for i in [Bl]*[1-3] ; do 
    echo $i
    # link to new name (Logical -dereference symlinks-, relative, symbolic) 
    ln -Lrs $i `basename $i \
	      | sed \
              -e 's/_.*R/_R/g' \
	      -e 's/lnc0/000.0pfu/g' \
	      -e 's/lncM/mock/g' \
	      -e 's/lnc1/100.0pfu/g' \
	      -e 's/lnc2/010.0pfu/g' \
	      -e 's/B3/001.0pfu/g' \
	      -e 's/B4/000.1pfu/g'` ; 
done
