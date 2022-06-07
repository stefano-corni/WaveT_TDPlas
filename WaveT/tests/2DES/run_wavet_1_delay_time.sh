#!/bin/bash

ICOUNT=1
ISTOP=100

while [  $ICOUNT -le $ISTOP ]; do

         dir=$ICOUNT'_0'

         cd $dir
  
         echo $dir

         ./script.sh

         cd ../

         let ICOUNT=ICOUNT+1

done

