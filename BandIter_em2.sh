#!/bin/bash

for (( i=6; i<9; i++ ))
do
    root -b -l -q "HoughTransform_BandIter_em.C($i)"
done
