#!/bin/bash

for (( i=3; i<7; i++ ))
do
    root -b -l -q "HoughTransform_BandIter_em.C($i)"
done
