#!/bin/bash

for (( i=4; i<10; i++ ))
do
    root -b -l -q "HoughTransform_BandIter.C($i)"
done
