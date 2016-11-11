#!/bin/bash

for (( i=7; i<10; i++ ))
do
    root -b -l -q "HoughTransform_BandIter.C($i)"
done
