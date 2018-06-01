#!/bin/bash

MULTISCATTER=../bin/multiscatter

if [ ! "$1" ]
then
    FILE=benchmark_profile100.in
    echo Using $FILE
else
    FILE=$1
fi

echo "*** O(N) QSA model (10000 iterations)"

time $MULTISCATTER -algorithms fast none -repeat 10000 < $FILE > /dev/null

echo
echo
echo "*** O(N^2) QSA model (Hogan 2006) (1000 iterations)"

time $MULTISCATTER -algorithms original none -repeat 1000 -quiet < $FILE > /dev/null

echo
echo
echo "*** Wide angle plus O(N) QSA models (100 iterations)"

time $MULTISCATTER -algorithms fast radar -repeat 100 --quiet < $FILE > /dev/null

echo
echo
echo "*** Explicit model 2nd order (1000 iterations)" 

time $MULTISCATTER -algorithms explicit none -explicit 2 -repeat 1000 -quiet < $FILE > /dev/null

echo
echo
echo "*** Explicit model 3rd order (100 iterations)" 

time $MULTISCATTER -algorithms explicit none -explicit 3 -repeat 100 -quiet < $FILE > /dev/null

echo
echo
echo "*** Explicit model 4th order (10 iterations)" 

time $MULTISCATTER -algorithms explicit none -explicit 4 -repeat 10 -quiet < $FILE > /dev/null

echo
echo
echo "*** Explicit model 5th order (1 iteration)" 

time $MULTISCATTER -algorithms explicit none -explicit 5 -repeat 1 -quiet < $FILE > /dev/null


#[1.227 1.865 1.619 1.248 4.215 10.977]./[10000 1000 100 1000 100 10]
#[0.621 0.494 0.414 0.330 0.570 0.768]./[10000 1000 100 1000 100 10]
