#!/bin/ksh
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

rm -f PET*
rm -f *.err
rm -f *.out

export XPROC=1
export YPROC=1

export NPANELS=6

NPES=$(($NPANELS * $XPROC * $YPROC ))

EXE=../bin/dynamo

echo "mpirun -np $NPES $EXE 1>$EXE.out 2>$EXE.err"

# Run the code
mpirun -np $NPES $EXE 1>$EXE.out 2>$EXE.err
