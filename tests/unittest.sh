#!/bin/bash

echo "***************************************************"
echo "Unit-testing script for flusi/mhd pseudospectators."
echo "***************************************************"

# this is to reduce the number of files in the repository
#tar xzf tests_data.tar.gz

# list all the test scrits you want, separated by spaces
tests=(sponge.sh dipole.sh cylinder.sh cylinder_resume.sh)

# link flusi and mhd from .. to . if this isn't already done.
if [ ! -f UP2D ]; then
    ln -s ../UP2D .
fi

numtests=0
numsuccess=0
numfail=0

# Get time as a UNIX timestamp (seconds elapsed since Jan 1, 1970 0:00 UTC)
T="$(date +%s)"

for test in ${tests[*]}
do
    numtests=$(($numtests + 1))

    logfile=${test%%.sh}.log
    rm -f $logfile
    touch $logfile

    # run the test, copy the output to the logfile
    ./$test > $logfile
    success=$?
    if [ $success == 0 ]
    then
	echo -e "OK\tRunning test: "${test}", log in: "${logfile}
	numsuccess=$(($numsuccess + 1))
    else
	echo -e "FAIL\tRunning test: "${test}", log in: "${logfile}
	numfail=$(($numfail + 1))
    fi

    # cleanup
    rm -f *.key
    rm -f *.h5
    rm -f drag_data
    rm -f *.t end IBES_iter
    rm -f runtime*.ini
    rm -f success
    rm -f deco*
done

echo
T="$(($(date +%s)-T))"
echo "Time used in tests: ${T} seconds"

echo
echo "Total number of tests: " $numtests
echo "Total number successful: " $numsuccess
echo "Total number failed: " $numfail

if [ $numfail -gt 0 ]
then
    echo "NOT ALL TESTS PASSED: DANGER! DANGER!"
fi
