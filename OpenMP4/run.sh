#! /bin/bash

# Set script directory for the check performed in openmp_run.inc
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Set a variable to ensure that openmp_run.inc is not ran directly
RUN=1

. $SCRIPTDIR/../openmp_run.inc
