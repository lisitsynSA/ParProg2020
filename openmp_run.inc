#! /bin/bash

# DO NOT RUN THIS FILE DIRECTLY, USE run.sh FILES IN CORRESPONDING DIRECTORIES
# Special check to ensure that
if [[ $RUN != 1 ]]; then
  exit 1
fi

# Abort when error occurs

  set -e

# Declare helper functions

  function pre_target() {
    TARGET=$1
    echo "[$TARGET]"
    trap '{ set +x; } 2>/dev/null; echo "[$TARGET] FAIL"' EXIT
  }

  function post_target() {
    printf "[$TARGET] $1\n\n"
    trap - EXIT
  }

# Process options and check for some errors

  while getopts ":fvt:" opt; do
    case ${opt} in
      f ) FORCE=1
        ;;
      v ) VERBOSE=1
        ;;
      t ) TESTS+="tests/$OPTARG "
          echo "TESTS= '$TESTS'"
        ;;
      \? ) echo "Usage: run.sh [-f] [-v] [-t TESTNAME] [iternum]"
           echo "       -f:      disable paranoid errors"
           echo "       -v:      verbose"
           echo "       -t:      only run specific test. This "
           echo "                option can be specified multiple times"
           echo "       iternum: number of times to run test"
           echo "                attempting to find a race"
           exit 0
        ;;
    esac
  done
  shift $((OPTIND-1))

  if [ -z $VERBOSE ]; then
    # Use relative paths
    CURDIR="."
  else
    # Use absolute paths
    # Note that this will force failure if current user
    # does not have appropriate permissions for parent directories
    CURDIR="`pwd`"
  fi

  # Compare directory of the bash script source to current
  # Bash script's directory is set in the including script
  MESSAGE="Script source path is not equal to current directory path"
  if [ $SCRIPTDIR != `pwd` ]; then
    if [ -z $FORCE ]; then
      echo "[ERROR] $MESSAGE"
      exit 1
    else
      echo "[FORCED][WARNING] $MESSAGE"
    fi
  fi

  # Find out if we are root
  if [[ $EUID -eq 0 ]]; then
    if [ -z $FORCE ]; then
      echo "[ERROR] Attempt to run as root"
      exit 1
    else
      echo "[FORCED][WARNING] Running as root"
    fi
  fi

  # Check if ITERNUM argument is passed, default to 1
  ITERNUM="$1"
  if  [ -z "$ITERNUM" ]; then
    ITERNUM=1
  fi
  # Check if passed argument is number
  MESSAGE="Argument \"$ITERNUM\" is not a number"
  if ! [[ "$ITERNUM" =~ ^[0-9]+$ ]]; then
    if [ -z $FORCE ]; then
      echo "[ERROR] $MESSAGE"
      exit 1
    else
      echo "[FORCED][WARNING] $MESSAGE, defaulting to 1"
    fi
  fi

# CLEAR and BUILD targets

  CC="g++"
  FLAGS="-fopenmp -std=c++11 -Wall -Wextra -Werror"
  BUILDDIR="$CURDIR/build"
  EXEC="$BUILDDIR/task"
  TESTSDIR="$CURDIR/tests"
  SRC="$CURDIR/src/main.cpp"

  pre_target "CLEAR"; set -x
  rm $BUILDDIR -rf
  { set +x; } 2>/dev/null; post_target "OK"

  pre_target "BUILD"; set -x
  mkdir $BUILDDIR
  $CC $SRC -o $EXEC $FLAGS
  { set +x; } 2>/dev/null; post_target "OK"

# TESTING

  SUCCESS_TESTS=()
  FAIL_TESTS=()

  if [ -z "$TESTS" ]; then
    TESTS=$TESTSDIR/*
  fi
  for CURTESTDIR in $TESTS; do
    CURTEST="$(basename $CURTESTDIR)"

    pre_target "TEST $CURTEST"

    TIMESUM=0
    TEST_FAILED=0
    for ((i=0; i < $ITERNUM; i++)); do
      printf "\rRunning iteration $(($i+1))/$ITERNUM"

      # Run test, measure execution time
      START=$(date +%s%N)
      $EXEC $CURTESTDIR/input.txt $BUILDDIR/${CURTEST}_$i.txt
      END=$(date +%s%N)
      TIMESUM=$(($TIMESUM + ($END - $START)/1000000))

      # Check output validity
      set +e
      cmp -s $BUILDDIR/${CURTEST}_$i.txt $CURTESTDIR/output.txt
      CMPRET=$?
      set -e
      if [[ $CMPRET -ne 0 ]]; then
        RESULT="Wrong output, check with 'vimdiff $BUILDDIR/${CURTEST}_$i.txt $CURTESTDIR/output.txt'"
        FAIL_TESTS+=($CURTEST)
        TEST_FAILED=1
        break
      fi
    done
    if [[ $TEST_FAILED -eq 0 ]]; then
      RESULT="OK (average $(($TIMESUM/$ITERNUM)) ms)"
      SUCCESS_TESTS+=($CURTEST)
    fi

    echo; { post_target "$RESULT"; } 2>/dev/null
  done

# Outro

  echo "==========="
  echo "SUCCESSFUL: ${SUCCESS_TESTS[@]}"
  echo "FAIL: ${FAIL_TESTS[@]}"
