#! /bin/bash
# run.sh <number of repetitions>
compiler="mpic++"
flags="-std=c++11 -Wall -Werror"
src="./src/main.cpp"
build="./build"
exe="$build/task"
tests_dir="./tests"
mpiexecflags=

# Create 01 test output file if it does not exist
testoutfilechecksum="8dbe67f912e0ccbc53411adfeaa42b9b"
testoutfilename="tests/01/output.txt"
testinfilename="tests/01/input.txt"
if [ ! -f $testoutfilename ]; then
  echo "Creating $testoutfilename file"
  mkdir $build 2> /dev/null
  $compiler $src -o $exe -DREFERENCE $flags
  if [ ! $? -eq 0 ]; then
    echo "[ERROR] Failed to compile $src with -DREFERENCE"
    exit 1
  fi
  mpiexec $mpiexecflags -np 1 $exe $testinfilename $testoutfilename
fi

# Assuming that md5sum is available out of the box
curchecksum=$(md5sum $testoutfilename | cut -d " " -f1)
if [ $curchecksum != $testoutfilechecksum ]; then
  echo "[WARNING] Weird, $testoutfilename checksum does not match." \
       "Did you edit test output/code inside of REFERENCE conditional" \
       "compilation block and then created test output? If you want" \
       "to fix it, try deleting $testoutfilename and running run.sh again"
fi

echo "[CLEAR]"
echo "  rm $build -r"
rm $build -rf

echo "[BUILD]"
echo "  mkdir $build"
mkdir $build
echo "  $compiler $src -o $exe $flags"
$compiler $src -o $exe $flags

if [ ! $? -eq 0 ]; then
  echo "[ERROR] Can't compile $src"
  exit 1
fi

SUCCESS_TESTS=()
FAIL_TESTS=()

for test_dir in $tests_dir/*; do
  for proc in {1..4}; do
    test="$(basename $test_dir)_p$proc"
    printf "\n[TEST $test]\n"
    echo "mpiexec -np $proc $exe $test_dir/input.txt $build/$test.txt"
    START=$(date +%s%N)
    mpiexec $mpiexecflags -np $proc $exe $test_dir/input.txt $build/$test.txt
    END=$(date +%s%N)
    DIFF=$((($END - $START)/1000000))
    if [ ! $? -eq 0 ]; then
      echo "[TEST $test] RUNTIME FAIL"
      continue;
    fi
    if cmp -s $build/$test.txt $test_dir/output.txt; then
      echo "[TEST $test] OK ($DIFF ms)"
      SUCCESS_TESTS+=($test)
    else
      echo "[TEST $test] DIFF FAIL($DIFF ms): vimdiff $build/$test.txt $test_dir/output.txt"
      FAIL_TESTS+=($test)
    fi
  done
done
echo "==========="
echo "SUCCESSFUL: ${SUCCESS_TESTS[@]}"
echo "FAIL: ${FAIL_TESTS[@]}"
