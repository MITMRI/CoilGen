#!/bin/zsh
# Run each coil of Examples/halbach_wound_gradient_set.m in its own MATLAB
# process with a hard time limit. Usage: run_jobs.sh "Z,AWG16" "Y,AWG16" ...
cd "$(dirname "$0")/.."
log=$(mktemp)
for job in "$@"; do
 WOUND_JOB=$job perl -e 'alarm shift; exec @ARGV' 300 matlab -batch "cd('Examples'); run('halbach_wound_gradient_set.m')" > $log 2>&1 < /dev/null
 st=$?
 grep -E "^DONE|^Error|Error in" $log | cut -c1-140
 [ $st -ne 0 ] && echo "FAILED $job (exit $st)"
done
echo JOBS_END
