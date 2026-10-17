#!/usr/bin/env bash
# Runs the single-threaded CI jobs' tests as independent processes, one per core.
#
# The single-threaded build can only keep one core busy, but its tests are independent of each
# other, so running them one after another left the rest of the runner idle for most of the job.
# This runs the same tests, split into tasks that run concurrently, longest first, each in its own
# directory under selftest-par/ so that their mlucas.cfg, fermat.cfg and logs do not collide.
#
# Usage, from the repository root after the builds:
#   bash .github/scripts/singlethread-tests.sh <Mlucas dir> <Fermat dir> <Mfactor mode suffix|->
# e.g. 'obj obj_NOSIMD _NOSIMD' on ARM, 'obj obj -' on x86. <Fermat dir> holds the Mlucas that
# config-fermat.sh runs. The Mfactor runs use obj_mfac{,_1word,...}<suffix>; '-' means no suffix.
#
# Written for bash 3.2 (the macOS /bin/bash): no wait -n, mapfile or associative arrays.
set -e -o pipefail

if [[ ${1:-} == --task ]]; then
	task=$2; mlucas_dir=$3; fermat_dir=$4; mfac_mode=$5
	[[ $mfac_mode == - ]] && mfac_mode=
	root=$PWD
	dir=selftest-par/$task
	mkdir -p "$dir"
	cd "$dir"
	exec >task.log 2>&1
	set -x
	case $task in
		s-m)
			ln -sf "$root/$mlucas_dir/Mlucas" Mlucas
			time ./Mlucas -s m
			;;
		s-small)
			ln -sf "$root/$mlucas_dir/Mlucas" Mlucas
			for s in tt t s; do time ./Mlucas -s "$s"; done
			;;
		prp)
			ln -sf "$root/$mlucas_dir/Mlucas" Mlucas
			for s in tt t s; do time ./Mlucas -prp -s "$s"; done
			;;
		fermat-top|fermat-rest)
			# config-fermat.sh runs ./Mlucas over Fermat numbers MIN..MAX. The top one is about as long
			# as all the others together, so it is a task of its own. Each part gets a copy of the
			# script with its own range.
			ln -sf "$root/$fermat_dir/Mlucas" Mlucas
			max=$(sed -n 's/^MAX=\([0-9]*\)$/\1/p' "$root/config-fermat.sh")
			if [[ $task == fermat-top ]]; then
				sed "s/^MIN=.*/MIN=$max/" "$root/config-fermat.sh" > config-fermat.sh
			else
				sed "s/^MAX=.*/MAX=$((max - 1))/" "$root/config-fermat.sh" > config-fermat.sh
			fi
			bash -e -o pipefail -- config-fermat.sh
			;;
		mfactor)
			# On ARM the native build is ASIMD, whose trial-factoring routines are stubs, so the
			# caller passes the NOSIMD suffix there.
			for word in '' 1word 2word 3word 4word nword; do
				(cd "$root/obj_mfac${word:+_$word}$mfac_mode" && time "./Mfactor${word:+_$word}" -m 521 -bmax 20)
			done
			;;
		*)
			echo "unknown task $task" >&2; exit 2
			;;
	esac
	exit 0
fi

mlucas_dir=${1:?usage: $0 <Mlucas dir> <Fermat dir> <Mfactor mode suffix|->}
fermat_dir=${2:?}
mfac_mode=${3:?usage: $0 <Mlucas dir> <Fermat dir> <Mfactor mode suffix|->}
ncpu=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)

rm -rf selftest-par
mkdir -p selftest-par
echo "Running the single-threaded tests as $ncpu concurrent processes"
# Longest first, so the short tasks fill in around the two long ones.
rc=0
printf '%s\n' s-m fermat-top fermat-rest prp s-small mfactor |
	xargs -I{} -P "$ncpu" bash "$0" --task {} "$mlucas_dir" "$fermat_dir" "$mfac_mode" || rc=$?
for t in s-m s-small prp fermat-top fermat-rest mfactor; do
	echo "::group::$t"
	cat "selftest-par/$t/task.log" 2>/dev/null || echo "(no log)"
	echo "::endgroup::"
done
if [[ $rc -ne 0 ]]; then
	echo "::error::At least one single-threaded test task failed (xargs exit status $rc); see the task logs above."
fi
exit "$rc"
