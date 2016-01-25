#! /bin/bash

set -u # Check for unbound variables
set -o noclobber

supp_file=${SPARX}/bin/memcheck.supp

if [[ ! -f "${supp_file}" ]]; then
	# Generate suppressions
	valgrind --gen-suppressions=all 2>&1 sparx2 t_nothing | grep -v '^==' | grep -v "^t_nothing" > ${supp_file}

	echo "Generated suppresions file \`${supp_file}'"
fi

# Do memory check with suppressions=supp_file
rank=$(printenv LAMRANK)
if [[ -n ${rank} ]]; then
	OUTFILE=${1}-rank${LAMRANK}-$(hostname)-$(date '+%y%b%d_%H%M').err
	valgrind --leak-check=full --show-reachable=yes --suppressions=${supp_file} "$@" 2>&1 1>${OUTFILE}
else
	valgrind --leak-check=full --show-reachable=yes --suppressions=${supp_file} "$@"
fi
