#!/bin/sh

if [ "$#" -ge 1 -a "Z$1" = "Z--parallel" ]; then
	echo "Will compute in parallel"
	PARALLEL='1'
fi

CONFIGURATIONS="cc_nu_iso.conf cc_nubar_iso.conf nc_nu_iso.conf nc_nubar_iso.conf"

for CONFIG in $CONFIGURATIONS; do
	if [ "$PARALLEL" ]; then
		../../DISintegrate --config $CONFIG &
		sleep 1
	else
		../../DISintegrate --config $CONFIG
	fi
done

if [ "$PARALLEL" ]; then
	wait
fi
echo "Done"
