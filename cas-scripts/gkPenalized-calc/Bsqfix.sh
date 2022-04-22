#!/bin/bash
sed -i '' 's/Bmag\[\(.\)\]\^2/Bmag[\1]*Bmag[\1]/g' $@
