#!/bin/bash

cat src/exact-solution-test.cpp | az containerapp exec -n amsc-container -g amsc --command "sh -c 'cat > /tmp/nmpde-projects/src/space-time-adaptivity/exact-solution-test.cpp'"

az containerapp exec -n amsc-container -g amsc --command \
"sh -c 'ls -l /tmp/nmpde-projects/src/space-time-adaptivity/exact-solution-test.cpp && head -n 3 /tmp/nmpde-projects/src/space-time-adaptivity/exact-solution-test.cpp'"

cat src/exact-solution-test.cpp | script -q /dev/null -c \
"az containerapp exec -n amsc-container -g amsc --command \"sh -c 'cat > /tmp/nmpde-projects/src/space-time-adaptivity/exact-solution-test.cpp'\"" \
/dev/null

az storage blob upload --account-name <storage-account-name> --container-name <container-name> --name <blob-name> --file <local-file-path>