set -ex

PROJECT_DIR=$(pwd)/test

mkdir -p ${PROJECT_DIR}
cd ${PROJECT_DIR}
curl -LkSs https://api.github.com/repos/mirnylab/distiller-test-data/tarball \
    | tar -zxf - --strip=1 --wildcards '*/genome/*' --wildcards '*/fastq/*'
