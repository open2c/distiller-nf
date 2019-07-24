set -ex

version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
echo "version: $version"

# tag it
docker tag mirnylab/distiller_env:latest mirnylab/distiller_env:$version
