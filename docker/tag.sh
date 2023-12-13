set -ex

version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
echo "version: $version"

# tag it
docker tag galitsyna/distiller_env:minimap2 galitsyna/distiller_env:$version
