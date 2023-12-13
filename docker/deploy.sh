echo "$DOCKER_PASSWORD" | docker login -u galitsyna --password-stdin

version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
echo "version: $version"

# push it
docker push galitsyna/distiller_env:latest
docker push galitsyna/distiller_env:$version
