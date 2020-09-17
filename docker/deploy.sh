echo "$DOCKER_PASSWORD" | docker login -u open2c --password-stdin

version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
echo "version: $version"

# push it
docker push open2c/distiller_env:latest
docker push open2c/distiller_env:$version
