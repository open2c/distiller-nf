set -euxo pipefail

echo "$DOCKER_PASSWORD" | docker login -u mirnylab --password-stdin

cp  ../VERSION ./VERSION
cp ../environment.yml ./environment.yml

function cleanup {
    rm  ./VERSION
    rm  ./environment.yml
}

trap cleanup EXIT

# bop it
docker build -t mirnylab/distiller_env:latest .
docker run -it mirnylab/distiller_env:latest apt list | sed 's/\x1b\[[0-9;]*m//g'
docker run -it mirnylab/distiller_env:latest conda list
docker images

# tag it
version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
echo $version

docker tag mirnylab/distiller_env:latest mirnylab/distiller_env:$version

# push it
docker push mirnylab/distiller_env:latest
docker push mirnylab/distiller_env:$version
