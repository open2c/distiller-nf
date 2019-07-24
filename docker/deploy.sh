set -ex

version=$(cat ./VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')
docker run -it mirnylab/distiller_env:latest apt list | sed 's/\x1b\[[0-9;]*m//g'
docker run -it mirnylab/distiller_env:latest conda list
docker images

# tag it
docker tag mirnylab/distiller_env:latest mirnylab/distiller_env:$version

# push it
docker push mirnylab/distiller_env:latest
docker push mirnylab/distiller_env:$version
