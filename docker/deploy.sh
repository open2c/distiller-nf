set -ex

# push it
docker push mirnylab/distiller_env:latest
docker push mirnylab/distiller_env:$version
