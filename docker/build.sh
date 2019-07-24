set -ex

USERNAME=mirnylab
IMAGE=distiller_env

docker build -t $USERNAME/$IMAGE:latest docker
docker run -it mirnylab/distiller_env:latest apt list | sed 's/\x1b\[[0-9;]*m//g' > ./apt.list
docker run -it mirnylab/distiller_env:latest conda list > ./conda.list
