set -ex

USERNAME=mirnylab
IMAGE=distiller_env

cp  ../VERSION ./VERSION
cp ../environment.yml ./environment.yml

function cleanup {
    rm  ./VERSION
    rm  ./environment.yml
}

trap cleanup EXIT

docker build -t $USERNAME/$IMAGE:latest . 
docker run -it mirnylab/distiller_env:latest apt list > ./apt.list
docker run -it mirnylab/distiller_env:latest conda list > ./conda.list
