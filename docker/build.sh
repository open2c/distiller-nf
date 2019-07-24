set -ex

USERNAME=mirnylab
IMAGE=distiller_env
version=$(cat ../VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')

cp  ../VERSION ./VERSION
cp ../environment.yml ./environment.yml

function cleanup {
    rm  ./VERSION
    rm  ./environment.yml
}

trap cleanup EXIT


# bop it
docker build -t $USERNAME/$IMAGE:latest . 
docker run -it mirnylab/distiller_env:latest apt list | sed 's/\x1b\[[0-9;]*m//g' | tee ./apt.list
docker run -it mirnylab/distiller_env:latest conda list | tee ./conda.list
docker images

# tag it
docker tag mirnylab/distiller_env:latest mirnylab/distiller_env:$version
