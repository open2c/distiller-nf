set -ex

USERNAME=dekkerlab
IMAGE=distiller_env

cp  ../VERSION ./VERSION
cp ../environment.yml ./environment.yml

function cleanup {
    rm  ./VERSION
    rm  ./environment.yml
}

trap cleanup EXIT

docker build -t $USERNAME/$IMAGE:latest . 
