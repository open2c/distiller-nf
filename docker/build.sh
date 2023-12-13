set -ex

cp  ../VERSION ./VERSION
cp ../environment.yml ./environment.yml

function cleanup {
    rm  ./VERSION
    rm  ./environment.yml
}

trap cleanup EXIT


# bop it
docker build -t galitsyna/distiller_env:minimap2 . 
docker run -it galitsyna/distiller_env:minimap2 apt list | sed 's/\x1b\[[0-9;]*m//g' > ./apt.list
docker run -it galitsyna/distiller_env:minimap2 conda list > ./conda.list
docker images
