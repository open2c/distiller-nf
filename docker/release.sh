set -ex

USERNAME=mirnylab
IMAGE=distiller_env

# ensure we're up to date
git pull

version=$(cat ../VERSION | sed -E "s/\s*version\s*=\s*'(.*)'\s*$/\1/")
echo "version: $version"

# run build
bash ./build.sh

# tag it
docker tag $USERNAME/$IMAGE:latest $USERNAME/$IMAGE:$version
# push it
docker push $USERNAME/$IMAGE:latest
docker push $USERNAME/$IMAGE:$version

