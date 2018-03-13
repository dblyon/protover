# build the image
docker build ./app -t protover:latest
# mount test data directory and run interactively
docker run -it --volume /Users/dblyon/modules/protover/test_data:/data protover bash
# run container and exit
docker run --volume /Users/dblyon/modules/protover/test_data:/data protover python run_experimentfile.py
# push to dockerhub
export DOCKER_ID_USER="dblyon"
docker login
docker tag protover dblyon/protover
docker push dblyon/protover

# run the image for the first time
docker pull dblyon/protover
docker run --volume /Users/dblyon/modules/protover/test_data:/data protover python run_experimentfile.py
