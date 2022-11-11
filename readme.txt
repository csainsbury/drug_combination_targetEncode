libre data pipeline

in R run process_libre_data
in debian python _ process_libre_data
then xgb classifiers etc

to start remote jupyter notebook

1. log in to remote via ssh
2. cd to Documents/code
3. run notebook wihtout browser / specify port:

jupyter notebook --no-browser --port=8888

4. in another shell set up port forwarding from remote 8888 to local 8888

ssh -N -f -L localhost:8888:localhost:8888 chris@64.71.146.66

5. start browser and go to:

localhost:8888