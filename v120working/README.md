This is a step by step guide to re-creating Brian's run of jwst pipeline v1.20 on leggos data.  Next is getting it to run w the new pipeline.

First, make a raw data dir:   ```mkdir /Users/jrrigby1/SCIENCE/JWST_Data/LEGGOS/Raw_data/PID04125/```       
and cd there

Second, run download_leggos_from_mast.py to download the data.  Run it as 
    ```python download_leggos_from_mast.py 04125 AWS --instrument nirspec```

The download script should move the files into directories by TARG_NAME.  However, the names in the headers are not quite what we want, especially if there were multiple IFU positions that we want to group together in the final reduction.  So, run this script to rename the Raw directory files:
      ```./rename_leggos_dirs.sh```

Next, make sure your python is pointing to a version built to work with v1.20 of the JWST pipeline, and run the this script:
```python -c "from run_pipeline_v1p21_IFS import run_pipeline; run_pipeline('SGAS1110')"```

This will call run_pipeline_v1p21_IFS.py, with is Codex's cleaned up version of Brian's v1.20 script to reduce the LEGGOS nirspec data.
The reduction will take several hours.


