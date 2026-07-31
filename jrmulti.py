from os.path import basename
import traceback
from jrr.jrjwst import    median_combine_level3_nirspecFS,  wrap_median_combine_level3_nirspecFS

def simple_task_to_multiprocess(arg1, arg2, arg3):
    # Trying to get https://jwst-pipeline.readthedocs.io/en/1.20.0/jwst/user_documentation/running_pipeline_python.html#multiprocessing to work
    return(arg1 + 1)


def run_jwst_det1(uncal_file, output_dir, paramdict):
    """  Run the Detector1 pipeline on the given file.
    Args:
        uncal_file: str, name of uncalibrated file to run
        output_dir: str, path of the output directory"""
    from jwst.pipeline.calwebb_detector1 import Detector1Pipeline  # Goes here to prevent a memory leak

    log_name = basename(uncal_file).replace('.fits', '')

    pipe_success = False
    try:
        # Run the pipeline, turning off terminal logging messages
        Detector1Pipeline.call(uncal_file, output_dir=output_dir, steps=paramdict, save_results=True, configure_log=False)
        pipe_success = True
        print('Pipeline ran: ', uncal_file)
    except Exception:
        print('\n *** OH NO! The detector1 pipeline crashed! *** \n')
        pipe_crash_msg = traceback.print_exc()
    if not pipe_success:
        crashfile = open(log_name+'_pipecrash.txt', 'w')
        print('Printing file with full traceback')
        print(pipe_crash_msg, file=crashfile)
    

def run_jwst_spec2(infile, output_dir, paramdict):
    """  Run the Spec2 pipeline on the given file."""
    from jwst.pipeline import Spec2Pipeline  # Goes here to prevent a memory leak
    log_name = basename(infile).replace('.fits', '')
    pipe_success = False
    try:
        # Run the pipeline, turning off terminal logging messages
        Spec2Pipeline.call(infile, output_dir=output_dir, steps=paramdict, save_results=True, configure_log=False)
        pipe_success = True
        print('Pipeline ran: ', infile)
    except Exception:
        print('\n *** OH NO! The Spec2 pipeline crashed! *** \n')
        pipe_crash_msg = traceback.print_exc()
    if not pipe_success:
        crashfile = open(log_name+'_pipecrash.txt', 'w')
        print('Printing file with full traceback')
        print(pipe_crash_msg, file=crashfile)


#def run_jwst_spec3(association_file, outdir, paramdict):
#    from jwst.pipeline import Spec3Pipeline  # goes here to prevent a memory leak
#    Spec3Pipeline.call(association_file, output_dir=outdir, steps=paramdict, save_results=True, configure_log=False)
#    spec3.outlier_detection.skip = False
#    spec3.save_results = True # DON'T FORGET THIS OR YOU'LL WASTE SEVERAL HOURS FOR NAUGHT!
#    spec3.run(association_file)


    
def run_jwst_custom_extraction(indir, outdir):  # Experimental, speeding up extraction
    spectra = wrap_median_combine_level3_nirspecFS(indir, outdir)
    return(spectra)
    #df = median_combine_level3_nirspecFS(indir, thisslit, outdir)
