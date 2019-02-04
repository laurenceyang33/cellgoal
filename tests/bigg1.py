#!/bin/python2
#============================================================
# Validate on validation/test data to test ability to 
# generalize to new conditions
# 
# Laurence Yang, UCSD
#
# 03 Oct 2018: first version
# 28 Nov 2018: applied to all BiGG models
# 20 Dec 2018: test with fully sanitized models
# 10 Jan 2019: meta_file bugfix. drop conds if model infeasible
#============================================================

from __future__ import division
from six import iteritems, string_types
from multiprocessing import Pool

from hydra.solvers import Solvers
from hydra.metrics import compute_perf
from cobra.io import load_json_model

from hydra.make_data import make_stacked_data

import hydra
import argparse
import numpy as np
import pandas as pd
import json
import itertools
import random

import warnings
import glob
import os
import time

MODEL_BASE = 'bigg_sanitized'
FBASE = 'bigg1'
OUTFILE = '%s.csv'%FBASE
SPEC_FILE = 'crossval.spec'
SOLVER_STATS_FILE = 'solver_stats.csv'
INFO_FILE = 'info.csv'
SKIP_IF_EXISTS = False

DELTA = 1e-3
MAX_ITER = 1000000
PRINT_FREQ = 1000 #MAX_ITER
VERBOSE    = 1  #0
BILINEAR_D = 1
REMOVE_BIOMASS = False  # Remove biomass completely? Already using sanitized models.
SUPPRESS_OUTPUT = False

#============================================================
# ADMM options
#============================================================
# Final convergence epsilon
ABS_CONV = 1e-9
# Over-relaxation parameter
GAM = 1.6
# ADMM penalty: initial values
RHO_LINSYS = 0.1
RHO_BND_L1 = 0.1
RHO_BI     = 0.1
ADAPT_RHO  = 2
ADAPT_FREQ = 1000
ADAPT_FOLD = 2
RES_THRESH = 1e-3

# Preconditioner: Ruiz (modified)
PRECOND    = True
# Add diagonal to linear system to ensure unique minimizer
SHIFT      = 1e-12

# Stepsize: fixed or adaptive
STEPSIZE   = 1.0      
STEP_MIN   = 0.1
STEP_MAX   = 1.0
ADAPT_STEP = 0.0
STEP_FOLD  = 2.0
STEP_THRESH= 1e-4
ADAPT_STEP_FREQ = 1000.

# Increase N_CPU if using the parallel (hydrax) solver
N_CPU      = 1
#------------------------------------------------------------
rho_set = (RHO_LINSYS, RHO_BND_L1, RHO_BI)

#============================================================
# Directory management
def pathto(*subdirs):
    p = os.path.join(hydra.__path__[0],'make_data','data', *subdirs)
    return p

class cd:
    """ Context manager: return to original dir even if exception"""
    def __init__(self, new_path):
        self.new_path = new_path

    def __enter__(self):
        self.saved_path = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.saved_path)

#============================================================

#------------------------------------------------------------
# Test missing and noisy measurements
#------------------------------------------------------------
def make_test_data(df_meas_all, df_meta_all, save_folder,
        model_file='e_coli_core.json',
        n_conds=[1], miss_fracs=[0.], nreps=1,
        need_agg=False, max_error=None, max_error_mode=None,
        elem='C', model_base=MODEL_BASE,
        skip_if_exists=False, meta_file=None):
    """
    Create test data, including missing and noisy measurements
    """
    try:
        os.mkdir(save_folder)
    except OSError:
        warnings.warn("Folder already exists: %s"%save_folder)
        if skip_if_exists:
            print("Skipping.")
            return 1

    df_meta = df_meta_all[ df_meta_all.elem==elem]

    rows_info = []
    count = 0

    model_path = os.path.join(model_base, model_file)
    #--------------------------------------------------------
    # Need to use the same random fluxes and change number of
    # conditions
    #--------------------------------------------------------
    nRxns  = df_meas_all.shape[0]
    cols_all = np.array([c for c in df_meas_all.columns if c not in ['rxn']])

    rind_dict = {(rep,miss):random.sample(range(nRxns), int(np.ceil((1-miss)*nRxns))) for
            rep in range(nreps) for miss in miss_fracs}

    for rep, miss, n_cond in itertools.product(
            range(nreps), miss_fracs, n_conds):

        fname, fext = os.path.splitext(model_file)
        out_file = '%s_%d_%d'%(fname, n_cond, count)

        #----------------------------------------------------
        # Conditions
        #----------------------------------------------------
        inds_meas = range(n_cond)
        cols_meas = cols_all[inds_meas]
        conds_include = df_meta.iloc[inds_meas]['met']

        df_meas_pert = df_meas_all[cols_meas]

        #----------------------------------------------------
        # Leave out random measurements: the same ones
        # across multiple conditions
        #----------------------------------------------------
        rinds = rind_dict[(rep,miss)]
        df_meas_pert = df_meas_pert.iloc[rinds]

        meta_path = os.path.join(model_base, meta_file)

        #----------------------------------------------------
        # Write test data
        with cd(save_folder):
            make_stacked_data.export_stacked(n_cond, model_file,
                conds_include=conds_include,
                meas_file=None, meta_file=meta_path, need_agg=need_agg,
                max_error=max_error, max_error_mode=max_error_mode,
                out_file=out_file, df_meas=df_meas_pert, verbosity=0, model_base=MODEL_BASE,
                bilinear_d=BILINEAR_D,
                exclude_biomass=True, remove_biomass=REMOVE_BIOMASS)

        count += 1
        #----------------------------------------------------
        # Write info
        for rind in rinds:
            rows_info.append({'out_file':out_file, 'miss':miss,
                'sample_ind':count, 'n_cond':n_cond, 'rep':rep, 'meas_ind':rind,
                'conds_include':','.join(conds_include)})

    mdl_name, _ = os.path.splitext(model_file)
    df_info = pd.DataFrame(rows_info)
    info_path = '%s_%s'%(mdl_name, INFO_FILE)
    df_info.to_csv(os.path.join(save_folder, info_path), index=False)

#------------------------------------------------------------
# Run tests
#------------------------------------------------------------
def run_tests(model_file, meas_file, meta_file, make_data=True,
        test_folders=[FBASE], solvers=['hydra'], model_base=MODEL_BASE,
        skip_if_exists=False):
    """
    Run all tests from files.
    """
    # All measurements:
    meas_path = pathto(meas_file)
    #--------------------------------------------------------
    try:
        df_meas_all = pd.read_csv(meas_path, index_col=0)
    except IOError:
        print("Measurement file %s not found. Skipping."%meas_path)
        return 1
    #--------------------------------------------------------

    inds_all = df_meas_all.columns

    meta_path = pathto(model_base, meta_file)
    df_meta_all = pd.read_csv(meta_path)

    # Base model
    model_path = pathto(model_base, model_file)

    solver = Solvers()

    for test_base in test_folders:

        mdl_name, _ = os.path.splitext(model_file)
        test_folder = os.path.join(test_base, mdl_name)
        #----------------------------------------------------
        # Create test data
        #----------------------------------------------------
        if make_data:
            make_test_data(df_meas_all, df_meta_all, test_folder, model_file,
                    skip_if_exists=skip_if_exists, meta_file=meta_file)

        #----------------------------------------------------
        # Load test data
        #----------------------------------------------------
        files = glob.glob(os.path.join(test_folder,'*.mtx'))
        files = [os.path.split(f)[1] for f in files]
        fbases = sorted(list(set(
            [ '_'.join(os.path.splitext(f)[0].split('_')[:-1]) 
                for f in files])))

        print("file bases: ", fbases)

        #----------------------------------------------------
        # Make results folder
        #----------------------------------------------------
        try:
            os.mkdir(os.path.join(test_folder,'results'))
        except OSError:
            pass

        print('Running %d tests in folder %s...'%(len(fbases), test_folder))

        print('%s'%('-'*170))
        print(
        '%-12.10s%-40.38s%-12.10s%-12.10s%-20.15s%-15.10s%-12.10s%-20.10s%-15.10s%-25.22s%-12.10s%-12.10s'%(
            'solver','problem','N','M','status','objval','iter','time(s)',
            'sec/iter','rho[ls,bl1,bi]','gamma','precond'))
        print('%s'%('-'*170))

        df_logs = []
        rows = []
        solver_rows = []

        for solver_id in solvers:

            sol_dict = {}
            ind_dict = {}

            for fbase in fbases:
                rho_linsys = rho_set[0]
                rho_bnd_l1 = rho_set[1]
                rho_bi = rho_set[2]

                meta_file = '%s_%s.csv'%(fbase,'meta')
                meta_path = os.path.join(test_folder, meta_file)
                df_meta = pd.read_csv(meta_path, index_col=0)
                train_inds = df_meta['cond_ind'].values
                train_inds = [str(c) for c in train_inds]
                val_inds = np.setdiff1d(inds_all, train_inds)

                ind_dict[fbase] = {'val':val_inds, 'train':train_inds}

                options = pd.Series({
                    'gam':GAM, 'delta':DELTA,
                    'max_iter':MAX_ITER, 'print_freq':PRINT_FREQ,
                    'rho_linsys':RHO_LINSYS, 'rho_bi':RHO_BI,
                    'rho_bnd_l1':RHO_BND_L1, 'precond':int(PRECOND),
                    'verbose':VERBOSE, 'abs_conv':ABS_CONV,
                    # Small diagonal to ensure full rank
                    'shift':SHIFT, 
                    # Adaptive penalty:
                    'adapt_rho':int(ADAPT_RHO), 'adapt_fold':ADAPT_FOLD, 
                    'adapt_rho_freq':ADAPT_FREQ, 'res_thresh':RES_THRESH,
                    # STEPSIZE
                    'adapt_step':ADAPT_STEP, 'adapt_step_freq':ADAPT_STEP_FREQ,
                    'step_fold':STEP_FOLD, 'step_min':STEP_MIN, 'step_max':STEP_MAX,
                    'stepsize':STEPSIZE, 'step_thresh':STEP_THRESH })

                #------------------------------------
                # Solve
                #------------------------------------
                status, result, toc = solver.solve(solver_id, test_folder, fbase, 
                        delta=DELTA, options=options, spec_file=SPEC_FILE, n_cpu=N_CPU,
                        suppress_output=SUPPRESS_OUTPUT)

                file_in  = solver.file_in
                sol_file = solver.file_out

                sol_dict[fbase] = sol_file

                m = solver.m
                n = solver.n
                objval = solver.objval
                n_iter = solver.n_iter
                try:
                    t_iter = toc / n_iter
                except:
                    print('fbase:',fbase)
                    print('toc:',toc)
                    print('n_iter:',n_iter)

                rho_str = ','.join([str(r) for r in rho_set])

                print(
            # solver   problem N     M    stat  objval  iter   time  sec/iter  rho gamma
        '%-12.10s%-40.38s%-12d%-12d%-20.17s%-15.4g%-12i%-20.4g%-15.3g%-25.22s%-12.4g%-12.10s'%(
                solver_id,file_in, int(n), int(m), status, objval, n_iter,toc,
                t_iter, rho_str, GAM, PRECOND))

                #--------------------------------------------
                # Save solver stats
                #--------------------------------------------
                solver_rows.append({'solver':solver_id, 'm':m, 'n':n,
                    'objval':objval, 'n_iter':n_iter, 't_iter':t_iter,
                    'time':toc, 'status':status, 'precond':PRECOND,
                    'rhos':rho_str, 'gamma':GAM, 'file_in':file_in})

                df_solver = pd.DataFrame(solver_rows)

                # Save now since next solve might take a long time
                solverpath = os.path.join(test_folder, SOLVER_STATS_FILE)
                df_solver.to_csv(solverpath, index=False)

            #------------------------------------------------
            # Do the cross-validation after solving all first
            #------------------------------------------------
            print('%s'%('-'*80))
            print('%-22.20s%-12.10s%-42.40s%-20.18s%-12.10s'%(
                'problem','ind','cond','perf','dtype'))
            print('%s'%('-'*80))
            for fbase, sol_file in iteritems(sol_dict):

                val_inds = ind_dict[fbase]['val']
                train_inds = ind_dict[fbase]['train']

                #--------------------------------------------
                # Read conditions
                #--------------------------------------------

                bound_file = '%s_%s.json'%(fbase,'bounds')
                bound_path = os.path.join(test_folder, bound_file)

                with open(bound_path) as f:
                    bound_dict = json.load(f)
                    
                model = load_json_model(model_path)
                #--------------------------------------------
                # Training perf
                #--------------------------------------------
                dtype = 'train'
                rows_train = []
                for cond_ind in train_inds:
                    ind = int(cond_ind)
                    df_meas = df_meas_all[cond_ind].to_frame()

                    cond = df_meta[ df_meta['cond_ind']==ind]['cond'].values[0]
                    bounds = bound_dict[cond]
                    try:
                        perf, df_sim = compute_perf(cond_ind, df_meas, model_path, sol_file, 
                                bounds=bounds, model=model)
                    except IOError:
                        perf = np.nan
                        df_sim = None

                    rows_train.append({'ind':ind, 'cond':cond,
                        'perf':perf, 'data':dtype, 'fbase':fbase})

                    # PRINT
                    print('%-22.20s%-12.10s%-42.40s%-20.4g%-12.10s'%(fbase,ind,cond,perf,dtype))

                df_perf_train = pd.DataFrame(rows_train)

                #--------------------------------------------
                # Cross-validate: on ALL other data conditions
                #--------------------------------------------
                dtype = 'validation'
                rows_val = []
                for cond_ind in val_inds:
                    ind = int(cond_ind)
                    df_meas = df_meas_all[cond_ind].to_frame()

                    try:
                        perf, df_sim = compute_perf(cond_ind, df_meas, model_path, sol_file, 
                                df_meta = df_meta_all, model=model)
                    except IOError:
                        perf = np.nan
                        df_sim = None

                    cond = df_meta_all.iloc[int(cond_ind)]['met'] 
                    rows_val.append({'ind':ind, 'cond':cond,
                        'perf':perf, 'data':dtype, 'fbase':fbase})
                    # PRINT
                    print('%-22.20s%-12.10s%-42.40s%-20.4g%-12.10s'%(fbase,ind,cond,perf,dtype))

                df_perf_val = pd.DataFrame(rows_val)

                #--------------------------------------------
                # Combine all results
                #--------------------------------------------
                df_logs.append(df_perf_train)
                df_logs.append(df_perf_val)

                df_log = pd.concat(df_logs)

                #------------------------------------
                # Record
                #------------------------------------
                outpath = os.path.join(test_folder, OUTFILE)
                df_log.to_csv(outpath,index=False)

                print('%s'%('-'*170))

def run_test_slice(model_file):
    # Create the synthetic measurements
    fname, _ = os.path.splitext(model_file)
    meas_file = '%s_ref_sims.csv'%fname
    meta_file = '%s_ref_meta.csv'%fname
    if os.path.isfile(pathto(model_base, meas_file)):
        print('model_file: %s'%model_file)
        print('meas_file: %s'%meas_file)
        print('meta_file: %s'%meta_file)
        print('Saving to: %s'%OUTFILE)
        run_tests(model_file, meas_file, meta_file, make_data=make_data,
                skip_if_exists=SKIP_IF_EXISTS)
    else:
        print('Measurement file does not exist: %s'%meas_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='meas_file, model_file')
    parser.add_argument('--model_file', metavar='model_file', type=str, default='iML1515.json') 
    parser.add_argument('--meas_file', metavar='meas_file', type=str, 
            default='iML1515_ref_sims.csv') 
    parser.add_argument('--make-data', dest='make_data', action='store_true',
            help='Make data')
    parser.add_argument('--no-make-data', dest='make_data', action='store_false',
            help='Make data')

    parser.set_defaults(make_data=True)
    args = parser.parse_args()

    model_file  = args.model_file
    meas_file   = args.meas_file
    make_data   = args.make_data

    model_base = MODEL_BASE

    p = pathto(model_base,'*.json')
    model_files = [os.path.split(g)[1] for g in glob.glob(p)]
    print('model_files:', model_files)

    #--------------------------------------------------------
    # Can run tests in parallel
    #--------------------------------------------------------
    #pool = Pool(N_CPU)
    #results = [pool.apply_async(run_test_slice, [f]) for f in model_files]
    #answer  = [result.get() for result in results]
    #--------------------------------------------------------
    # OR serially
    #--------------------------------------------------------
    for f in model_files:
        run_test_slice(f)
