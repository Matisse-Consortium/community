# from mat_tools import mat_autoPipeline as mat_autoPipeline
# from mat_tools import libAutoPipeline as libAutoPipeline
import libAutoPipeline as libAutoPipeline
import matplotlib
from my_mat_autoPipeline import mat_autoPipeline

matplotlib.use("Agg")
import glob
import importlib
import os
import shutil

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import getheader
from astroquery.vizier import Vizier

# EWS
# import VioletasExp as ve
# import matutil as mu
# import wutil as wu
from ..plotting.show_allred import create_obs_lst, show_allred_mosaic
from .dependencies.pk2oifits import pk_2_oifits

# print(mat_autoPipeline.__file__)

jsdc_v2 = Vizier(catalog="II/346/jsdc_v2")

# location of the python tools: /home/varga/.local/lib/python3.6/site-packages/mat_tools/
RESDIR_MAIN = "/allegro6/matisse/varga/matisse_red9/"  # N_K2N/' #matisse_red9N/' #matisse_red9N_K2N
PROGRAMDIR = "/allegro6/matisse/varga/pro/"
dir_calib = "/allegro6/matisse/rawdata/calibrations/"  # default: '' #previously it was alt_calibdir

red_opt = "tpl"  #'tpl' #'night', 'tpl'

if red_opt == "night":
    # Reduction option 1) Night reduction
    # do night reduction for the nights listed in nightdirs
    # Note: obs_lst should be empty in this case
    obs_lst = []
    nightdir_list_fname = "nightdir_lst_p114"  # "nightdir_lst_p103_p104" #"nightdir_lst_p105_p106_p107" #"nightdir_lst_p108_p109_p110" #"nightdir_lst_p111_p112_p113"
    module = importlib.import_module(nightdir_list_fname)
    nightdirs = module.nightdirs
    DATADIR = module.DATADIR
    # note: closing '/' in the paths are important! without it create_obs_lst does not work properly (os.path.basename won't return the night folder but the folder above it)
    # nightdirs = []
    # DATADIR = '/allegro6/matisse/rawdata'
    # nightdirs = [DATADIR +'/2022_missing/'] # [DATADIR +'/2024-06-04/',DATADIR +'/2024-06-05/',]
if red_opt == "tpl":
    # Reduction option 2) Template reduction
    # reduce the data specified in obs_lst (nightdirs is not taken into account)
    nightdirs = []
    DATADIR = "/allegro6/matisse/rawdata/"  #'/allegro6/matisse/varga/rawdata/Z_CMa_img/' #'/allegro6/matisse/rawdata/'
    obs_list_fname = "obs_lst_red_p108_p109_p110"  # "obs_lst_red_p111_p112_p113" #
    # "obs_lst_red_p105_p106_p107"  #"obs_lst_red_p103_p104"
    # "obs_lst_red_p114" "obs_lst_red_p103_p104" #
    module = importlib.import_module(obs_list_fname)
    obs_lst = module.obs_lst

bands = ["L"]  # ['L','N']
redo_reduction = True  # if True, execute reduction even if there are reduced files present from earlier reduction
redo_plots = redo_reduction

if "K2N" in RESDIR_MAIN:
    try_K2N_cophasing = True  # default: True, only applies to the N band
    do_reduce_only_GRA4MAT = True  # default: False
else:
    try_K2N_cophasing = False  # default: True, only applies to the N band
    do_reduce_only_GRA4MAT = False  # default: False

pipeline = "drs"  #'drs' #'drs', 'ews'
maxIter = 1  # default: 1

do_reduction = True
do_plot = True

spectralbinning_L = "5"  #'70' #'5' #DRS option (tried value = 5 and 10 - matisse_redM), default for LR: '1'
spectralbinning_N_nominal = "7"
spectralbinning_N_low_UT = "7"  # DRS option, default for LR: '7' (tried value = 49 - matisse_redM and some matisse_red6)
spectralbinning_N_low_AT = "21"
spectralbinning_N_high_UT = "49"  # (l/dl=40) it was 49 and 35 before
spectralbinning_N_high_AT = "49"  # (l/dl=15) it was 49 and 98 before

main_JSDC_path = "/allegro6/matisse/rawdata/calibrations/JSDC/jsdc_2017_03_03.fits"
alt_JSDC_path = "/allegro6/matisse/rawdata/calibrations/JSDC/cal_catalog_for_cals_not_in_jsdc_2025-01.fits"

ews_modes = ["ews"]  # try averaging N band frames for longer in EWS  !!!!!!!!!!!!!!!!!!

# incoherent processing (default)
# coherent processing: corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/

drs_modes = [
    "coherent",
    "incoherent",
]

if pipeline == "drs":
    modes = drs_modes
if pipeline == "ews":
    modes = ews_modes

paramL_list = [
    "/spectralBinning="
    + spectralbinning_L
    + "/corrFlux=TRUE/useOpdMod=FALSE/coherentAlgo=2",
    "/spectralBinning=" + spectralbinning_L,
]

paramN_list_AT = [
    "/replaceTel=3/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/spectralBinning="
    + spectralbinning_N_nominal,
    "/replaceTel=3/spectralBinning=" + spectralbinning_N_nominal,
]
paramN_list_UT = [
    "/replaceTel=0/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/spectralBinning="
    + spectralbinning_N_nominal,
    "/replaceTel=0/spectralBinning=" + spectralbinning_N_nominal,
]

# old settings
# paramL_list = [
#    '/spectralBinning='+spectralbinning_L+'/corrFlux=TRUE/useOpdMod=TRUE/coherentAlgo=2/compensate="[pb,rb,nl,if,bp,od]"', #for coherent visibilities: '/corrFlux=FALSE/useOpdMod=FALSE/coherentAlgo=2/compensate="[pb,rb,nl,if,bp,od]"/'
#    '/spectralBinning='+spectralbinning_L+'/compensate="pb,rb,nl,if,bp,od"'
#    ]

cal_txt = ""
if len(obs_lst) == 0:
    obs_lst = []
    print("Start night reduction: ", nightdirs)
    for nightdir in nightdirs:
        # chmod_cmd = 'chmod -R 777 '+nightdir
        # print(chmod_cmd)
        # os.system(chmod_cmd)
        obs_lst_temp, cal_txt_temp = create_obs_lst(
            nightdir, write_output_file=True, outfile_path=nightdir + "/fits_cat.txt"
        )
        obs_lst += obs_lst_temp
        cal_txt += cal_txt_temp

print(" ")
print(cal_txt)
N_obs = len(obs_lst)

if "L" in bands:
    do_L = True
else:
    do_L = False
if "N" in bands:
    do_N = True
else:
    do_N = False

if do_L:
    skip_L = 0
else:
    skip_L = 1
if do_N:
    skip_N = 0
else:
    skip_N = 1
skip_L_orig = skip_L
skip_N_orig = skip_N


def check_tag(dirRaw, requested_tag, tplstart, break_when_found=False):
    hdr_lst = []
    fpath_lst = []
    # listRaw = sorted(glob.glob(dirRaw+"/MATIS*.fits"))
    listRaw = sorted(glob.glob(dirRaw + "/*.fits"))
    for filename in listRaw:
        try:
            hdr = getheader(filename, 0)
            if tplstart is None or hdr["HIERARCH ESO TPL START"] == tplstart:
                tag = libAutoPipeline.matisseType(hdr)
                if tag == requested_tag:
                    hdr_lst.append(hdr)
                    fpath_lst.append(filename)
                    if break_when_found:
                        break
        except:
            print("\nWARNING: unable to get header data.")
    return fpath_lst, hdr_lst


TPL_start_lst = []
GRA4MAT_TPL_start_lst = []

for i in range(N_obs):
    TPL_START = obs_lst[i]["tpl_start"]
    NIGHT = obs_lst[i]["night"]
    TEL = obs_lst[i]["tel"]
    DIL = obs_lst[i]["diL"]
    DIN = obs_lst[i]["diN"]
    # ----------define the folders for saving data-----------------
    RAWDIR = DATADIR + "/" + NIGHT
    # CALIBDIR = '/allegro6/matisse/varga/CalibMap'
    # CALIBDIR2 = '/allegro6/matisse/drs/calibmaps/CalibMapMayJuneJuly2019'

    # ----------run the pipeline-------------------------------
    if pipeline == "drs":
        for j in range(len(modes)):
            RESDIR = (
                RESDIR_MAIN
                + "/"
                + modes[j]
                + "/"
                + NIGHT
                + "/"
                + TPL_START.replace(":", "_")
            )
            skip_L = skip_L_orig
            skip_N = skip_N_orig
            # check for GRA4MAT observations
            if do_reduce_only_GRA4MAT == True:
                GRA4MAT_found = 0
                if TPL_START in TPL_start_lst:
                    if TPL_START in GRA4MAT_TPL_start_lst:
                        GRA4MAT_found = 1
                else:
                    rawfile_lst = glob.glob(RAWDIR + "/*")
                    for rawfile in rawfile_lst:
                        try:
                            hdul = fits.open(rawfile, ignore_missing_end=True)
                            hdr = hdul[0].header
                            key = "HIERARCH ESO TPL START"
                            if key in hdr:
                                if hdr[key] not in TPL_start_lst:
                                    TPL_start_lst.append(hdr[key])
                                key = "HIERARCH ESO DPR TYPE"
                                if key in hdr:
                                    if "RMNREC" in hdr[key]:
                                        key = "HIERARCH ESO DEL FT SENSOR"
                                        if key in hdr:
                                            if "GRAVITY" in hdr[key]:
                                                GRA4MAT_TPL_start_lst.append(
                                                    hdr["HIERARCH ESO TPL START"]
                                                )
                            hdul.close()
                        except OSError as e:
                            pass
                            # print(os.path.basename(item))
                            # print(e)
                    if TPL_START in GRA4MAT_TPL_start_lst:
                        GRA4MAT_found = 1

            if do_reduce_only_GRA4MAT == True and GRA4MAT_found == 0:
                print("Skipping " + TPL_START + "(" + modes[j] + ")")
                continue

            if not os.path.exists(RESDIR):
                os.makedirs(RESDIR)

            if do_reduction:
                if do_L:
                    # RESDIR_L = glob.glob(RESDIR + '/Iter1/*HAWAII*/')
                    RESDIR_L = (
                        RESDIR
                        + "/Iter1/mat_raw_estimates."
                        + TPL_START.replace(":", "_")
                        + ".HAWAII-2RG.rb/"
                    )
                    # check if reduced files already exist
                    reduction_success_L = True
                    n_exp = 6  # TODO: first check if there are 4 (nophot) or 6 (phot) exposures
                    for i in range(n_exp):
                        if not (
                            os.path.exists(
                                RESDIR_L + "/TARGET_RAW_INT_%04d.fits" % (i + 1)
                            )
                            or os.path.exists(
                                RESDIR_L + "/CALIB_RAW_INT_%04d.fits" % (i + 1)
                            )
                        ):
                            reduction_success_L = False
                    if os.path.isdir(RESDIR_L) and redo_reduction:
                        print("Delete existing reduction folder: " + RESDIR_L)
                        try:
                            shutil.rmtree(RESDIR_L)
                        except OSError as e:
                            print(e)
                        reduction_success_L = False
                        soffiles = glob.glob(RESDIR + "/Iter1/*HAWAII*.sof*")
                        if soffiles:
                            for file in soffiles:
                                os.remove(file)
                    if reduction_success_L and not (redo_reduction):
                        skip_L = 1
                    else:
                        if os.path.isdir(RESDIR_L):
                            print("Delete existing reduction folder: " + RESDIR_L)
                            try:
                                shutil.rmtree(RESDIR_L)
                            except OSError as e:
                                print(e)
                            soffiles = glob.glob(RESDIR + "/Iter1/*HAWAII*.sof*")
                            if soffiles:
                                for file in soffiles:
                                    os.remove(file)
                if do_N:
                    # RESDIR_N = glob.glob(RESDIR + '/Iter1/*AQUARIUS*/')
                    RESDIR_N = (
                        RESDIR
                        + "/Iter1/mat_raw_estimates."
                        + TPL_START.replace(":", "_")
                        + ".AQUARIUS.rb/"
                    )
                    reduction_success_N = True
                    n_exp = 4
                    for i in range(n_exp):
                        if not (
                            os.path.exists(
                                RESDIR_N + "/TARGET_RAW_INT_%04d.fits" % (i + 1)
                            )
                            or os.path.exists(
                                RESDIR_N + "/CALIB_RAW_INT_%04d.fits" % (i + 1)
                            )
                        ):
                            reduction_success_N = False
                    if os.path.isdir(RESDIR_N) and redo_reduction:
                        print("Delete existing reduction folder: " + RESDIR_N)
                        try:
                            shutil.rmtree(RESDIR_N)
                        except OSError as e:
                            print(e)
                        soffiles = glob.glob(RESDIR + "/Iter1/*AQUARIUS*.sof*")
                        if soffiles:
                            for file in soffiles:
                                os.remove(file)
                    if reduction_success_N and not (redo_reduction):
                        skip_N = 1
                    else:
                        if os.path.isdir(RESDIR_N):
                            print("Delete existing reduction folder: " + RESDIR_N)
                            try:
                                shutil.rmtree(RESDIR_N)
                            except OSError as e:
                                print(e)
                            soffiles = glob.glob(RESDIR + "/Iter1/*AQUARIUS*.sof*")
                            if soffiles:
                                for file in soffiles:
                                    os.remove(file)

                if not (skip_L) or not (skip_N):
                    if TEL == "UTs":
                        paramN_list = paramN_list_UT
                    if TEL == "ATs":
                        paramN_list = paramN_list_AT

                    if "HIGH" in DIN:
                        paramN_list_local = []
                        if TEL == "UTs":
                            spectralbinning_N_high = spectralbinning_N_high_UT
                        if TEL == "ATs":
                            spectralbinning_N_high = spectralbinning_N_high_AT
                        for ih in range(len(paramN_list)):
                            paramN_list_local.append(
                                paramN_list[ih].replace(
                                    "spectralBinning=" + spectralbinning_N_nominal,
                                    "spectralBinning=" + spectralbinning_N_high,
                                )
                            )  # it was 49 before
                    else:
                        paramN_list_local = []
                        if TEL == "UTs":
                            spectralbinning_N_low = spectralbinning_N_low_UT
                        if TEL == "ATs":
                            spectralbinning_N_low = spectralbinning_N_low_AT
                        for ih in range(len(paramN_list)):
                            paramN_list_local.append(
                                paramN_list[ih].replace(
                                    "spectralBinning=" + spectralbinning_N_nominal,
                                    "spectralBinning=" + spectralbinning_N_low,
                                )
                            )

                    """
                    # first try to find calibration files in RAWDIR+'/calibration_files/'
                    
                    if alt_calibdir == '':
                        if os.path.exists(RAWDIR+'/calibration_files/'):
                            dir_calib = RAWDIR+'/calibration_files/'
                        else:
                            dir_calib = RAWDIR+'/'
                    else:
                        dir_calib = alt_calibdir
                    """

                    # check whether object is CAL
                    fpath_CAL_lst, hdr_lst = check_tag(
                        RAWDIR, "CALIB_RAW", TPL_START, break_when_found=True
                    )
                    jsdc_path_match = ""
                    if fpath_CAL_lst:  # if object is a calibrator
                        hdr = hdr_lst[0]
                        # check whether calibrator is in JSDC
                        match_radius = 20.0  # arcsec
                        jsdc_match = False
                        jsdc_path_match = ""
                        target_name = hdr["HIERARCH ESO OBS TARG NAME"]
                        target_ra = hdr["RA"]  # J2000
                        target_dec = hdr["DEC"]  # J2000
                        # print(target_ra,target_dec)
                        print(
                            "Check whether calibrator "
                            + target_name
                            + " is in JSDC v2:"
                        )
                        result_jsdc = jsdc_v2.query_object(
                            target_name,
                            catalog="II/346/jsdc_v2",
                            radius=match_radius * u.arcsec,
                        )
                        # print(result_jsdc)
                        if result_jsdc != []:
                            if len(result_jsdc[0]) > 0:
                                # match
                                jsdc_match = True
                                jsdc_path_match = main_JSDC_path
                                print(
                                    "Calibrator found in JSDC v2: "
                                    + result_jsdc[0]["Name"][0]
                                )
                                #    +', separation: %.2f arcsec'%(3600.0*min_sep.value))
                        else:
                            # check whether the calibrator is in the supplement catalog
                            c_cal = SkyCoord(
                                target_ra * u.deg, target_dec * u.deg, frame="icrs"
                            )
                            caldb = fits.open(alt_JSDC_path)
                            cal_name_lst = caldb[1].data["NAME"]
                            cal_ra_lst = caldb[1].data["RAJ2000"]
                            cal_dec_lst = caldb[1].data["DEJ2000"]
                            # print(cal_ra_lst[0:10], cal_dec_lst[0:10])
                            c_lst = SkyCoord(
                                ra=cal_ra_lst,
                                dec=cal_dec_lst,
                                unit=(u.hourangle, u.deg),
                                frame="icrs",
                            )
                            # search for the calibrator in the calibrator database
                            sep = c_cal.separation(c_lst)
                            min_sep_idx = np.nanargmin(sep)
                            min_sep = sep[min_sep_idx]
                            caldb.close()
                            if (
                                min_sep < match_radius * u.deg / 3600.0
                            ):  # match_radius = 20 arcsec
                                # match
                                jsdc_match = True
                                jsdc_path_match = alt_JSDC_path
                                print(
                                    "Calibrator found in the supplement catalog: "
                                    + os.path.basename(jsdc_path_match)
                                    + ": "
                                    + cal_name_lst[min_sep_idx]
                                    + ", separation: %.2f arcsec"
                                    % (3600.0 * min_sep.value)
                                )
                        """
                        if jsdc_match:
                           
                            cats_to_remove = []
                            #check for an existing JSDC catalog in dir_calib
                            fpath_JSDC_lst,hdr_lst = check_tag(dir_calib,'JSDC_CAT',None)
                            replace_JSDC_cat = True
                            if len(fpath_JSDC_lst) > 0:
                                for fpath_JSDC,hdr in zip(fpath_JSDC_lst,hdr_lst):
                                    #check if the local JSDC catalog is the same as the reference JSDC catalog
                                    hdr_match = getheader(jsdc_path_match)
                                    if 'DATE' in hdr_match and 'DATE' in hdr:
                                        if hdr_match['DATE'] == hdr['DATE']:
                                            #if the JSDC catalog in dir_calib is the same as the one in the list of JSDC catalogs, then keep it
                                            replace_JSDC_cat = False
                                    else:
                                        cats_to_remove.append(fpath_JSDC)
                                        
                            # remove unwanted JSDC cats
                            for fpath in cats_to_remove:
                                try:
                                    print('Removing unwanted JSDC cat: '+fpath)
                                    os.remove(fpath)
                                except FileNotFoundError as e:
                                    print('File does not exist.')
                            if replace_JSDC_cat:
                                # copy the wanted JSDC catalog to dir_calib
                                print('Copying JSDC cat '+os.path.basename(jsdc_path_match)+' to '+dir_calib)
                                shutil.copyfile(jsdc_path_match,os.path.join(dir_calib,os.path.basename(jsdc_path_match)))
                        else:
                            print('Calibrator not found in JSDC, reduced data will not contain TF2, thus visibility calibration (mat_cal_oifits) will fail.')
                            #make sure that the calibration folder contain a JSDC catalog
                            
                            #check for an existing JSDC catalog in dir_calib
                            fpath_JSDC_lst,hdr_lst = check_tag(dir_calib,'JSDC_CAT',None)
                            # if there is no catalog then put a copy there
                            if len(fpath_JSDC_lst) == 0:
                                print('Copying JSDC cat '+os.path.basename(main_JSDC_path)+' to '+dir_calib)
                                shutil.copyfile(main_JSDC_path,os.path.join(dir_calib,os.path.basename(main_JSDC_path)))
                        """
                    # mat_autoPipeline.mat_autoPipeline(...)
                    mat_autoPipeline(
                        dirRaw=RAWDIR,
                        dirResult=RESDIR,
                        dirCalib=dir_calib,
                        nbCore=6,
                        tplstartsel=TPL_START,
                        resol="",
                        paramL=paramL_list[j],
                        paramN=paramN_list_local[j],
                        overwrite=0,
                        maxIter=maxIter,
                        skipL=skip_L,
                        skipN=skip_N,
                        try_K2N_cophasing=try_K2N_cophasing,
                        jsdc_path=jsdc_path_match,
                    )

                # if res == 2: #if missing calibration
                # if calibration files were not found, then use general calibration directory (CALIBDIR)
                # res = mat_autoPipeline.mat_autoPipeline(dirRaw=RAWDIR, dirResult=RESDIR, dirCalib=CALIBDIR, nbCore=6, tplstartsel=TPL_START,
                #                              resol='', paramL=paramL_list[j], paramN=paramN_list[j], overwrite=0, maxIter=1,
                #                              skipL = skip_L, skipN = skip_N)

    if pipeline == "ews":
        for j in range(len(modes)):
            RESDIR = (
                RESDIR_MAIN
                + "/"
                + modes[j]
                + "/"
                + NIGHT
                + "/"
                + TPL_START.replace(":", "_")
            )
            if not os.path.exists(RESDIR):
                os.makedirs(RESDIR)
            current_dir = os.getcwd()
            os.chdir(RAWDIR)
            if not os.path.exists(RAWDIR + "/Tpl"):
                try:
                    mu.createTplDirectories(mother="Tpl")
                except FileExistsError as e:
                    print(e)

            if do_reduction:
                if do_L:
                    # create a directory for the reduced data
                    RESDIR_BAND = (
                        RESDIR
                        + "/Iter1/mat_raw_estimates."
                        + TPL_START.replace(":", "_")
                        + ".HAWAII-2RG.rb/"
                    )
                    if os.path.exists(RESDIR_BAND):
                        try:
                            shutil.rmtree(RESDIR_BAND)
                        except OSError as e:
                            print(e)
                    os.makedirs(RESDIR_BAND)
                    os.chdir(RAWDIR + "/Tpl")

                    try:
                        try:
                            rdata = ve.processTemplate(
                                TPL_START + ".tpl", band="LM", dophot=True, ws=0.05
                            )
                            if rdata is not None:
                                ofile = (
                                    RESDIR_BAND
                                    + "/"
                                    + TPL_START.replace(":", "_")
                                    + ".tpl.pk"
                                )
                                wu.msave(ofile, rdata)
                                pk_2_oifits(
                                    RESDIR_BAND
                                    + "/"
                                    + TPL_START.replace(":", "_")
                                    + ".tpl.pk",
                                    RESDIR_BAND,
                                    oifits_template=PROGRAMDIR
                                    + "/template_oifits.fits",
                                )
                        except FileNotFoundError as e:
                            print(e)
                    except IndexError as e:
                        print(e)

                if do_N:
                    # create a directory for the reduced data
                    RESDIR_BAND = (
                        RESDIR
                        + "/Iter1/mat_raw_estimates."
                        + TPL_START.replace(":", "_")
                        + ".AQUARIUS.rb/"
                    )

                    if os.path.exists(RESDIR_BAND):
                        try:
                            shutil.rmtree(RESDIR_BAND)
                        except OSError as e:
                            print(e)
                    os.makedirs(RESDIR_BAND)

                    os.chdir(RAWDIR + "/Tpl")

                    if "HIGH" in DIN:
                        # ws = 21
                        # ws = 11 #!new value (from Jul 2021)!
                        # ws = 8 #experimental! (2021 Oct)
                        ws = 0.025  # 0.05 #dlambda: 0.03-0.06 um in HIGH (um)
                    else:
                        # ws = 3
                        ws = 0.05  # try to better match the results from DRS, old value 0.1 #dlambda: 0.3 um in LOW (um)
                    try:
                        # print('ve.processTemplate('+TPL_START+'.tpl,band="'"N"'"'+',gsmooth=0.15, wsmooth = '+'%d'%ws+',dophot=True,nbessel=0,resample=2,opdcutn=3)') #wsmmoth=1: no smoothing, 3: smooth without losing spectral resolution,  wsmooth=6 small loss of resolution
                        ###rdata = ve.processTemplate(TPL_START+'.tpl',band='N',gsmooth=0.15, wsmooth = ws,dophot=True,nbessel=0,resample=2,opdcutn=3) #old #wsmmoth=1: no smoothing, 3: smooth without losing spectral resolution,  wsmooth=6 small loss of resolution

                        try:
                            rdata = ve.processTemplate(
                                TPL_START + ".tpl",
                                band="N",
                                gsmooth=0.15,
                                wsmooth=ws,
                                dophot=True,
                                opdcutn=3,
                            )  # new wsmooth = 0.025
                            ###rdata = wu.mrestore(RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk')
                            if rdata is not None:
                                ofile = (
                                    RESDIR_BAND
                                    + "/"
                                    + TPL_START.replace(":", "_")
                                    + ".tpl.pk"
                                )
                                wu.msave(ofile, rdata)

                                # print(len(rdata),list(rdata[0].keys()))
                                # pk_2_oifits(RESDIR_BAND+'/'+TPL_START.replace(':', '_')+'.tpl.pk', RESDIR_BAND,oifits_template=PROGRAMDIR+'/template_oifits.fits',version='alpha') #old
                                print("Convert pk to oifits.")
                                pk_2_oifits(
                                    RESDIR_BAND
                                    + "/"
                                    + TPL_START.replace(":", "_")
                                    + ".tpl.pk",
                                    RESDIR_BAND,
                                    oifits_template=PROGRAMDIR
                                    + "/template_oifits.fits",
                                )  # new
                                # print(rdata[0]['bcd1name'])
                                # fluxes = ve.templateAverageFlux(rdata, True)
                                # fluxes['header'] = rdata[0]['header']
                                # fluxes['bcd1name'] = ''
                                # fluxes['bcd2name'] = ''
                                # pk_2_oifits([fluxes], RESDIR_BAND,oifits_template=PROGRAMDIR+'/template_oifits.fits')
                        except FileNotFoundError as e:
                            print(e)
                    except IndexError as e:
                        print(e)
                os.chdir(current_dir)

    # ------------------------ make plots -----------------------
    if do_plot:
        for j in range(len(modes)):
            RESDIR = (
                RESDIR_MAIN
                + "/"
                + modes[j]
                + "/"
                + NIGHT
                + "/"
                + TPL_START.replace(":", "_")
            )
            if do_L:
                inputdir = (
                    RESDIR
                    + "/Iter1/mat_raw_estimates."
                    + TPL_START.replace(":", "_")
                    + ".HAWAII-2RG.rb/"
                )
                # show_allred(inputdir,outputdir=inputdir+'/plots/',verbose=False,nbProc=6)
                show_allred_mosaic(
                    inputdir,
                    outputdir=inputdir + "/plots/",
                    fn_pattern="RAW_INT",
                    fit_model=False,
                    annotate=False,
                    redo_overwrite=redo_plots,
                )
            if do_N:
                inputdir = (
                    RESDIR
                    + "/Iter1/mat_raw_estimates."
                    + TPL_START.replace(":", "_")
                    + ".AQUARIUS.rb/"
                )
                # show_allred(inputdir, outputdir=inputdir + '/plots/',fn_pattern='', verbose=False,file_type='oifits',pro_catg='TARGET_RAW_INT',nbProc=6)
                # show_allred(inputdir, outputdir=inputdir + '/plots/',fn_pattern='', verbose=False,file_type='oifits',pro_catg='CALIB_RAW_INT',nbProc=6)
                if pipeline == "ews":
                    show_allred_mosaic(
                        inputdir,
                        outputdir=inputdir + "/plots/",
                        fn_pattern="RAW_INT",
                        wl_lim=(7.3, 13.4),
                        fit_model=False,
                        annotate=False,
                        redo_overwrite=redo_plots,
                    )
                else:
                    show_allred_mosaic(
                        inputdir,
                        outputdir=inputdir + "/plots/",
                        fn_pattern="RAW_INT",
                        fit_model=False,
                        annotate=False,
                        redo_overwrite=redo_plots,
                    )

print("EXTERMINATE!")
