import glob
import importlib
import os
import sys
from shutil import copyfile

import numpy as np
import pandas as pd
import show_allred
from astropy.io import fits
from avg_oifits import avg_oifits, oifits_patchwork
from fluxcal import fluxcal
from .dependencies import calib_BCD2


def get_coh_short(coh):
    if coh == "coherent":
        return "coh"
    if coh == "incoherent":
        return "incoh"
    if coh == "ews":
        return "ews"


# import multiprocessing as mp

# example useage on helada:
# > python3.6 mat_calibrate_all_drs.py | tee mat_calibrate_log_N.txt

basedir = "/allegro6/matisse/varga/"
RESDIR_MAIN = (
    basedir + "/matisse_red9N/"
)  # /matisse_red9N/' #/matisse_red9N_K2N/' #8 basedir+'/matisse_red_Luna/' #'/test_red/matisse_red8/' #matisse_red7_test #matisse_red8 matisse_red6
targetdir = (
    basedir + "/targets6/"
)  #'/targets5/' #basedir+'/targets_ews3/'  #basedir+'/targets3_Luna/' #'/targets3v2/' #'/test_red/targets_test_spbinN7/' #'/targets_ews2/' #_test_spbinN7/' #/test_red/targets2_test' #targets_ews3_test targets_ews targets2
obs_list_fname = "obs_lst_cal_p108_p109_p110"  # "obs_lst_cal_p105_p106_p107" #"obs_lst_cal_p103_p104"
# "obs_lst_cal_p111_p112_p113" # "obs_lst_cal_p103_p104" #"obs_lst_cal_p114"
docal_mode = True  # default: False, if True, calibrate only those observations which have a 'docal':1 entry in their obs_lst dictionaries
bands = ["N"]  # ['L','N']

redo_calibration = True  # if True, execute calibration even if there are calibrated files present from earlier calibrations
redo_plots = redo_calibration

# targetlist_path = basedir+'/mat_target_list.txt'
# matisse_red6 -> targets2
# matisse_red7 N band  -> targets2_5
# matisse_red7 and matisse_red8  -> targets3
# matisse_redt -> targets3 T Tau S
# part of matisse_red8  -> targets4
# matisse_red9  -> targets5,targets6 # with DRS ver. 2.0.2
# target_cal_match_filepath =  basedir + '/mat_target_list.xlsx'
# target_cal_match_sheetname = 'sci_cal_match'
# esorex_cmd = '/allegro6/matisse/drs/matissorex/matisse-kit-1.7.5-1/bin/esorex'
# esorex_cmd = '/allegro6/matisse/drs/matissorex/matisse-kit-1.7.8-5/bin/esorex'
esorex_cmd = "/allegro6/matisse/drs/matissorex/matisse-kit-2.0.2-1/bin/esorex"
# nbProc = 6
calstat_file = targetdir + "/calstat.dat"

cohs = ["coherent", "incoherent"]
# cohs= ['ews']

BCDs = [1, 2, 3, 4]
chopped_BCDs = [5, 6]

do_airmass_correction_L = True  # default: True
do_airmass_correction_N = True  # default: True
do_pwv_correction_L = False  # default: False
do_pwv_correction_N = True  # default: True
cumulBlock = True  # default: True

do_visibility_calibration = True
do_flux_calibration = True
do_BCD_calibration = True
do_averaging = True
do_make_final_oifits = True
do_make_viscorrflux = True
do_make_all_plots = True
do_make_final_plots = True
do_copy = True
do_copy_plots = True
bypass_cal_status_check = True

"""
do_visibility_calibration = False
do_flux_calibration = False
do_BCD_calibration = False
do_averaging = False
do_make_final_oifits = False
do_make_all_plots = False
do_make_final_plots = True
do_copy = False
do_copy_plots = True
bypass_cal_status_check = False
"""

module = importlib.import_module(obs_list_fname)
obs_lst = module.obs_lst

# read target list
# data = np.genfromtxt(targetlist_path,comments="#",delimiter=';',encoding=None,
#     dtype=[('f0', '<U32'), ('f1', '<U128'), ('f2', '<f8'), ('f3', '<f8'), ('f4', '<f8'), ('f5', '<f8'),
#     ('f6', '<U32'), ('f7', '<f8'), ('f8', '<f8'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'),
#     ('f12', '<f8'), ('f13', '<f8'), ('f14', '<f8'), ('f15', '<U128'), ('f16', '<U32')])
# target_names = data['f0']
# target_tags = data['f1']

if "ews" in cohs:
    if not ("ews" in targetdir):
        print("ERROR. Check targetdir, and rerun.")
        sys.exit()


def check_files_exist(input_file_list):
    for filepath in input_file_list:
        if not os.path.exists(filepath):
            return False
    return True


# open output file: calibration status
if os.path.exists(calstat_file):
    os.remove(calstat_file)
stfile = open(calstat_file, "w")
stfile.write("%-26s %1s %-19s " % ("# SCI Target name", "B", "TPL start"))
for i in range(len(cohs)):
    stfile.write("%-5s %-6s " % ("reduc", "RAWSCI"))
stfile.write("%-26s %-19s " % ("Calibrator name", "TPL start"))
for i in range(len(cohs)):
    stfile.write("%-5s %-6s " % ("reduc", "RAWCAL"))
for i in range(len(cohs)):
    stfile.write("%-5s %-6s %-6s " % ("reduc", "CALVIS", "CALFLX"))
stfile.write("%-6s " % ("CALFIN"))
if bands[0] == "L":
    stfile.write("%-3s %-2s " % ("NCh", "Ch"))
if bands[0] == "N":
    stfile.write("%-3s" % ("AVG"))
stfile.write("\n")

# read target-calibrator match excel table
if len(obs_lst) == 0:
    dfs = pd.read_excel(
        target_cal_match_filepath, sheet_name=target_cal_match_sheetname
    )
    n_obs = len(dfs)
else:
    n_obs = len(obs_lst)

for i in range(n_obs):
    if len(obs_lst) == 0:
        target_name_sci = dfs["target_name_sci"][i]
        night = dfs["night_sci"][i]
        obs_targ_name_sci = dfs["obs_targ_name_sci"][i]
        obs_targ_name_Lcal = dfs["obs_targ_name_Lcal"][i]
        obs_targ_name_Ncal = dfs["obs_targ_name_Ncal"][i]
        tpl_start_sci = dfs["tpl_start_sci"][i]
        tpl_start_Lcal = dfs["tpl_start_Lcal"][i]
        tpl_start_Ncal = dfs["tpl_start_Ncal"][i]
        skip_cal = dfs["skip_calibration"][i]
    else:
        target_name_sci = obs_lst[i]["sci_name"]
        night = obs_lst[i]["sci"]["night"]
        tpl_start_sci = obs_lst[i]["sci"]["tpl_start"]
        if "Lcal" in obs_lst[i]:
            tpl_start_Lcal = obs_lst[i]["Lcal"]["tpl_start"]
        else:
            tpl_start_Lcal = obs_lst[i]["cal"]["tpl_start"]
        if "Ncal" in obs_lst[i]:
            tpl_start_Ncal = obs_lst[i]["Ncal"]["tpl_start"]
        else:
            if "Lcal" in obs_lst[i]:
                tpl_start_Ncal = obs_lst[i]["Lcal"]["tpl_start"]
            else:
                tpl_start_Ncal = obs_lst[i]["cal"]["tpl_start"]
        skip_cal = 0
        # obs_targ_name_sci = dfs['obs_targ_name_sci'][i]
        # obs_targ_name_Lcal = dfs['obs_targ_name_Lcal'][i]
        # obs_targ_name_Ncal = dfs['obs_targ_name_Ncal'][i]
        if docal_mode == True:
            if "docal" not in obs_lst[i]:
                skip_cal = 4

    if skip_cal != 1 and skip_cal != 2 and skip_cal != 3 and skip_cal != 4:
        for band in bands:
            skip_target_flag = False
            missing_cal_flag = False

            stfile.write("%-26s %1s %-19s " % (target_name_sci, band, tpl_start_sci))
            print(
                "Calibrating "
                + target_name_sci
                + ", "
                + band
                + " band, TPL_START: "
                + tpl_start_sci
            )
            if band == "L":
                try:
                    tpl_starts_cal = tpl_start_Lcal.replace(":", "_").split(",")
                    # obs_targ_name_cal = obs_targ_name_Lcal
                except AttributeError as e:
                    # if cell is empty, nan value will be given. The replace will give an attribute error
                    try:
                        # try to use the N band calibrator
                        tpl_starts_cal = tpl_start_Ncal.replace(":", "_").split(",")
                        # obs_targ_name_cal = obs_targ_name_Ncal
                        print(
                            "Warning: N band calibrator is used for L band calibration."
                        )
                    except AttributeError as e:
                        # if there is no calibrator at all, move on to the next target
                        missing_cal_flag = True
            if band == "N":
                try:
                    tpl_starts_cal = tpl_start_Ncal.replace(":", "_").split(",")
                    # obs_targ_name_cal = obs_targ_name_Ncal
                except AttributeError as e:
                    # if cell is empty, nan value will be given. The replace will give an attribute error
                    try:
                        # try to use L band calibrator
                        tpl_starts_cal = tpl_start_Lcal.replace(":", "_").split(",")
                        # obs_targ_name_cal = obs_targ_name_Lcal
                        # print('Warning: L band calibrator is used for N band calibration.')
                    except AttributeError as e:
                        # if there is no calibrator at all, move on to the next target
                        missing_cal_flag = True

            obs_targ_name_cal = ""
            if missing_cal_flag:
                print(
                    "No calibrator found in "
                    + band
                    + " band for "
                    + target_name_sci
                    + "."
                )
                stfile.write("\n")
                continue

            # there may be multiple calibrator observations given
            n_cals = len(tpl_starts_cal)
            for j in range(n_cals):
                tpl_start_cal = tpl_starts_cal[j]
                # print('  Calibrator '+obs_targ_name_cal +', TPL_START: '+tpl_start_cal)
                if band == "L":
                    tpl_start_sci = tpl_start_sci.replace(":", "_")
                    tpldir_sci = (
                        "/"
                        + tpl_start_sci
                        + "/Iter1/mat_raw_estimates."
                        + tpl_start_sci
                        + ".HAWAII-2RG.rb/"
                    )
                    tpldir_cal = (
                        "/"
                        + tpl_start_cal
                        + "/Iter1/mat_raw_estimates."
                        + tpl_start_cal
                        + ".HAWAII-2RG.rb/"
                    )
                if band == "N":
                    tpl_start_sci = tpl_start_sci.replace(":", "_")
                    tpldir_sci = (
                        "/"
                        + tpl_start_sci
                        + "/Iter1/mat_raw_estimates."
                        + tpl_start_sci
                        + ".AQUARIUS.rb/"
                    )
                    tpldir_cal = (
                        "/"
                        + tpl_start_cal
                        + "/Iter1/mat_raw_estimates."
                        + tpl_start_cal
                        + ".AQUARIUS.rb/"
                    )

                cal_database_dir = basedir + "/CalibMap/"  # calib_spec_db.fits'
                if band == "L":
                    cal_database_paths = [
                        cal_database_dir + "vBoekelDatabase.fits",
                        cal_database_dir + "calib_spec_db_v10.fits",
                        cal_database_dir + "calib_spec_db_v10_supplement.fits",
                        cal_database_dir + "calib_spec_db_supplement3.fits",
                        cal_database_dir + "calib_spec_db_supplement4.fits",
                        cal_database_dir + "calib_spec_db_v10_supplement5.fits",
                        cal_database_dir + "calib_spec_db_v10_supplement6.fits",
                    ]
                if band == "N":
                    cal_database_paths = [
                        cal_database_dir + "vBoekelDatabase.fits",
                        cal_database_dir + "vBoekelDatabase.fitsold",
                        cal_database_dir + "calib_spec_db_v10.fits",
                        cal_database_dir + "calib_spec_db_v10_supplement.fits",
                        cal_database_dir + "calib_spec_db_supplement3.fits",
                        cal_database_dir + "calib_spec_db_supplement4.fits",
                        cal_database_dir + "calib_spec_db_v10_supplement5.fits",
                        cal_database_dir + "calib_spec_db_v10_supplement6.fits",
                    ]

                for coh in cohs:
                    print(coh)
                    inputdir = (
                        RESDIR_MAIN + "/" + coh + "/" + night + "/" + tpldir_sci + "/"
                    )
                    flst_sci = glob.glob(inputdir + "/TARGET_RAW_INT_*.fits")
                    target_tag = "TARGET"
                    if len(flst_sci) == 0:
                        flst_sci = glob.glob(inputdir + "/CALIB_RAW_INT_*.fits")
                        target_tag = "CALIB"
                    flst_sci = sorted(flst_sci)
                    # print(flst_sci)
                    # flst_cal = glob.glob(RESDIR_MAIN+coh+'/'+night+tpldir_cal+'/*_RAW_INT_*.fits')
                    # flst_cal = sorted(flst_cal)
                    # print(flst_sci)

                    if not flst_sci:
                        print(
                            "WARNING: No reduced oifits files found in sci directory."
                        )
                    else:
                        flst_cal = []
                        for i in range(len(flst_sci)):
                            inputfile_sci = flst_sci[i]
                            inputfile_cal = (
                                RESDIR_MAIN
                                + coh
                                + "/"
                                + night
                                + tpldir_cal
                                + "/"
                                + os.path.basename(flst_sci[i]).replace(
                                    target_tag, "CALIB"
                                )
                            )
                            flst_cal.append(inputfile_cal)
                            if os.path.exists(inputfile_cal):
                                obs_targ_name_cal = fits.getval(
                                    inputfile_cal, "HIERARCH ESO OBS TARG NAME", 0
                                )
                            print(
                                "  Calibrator "
                                + obs_targ_name_cal
                                + ", TPL_START: "
                                + tpl_start_cal
                            )

                            if do_visibility_calibration:
                                if redo_calibration:
                                    do_cal = True
                                else:
                                    # check if calibrated file already exists
                                    calpath = inputfile_sci.replace(
                                        "RAW_INT_0", "CAL_INT_0"
                                    )
                                    if os.path.exists(calpath):
                                        do_cal = False
                                        print(
                                            os.path.basename(calpath)
                                            + " already exists. Skip visibility calibration."
                                        )
                                    else:
                                        do_cal = True

                            if do_flux_calibration:
                                outputfile = (
                                    inputdir
                                    + "/TARGET_FLUXCAL_INT_"
                                    + os.path.basename(flst_sci[i])[-9:-5]
                                    + ".fits"
                                )
                                output_fig_dir = inputdir
                                if redo_calibration:
                                    do_cal = True
                                else:
                                    # check if calibrated file already exists
                                    calpath = outputfile
                                    if os.path.exists(calpath):
                                        do_cal = False
                                        print(
                                            os.path.basename(calpath)
                                            + " already exists. Skip flux calibration."
                                        )
                                    else:
                                        do_cal = True
                                if do_cal:
                                    print("Flux calibration.")
                                    if coh == "coherent":
                                        mode = "corrflux"
                                    if coh == "incoherent":
                                        mode = "flux"
                                    if coh == "ews":
                                        mode = "both"
                                    # if (coh == 'coherent' or coh == 'incoherent'):
                                    if band == "L":
                                        do_airmass_correction = do_airmass_correction_L
                                        do_pwv_correction = do_pwv_correction_L
                                    if band == "N":
                                        do_airmass_correction = do_airmass_correction_N
                                        do_pwv_correction = do_pwv_correction_N
                                    # else:
                                    #    do_airmass_correction=False
                                    #    do_pwv_correction=False
                                    print(inputfile_sci, inputfile_cal)
                                    try:
                                        fluxcal(
                                            inputfile_sci,
                                            inputfile_cal,
                                            outputfile,
                                            cal_database_paths,
                                            mode=mode,
                                            output_fig_dir=output_fig_dir,
                                            do_airmass_correction=do_airmass_correction,
                                            do_pwv_correction=do_pwv_correction,
                                        )
                                    except KeyError as e:
                                        print(
                                            "Missing files? ("
                                            + inputfile_sci
                                            + " / "
                                            + inputfile_cal
                                        )
                                        print(e)

                        if do_visibility_calibration:
                            if do_cal:
                                print("Visibility calibration.")
                                # ----------create an SOF file ----------------
                                SOF_PATH = (
                                    RESDIR_MAIN
                                    + "/"
                                    + coh
                                    + "/"
                                    + night
                                    + "/"
                                    + tpldir_sci
                                    + "/viscal.sof"
                                )
                                soffile = open(SOF_PATH, "w")
                                for inputfile_sci in flst_sci:
                                    soffile.write(inputfile_sci + " TARGET_RAW_INT\n")
                                for inputfile_cal in flst_cal:
                                    soffile.write(inputfile_cal + " CALIB_RAW_INT\n")
                                soffile.close()

                                if cumulBlock == True:  # default
                                    os.system(
                                        esorex_cmd
                                        + " --check-sof-exist=false --output-dir="
                                        + inputdir
                                        + " mat_cal_oifits "
                                        + SOF_PATH
                                    )
                                else:
                                    os.system(
                                        esorex_cmd
                                        + " --check-sof-exist=false --output-dir="
                                        + inputdir
                                        + " mat_cal_oifits --cumulBlock=FALSE "
                                        + SOF_PATH
                                    )

                # check if calibration succeeded
                if cohs[0] != "ews":
                    calib_ok = [False] * 4
                else:
                    calib_ok = [True] * 2
                k = 0
                for coh in cohs:
                    inputdir = (
                        RESDIR_MAIN + "/" + coh + "/" + night + "/" + tpldir_sci + "/"
                    )
                    if target_tag == "TARGET":
                        len_flst_raw_sci = len(
                            glob.glob(inputdir + "/TARGET_RAW_INT_0*.fits")
                        )
                    else:
                        len_flst_raw_sci = len(
                            glob.glob(inputdir + "/CALIB_RAW_INT_0*.fits")
                        )
                    len_flst_cal_sci = len(
                        glob.glob(inputdir + "/TARGET_CAL_INT_0*.fits")
                    )
                    if (
                        len_flst_cal_sci >= 4
                    ):  # and len_flst_cal_sci == len_flst_raw_sci:
                        calib_ok[k] = True
                    else:
                        print(
                            coh
                            + " visibility calibration error: %d CAL_INT files instead of %d."
                            % (len_flst_cal_sci, len_flst_raw_sci)
                        )
                    k = k + 1
                    len_flst_fluxcal_sci = len(
                        glob.glob(inputdir + "/TARGET_FLUXCAL_INT_0*.fits")
                    )
                    if (
                        len_flst_fluxcal_sci >= 4
                    ):  # and len_flst_fluxcal_sci == len_flst_raw_sci:
                        calib_ok[k] = True
                    else:
                        print(
                            coh
                            + " flux calibration error: %d CAL_INT files instead of %d."
                            % (len_flst_fluxcal_sci, len_flst_raw_sci)
                        )
                    k = k + 1
                if np.all(calib_ok):
                    # if the calibration is OK, exit the loop and continue with BCD calibration
                    break
                else:
                    # if some calibrated files missing, check if there is an alternative calibrator
                    if j < (n_cals - 1):
                        # try the next calibrator
                        pass
                    else:
                        # there are no alternative calibrators left
                        # set a flag to move on to next band (or next target)
                        skip_target_flag = True

            # check if reduced files exist
            reduction_success = True
            if band == "L":
                # TODO: first check if there are 4 (nophot) or 6 (phot) exposures
                n_exp = 6
            if band == "N":
                n_exp = 4
            for coh in cohs:
                # check SCI
                stfile.write("%-5s " % (get_coh_short(coh)))
                st_sci_raw = np.zeros(n_exp, dtype=bool)
                inputdir = RESDIR_MAIN + coh + "/" + night + tpldir_sci + "/"
                for i in range(n_exp):
                    if os.path.exists(
                        inputdir + "/TARGET_RAW_INT_%04d.fits" % (i + 1)
                    ) or os.path.exists(
                        inputdir + "/CALIB_RAW_INT_%04d.fits" % (i + 1)
                    ):
                        st_sci_raw[i] = 1
                    else:
                        reduction_success = False
                    stfile.write("%1d" % st_sci_raw[i])
                if band == "N":
                    stfile.write("  ")
                stfile.write(" ")
            stfile.write("%-26s %-19s " % (obs_targ_name_cal, tpl_start_cal))
            for coh in cohs:
                # check CAL
                stfile.write("%-5s " % (get_coh_short(coh)))
                inputdir_cal = RESDIR_MAIN + coh + "/" + night + tpldir_cal + "/"

                st_calib_raw = np.zeros(n_exp, dtype=bool)
                for i in range(n_exp):
                    if os.path.exists(
                        inputdir_cal + "/CALIB_RAW_INT_%04d.fits" % (i + 1)
                    ):
                        st_calib_raw[i] = 1
                    else:
                        reduction_success = False
                    stfile.write("%1d" % st_calib_raw[i])
                if band == "N":
                    stfile.write("  ")
                stfile.write(" ")
            # check if calibrated files exist
            for coh in cohs:
                stfile.write("%-5s " % (get_coh_short(coh)))
                inputdir = RESDIR_MAIN + coh + "/" + night + tpldir_sci + "/"

                st_sci_calvis = np.zeros(n_exp, dtype=bool)
                st_sci_calflux = np.zeros(n_exp, dtype=bool)
                for i in range(n_exp):
                    if os.path.exists(inputdir + "/TARGET_CAL_INT_%04d.fits" % (i + 1)):
                        st_sci_calvis[i] = 1
                    if os.path.exists(
                        inputdir + "/TARGET_FLUXCAL_INT_%04d.fits" % (i + 1)
                    ):
                        st_sci_calflux[i] = 1
                for i in range(n_exp):
                    stfile.write("%1d" % st_sci_calvis[i])
                if band == "N":
                    stfile.write("  ")
                stfile.write(" ")
                for i in range(n_exp):
                    stfile.write("%1d" % st_sci_calflux[i])
                if band == "N":
                    stfile.write("  ")
                stfile.write(" ")

            if bypass_cal_status_check:
                skip_target_flag = False

            if skip_target_flag:
                print(
                    "Calibration failed in "
                    + band
                    + " band for "
                    + target_name_sci
                    + "."
                )
                if band == "L":
                    stfile.write("000000   0  0\n")
                if band == "N":
                    stfile.write("0000     0\n")
                continue
            else:
                print("Calibration OK " + band + " band for " + target_name_sci + ".")

            # if calibration is ok, continue execution
            for coh in cohs:
                if cohs[0] != "ews":
                    inputdir = (
                        RESDIR_MAIN + "/" + coh + "/" + night + "/" + tpldir_sci + "/"
                    )
                    if do_BCD_calibration == True:
                        print("Do BCD calibration.")
                        # first calibrate non-chopped files
                        k_ii = 0
                        k_io = 0
                        k_oi = 0
                        k_oo = 0
                        # check BCDs
                        for k in BCDs:
                            try:
                                hdu = fits.open(
                                    inputdir + "/TARGET_CAL_INT_%04d.fits" % k,
                                    memmap=True,
                                    ignore_missing_end=True,
                                )
                            except FileNotFoundError as e:
                                print(
                                    inputdir
                                    + "/TARGET_CAL_INT_%04d.fits not found." % k
                                )
                                print(e)
                                skip_target_flag = True
                                break

                            hdr = hdu[0].header

                            bcd1 = hdr["HIERARCH ESO INS BCD1 ID"]
                            bcd2 = hdr["HIERARCH ESO INS BCD2 ID"]
                            # print(bcd1,bcd2)
                            if bcd1 == "OUT":
                                if bcd2 == "OUT":
                                    # OO
                                    k_oo = k
                                else:
                                    # OI
                                    k_oi = k
                            else:
                                if bcd2 == "OUT":
                                    # IO
                                    k_io = k
                                else:
                                    # II
                                    k_ii = k
                            hdu.close()
                        if skip_target_flag:
                            print(
                                "BCD calibration failed in "
                                + band
                                + " band for "
                                + target_name_sci
                                + "."
                            )
                        else:
                            inin = inputdir + "/TARGET_CAL_INT_%04d.fits" % k_ii
                            inout = inputdir + "/TARGET_CAL_INT_%04d.fits" % k_io
                            outin = inputdir + "/TARGET_CAL_INT_%04d.fits" % k_oi
                            outout = inputdir + "/TARGET_CAL_INT_%04d.fits" % k_oo
                            # print(inin,inout,outin,outout)
                            bcdcal_file = inputdir + "/TARGET_BCDCAL_INT.fits"
                            if band == "L":
                                try:
                                    calib_BCD2.calib_BCD(
                                        inin,
                                        inout,
                                        outin,
                                        outout,
                                        outputfile=bcdcal_file,
                                        plot=0,
                                    )
                                except ValueError as e:
                                    print(e)
                                    print("BCD calibration failed.")
                                    try:
                                        os.remove(bcdcal_file)
                                    except OSError as e:
                                        print(e)
                            if band == "N":
                                try:
                                    calib_BCD2.calib_BCD(
                                        inin,
                                        inout,
                                        outin,
                                        outout,
                                        outputfile=bcdcal_file,
                                        plot=0,
                                    )
                                except ValueError as e:
                                    print(e)
                                    print("BCD calibration failed.")
                                    try:
                                        os.remove(bcdcal_file)
                                    except OSError as e:
                                        print(e)

                        skip_target_flag = False
                        # check if chopped files exist
                        if os.path.exists(
                            inputdir + "/TARGET_CAL_INT_0005.fits"
                        ) and os.path.exists(inputdir + "/TARGET_CAL_INT_0006.fits"):
                            # do BCD calibration for chopped files
                            k_ii = 0
                            k_oo = 0
                            # check BCDs
                            for k in chopped_BCDs:
                                try:
                                    hdu = fits.open(
                                        inputdir + "/TARGET_CAL_INT_%04d.fits" % k,
                                        memmap=True,
                                        ignore_missing_end=True,
                                    )
                                except FileNotFoundError as e:
                                    print(
                                        inputdir
                                        + "/TARGET_CAL_INT_%04d.fits not found." % k
                                    )
                                    print(e)
                                    skip_target_flag = True
                                    break

                                hdr = hdu[0].header

                                bcd1 = hdr["HIERARCH ESO INS BCD1 ID"]
                                bcd2 = hdr["HIERARCH ESO INS BCD2 ID"]
                                # print(bcd1,bcd2)
                                if bcd1 == "OUT":
                                    if bcd2 == "OUT":
                                        # OO
                                        k_oo = k
                                else:
                                    if bcd2 == "IN":
                                        # II
                                        k_ii = k
                                hdu.close()
                            if skip_target_flag:
                                print(
                                    "BCD calibration failed in "
                                    + band
                                    + " band for "
                                    + target_name_sci
                                    + "."
                                )
                            else:
                                inin = inputdir + "/TARGET_CAL_INT_%04d.fits" % k_ii
                                outout = inputdir + "/TARGET_CAL_INT_%04d.fits" % k_oo
                                bcdcal_file = (
                                    inputdir + "/TARGET_CHOPPED_BCDCAL_INT.fits"
                                )
                                try:
                                    calib_BCD2.calib_BCD(
                                        inin,
                                        "",
                                        "",
                                        outout,
                                        outputfile=bcdcal_file,
                                        plot=0,
                                    )
                                except ValueError as e:
                                    print(e)
                                    print("BCD calibration failed.")
                                    try:
                                        os.remove(bcdcal_file)
                                    except OSError as e:
                                        print(e)
                if do_averaging:
                    print("Do averaging.")
                    headerval = [
                        {"key": "HIERARCH ESO INS BCD1 ID", "value": " "},
                        {"key": "HIERARCH ESO INS BCD2 ID", "value": " "},
                        {"key": "HIERARCH ESO INS BCD1 NAME", "value": " "},
                        {"key": "HIERARCH ESO INS BCD2 NAME", "value": " "},
                    ]
                    # average non chopped files
                    infile_list = [
                        inputdir + "TARGET_CAL_INT_%04d.fits" % (n + 1)
                        for n in range(4)
                    ]
                    # del infile_list[1]
                    # infile_list = sorted(infile_list)
                    outfile_path = inputdir + "TARGET_AVGCAL_INT.fits"
                    avg_oifits(infile_list, outfile_path, headerval=headerval)

                    infile_list = [
                        inputdir + "TARGET_FLUXCAL_INT_%04d.fits" % (n + 1)
                        for n in range(4)
                    ]
                    # del infile_list[1]
                    # infile_list = sorted(infile_list)
                    outfile_path = inputdir + "TARGET_AVGFLUXCAL_INT.fits"
                    avg_oifits(infile_list, outfile_path, headerval=headerval)

                    # check if chopped files exist
                    if os.path.exists(
                        inputdir + "/TARGET_CAL_INT_0005.fits"
                    ) and os.path.exists(inputdir + "/TARGET_CAL_INT_0006.fits"):
                        # average chopped files
                        infile_list = [
                            inputdir + "TARGET_CAL_INT_%04d.fits" % (n + 1)
                            for n in range(4, 6)
                        ]
                        outfile_path = inputdir + "TARGET_CHOPPED_AVGCAL_INT.fits"
                        avg_oifits(infile_list, outfile_path, headerval=headerval)

                        infile_list = [
                            inputdir + "TARGET_FLUXCAL_INT_%04d.fits" % (n + 1)
                            for n in range(4, 6)
                        ]
                        # infile_list = sorted(infile_list)
                        outfile_path = inputdir + "TARGET_CHOPPED_AVGFLUXCAL_INT.fits"
                        avg_oifits(infile_list, outfile_path, headerval=headerval)

            if band == "L":
                nochop_txt = "_NOCHOP"
            else:
                nochop_txt = ""

            if do_make_final_oifits:
                print("Make final OIFITS files.")
                headerval = [
                    {"key": "HIERARCH ESO INS BCD1 ID", "value": " "},
                    {"key": "HIERARCH ESO INS BCD2 ID", "value": " "},
                    {"key": "HIERARCH ESO INS BCD1 NAME", "value": " "},
                    {"key": "HIERARCH ESO INS BCD2 NAME", "value": " "},
                ]
                if cohs[0] == "ews":
                    inputdir = (
                        RESDIR_MAIN + "/ews" + "/" + night + "/" + tpldir_sci + "/"
                    )
                    for k in range(6):
                        cal = inputdir + "/TARGET_CAL_INT_%04d.fits" % (k + 1)
                        if os.path.exists(cal):
                            flux_cal = inputdir + "/TARGET_FLUXCAL_INT_%04d.fits" % (
                                k + 1
                            )
                            vis_cal = inputdir + "/TARGET_CAL_INT_%04d.fits" % (k + 1)
                            input_fitsfile_list = [flux_cal, vis_cal]
                            oi_types_list = [
                                ["flux", "visamp"],
                                ["visphi", "vis2", "t3"],
                            ]
                            ouput_file_path = (
                                inputdir + "/TARGET_MERGEDCAL_INT_%04d.fits" % (k + 1)
                            )
                            st = oifits_patchwork(
                                input_fitsfile_list,
                                ouput_file_path,
                                oi_types_list=oi_types_list,
                            )

                    flux_cal = (
                        inputdir + "/TARGET" + chop_st_txt + "_AVGFLUXCAL_INT.fits"
                    )
                    vis_cal = inputdir + "/TARGET" + chop_st_txt + "_AVGCAL_INT.fits"
                    input_fitsfile_list = [flux_cal, vis_cal]
                    oi_types_list = [["flux", "visamp"], ["visphi", "vis2", "t3"]]
                    ouput_file_path = (
                        inputdir + "/TARGET" + nochop_txt + "_FINALCAL_INT.fits"
                    )
                    if check_files_exist(input_fitsfile_list) == True:
                        st = oifits_patchwork(
                            input_fitsfile_list,
                            ouput_file_path,
                            oi_types_list=oi_types_list,
                            headerval=headerval,
                        )
                    else:
                        print("Missing files, oifits_patchwork() cannot be executed.")

                    if os.path.exists(
                        inputdir + "/TARGET_CAL_INT_0005.fits"
                    ) and os.path.exists(inputdir + "/TARGET_CAL_INT_0006.fits"):
                        flux_cal = inputdir + "/TARGET_CHOPPED_AVGFLUXCAL_INT.fits"
                        vis_cal = inputdir + "/TARGET_CHOPPED_AVGCAL_INT.fits"
                        input_fitsfile_list = [flux_cal, vis_cal]
                        oi_types_list = [["flux", "visamp"], ["visphi", "vis2", "t3"]]
                        ouput_file_path = inputdir + "/TARGET_CHOPPED_FINALCAL_INT.fits"
                        st = oifits_patchwork(
                            input_fitsfile_list,
                            ouput_file_path,
                            oi_types_list=oi_types_list,
                            headerval=headerval,
                        )

                else:
                    # DRS processing
                    inputdir_coh = (
                        RESDIR_MAIN + "/coherent" + "/" + night + "/" + tpldir_sci + "/"
                    )
                    inputdir_incoh = (
                        RESDIR_MAIN
                        + "/incoherent"
                        + "/"
                        + night
                        + "/"
                        + tpldir_sci
                        + "/"
                    )

                    # delete finalcal result files from earlier reductions:
                    cal_int_files = glob.glob(
                        inputdir_coh + "TARGET_*FINALCAL_INT.fits"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "TARGET_*FINALCAL_INT_K2N.fits"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "TARGET_*FINALCAL_INT_0*.fits"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "TARGET_*FINALCAL_INT_K2N_0*.fits"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "TARGET_MERGEDCAL_INT_0*.fits"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "TARGET_MERGEDCAL_INT_K2N_0*.fits"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "/plots/TARGET_*FINALCAL_INT_mosaic.png"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "/plots/TARGET_*FINALCAL_INT_K2N_mosaic.png"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "/plots/TARGET_*FINALCAL_INT_0*_mosaic.png"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "/plots/TARGET_*FINALCAL_INT_K2N_0*_mosaic.png"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "/plots/TARGET_MERGEDCAL_INT_0*_mosaic.png"
                    )
                    cal_int_files += glob.glob(
                        inputdir_coh + "/plots/TARGET_MERGEDCAL_INT_K2N_0*_mosaic.png"
                    )

                    for cal_int_file in cal_int_files:
                        try:
                            os.remove(cal_int_file)
                        except OSError as e:
                            print(e)

                    for k in range(6):
                        incoh_cal = inputdir_incoh + "/TARGET_CAL_INT_%04d.fits" % (
                            k + 1
                        )
                        if os.path.exists(incoh_cal):
                            coh_fluxcal = (
                                inputdir_coh + "/TARGET_FLUXCAL_INT_%04d.fits" % (k + 1)
                            )
                            incoh_fluxcal = (
                                inputdir_incoh
                                + "/TARGET_FLUXCAL_INT_%04d.fits" % (k + 1)
                            )
                            coh_cal = inputdir_coh + "/TARGET_CAL_INT_%04d.fits" % (
                                k + 1
                            )

                            k2n_tag = ""
                            if band == "L":
                                if os.path.exists(coh_cal):
                                    input_fitsfile_list = [
                                        incoh_fluxcal,
                                        coh_fluxcal,
                                        coh_cal,
                                        incoh_cal,
                                        incoh_cal,
                                    ]
                                else:
                                    input_fitsfile_list = [
                                        incoh_fluxcal,
                                        coh_fluxcal,
                                        incoh_cal,
                                        incoh_cal,
                                        incoh_cal,
                                    ]
                                oi_types_list = [
                                    ["flux"],
                                    ["visamp"],
                                    ["visphi"],
                                    ["vis2"],
                                    ["t3"],
                                ]
                            if band == "N":
                                try:
                                    hdul = fits.open(incoh_cal, ignore_missing_end=True)
                                    i = 1
                                    found = 1
                                    while found == 1:
                                        key = "HIERARCH ESO PRO REC1 CAL%d CATG" % i
                                        if key in hdul[0].header:
                                            if "RMNREC" in hdul[0].header[key]:
                                                k2n_tag = "_K2N"
                                                break
                                        else:
                                            found = 0
                                        i += 1
                                    hdul.close()
                                except OSError as e:
                                    print(os.path.basename(item))
                                    print(e)
                                if os.path.exists(coh_cal):
                                    input_fitsfile_list = [
                                        incoh_fluxcal,
                                        coh_fluxcal,
                                        coh_cal,
                                        incoh_cal,
                                        coh_cal,
                                    ]
                                else:
                                    input_fitsfile_list = [
                                        incoh_fluxcal,
                                        coh_fluxcal,
                                        incoh_cal,
                                        incoh_cal,
                                        incoh_cal,
                                    ]
                                oi_types_list = [
                                    ["flux"],
                                    ["visamp"],
                                    ["visphi"],
                                    ["visamp_to_vis2"],
                                    ["t3"],
                                ]
                            ouput_file_path = (
                                inputdir_coh
                                + "/TARGET_MERGEDCAL_INT"
                                + k2n_tag
                                + "_%04d.fits" % (k + 1)
                            )
                            st = oifits_patchwork(
                                input_fitsfile_list,
                                ouput_file_path,
                                oi_types_list=oi_types_list,
                            )

                    coh_fluxcal = inputdir_coh + "/TARGET_AVGFLUXCAL_INT.fits"
                    incoh_fluxcal = inputdir_incoh + "/TARGET_AVGFLUXCAL_INT.fits"
                    incoh_cal = inputdir_incoh + "/TARGET_AVGCAL_INT.fits"
                    coh_bcdcal = inputdir_coh + "/TARGET_BCDCAL_INT.fits"
                    incoh_bcdcal = inputdir_incoh + "/TARGET_BCDCAL_INT.fits"
                    k2n_tag = ""
                    if band == "L":
                        if os.path.exists(coh_bcdcal):
                            input_fitsfile_list = [
                                incoh_fluxcal,
                                coh_fluxcal,
                                coh_bcdcal,
                                incoh_cal,
                                incoh_bcdcal,
                            ]
                        else:
                            input_fitsfile_list = [
                                incoh_fluxcal,
                                coh_fluxcal,
                                incoh_bcdcal,
                                incoh_cal,
                                incoh_bcdcal,
                            ]
                        oi_types_list = [
                            ["flux"],
                            ["visamp"],
                            ["visphi"],
                            ["vis2"],
                            ["t3"],
                        ]
                    if band == "N":
                        try:
                            hdul = fits.open(incoh_cal, ignore_missing_end=True)
                            i = 1
                            found = 1
                            while found == 1:
                                key = "HIERARCH ESO PRO REC1 CAL%d CATG" % i
                                if key in hdul[0].header:
                                    if "RMNREC" in hdul[0].header[key]:
                                        k2n_tag = "_K2N"
                                        break
                                else:
                                    found = 0
                                i += 1
                            hdul.close()
                        except OSError as e:
                            try:
                                print(os.path.basename(item))
                                print(e)
                            except NameError as e:
                                pass
                        if os.path.exists(coh_bcdcal):
                            input_fitsfile_list = [
                                incoh_fluxcal,
                                coh_fluxcal,
                                coh_bcdcal,
                                incoh_cal,
                                coh_bcdcal,
                            ]
                        else:
                            input_fitsfile_list = [
                                incoh_fluxcal,
                                coh_fluxcal,
                                incoh_bcdcal,
                                incoh_cal,
                                incoh_bcdcal,
                            ]
                        oi_types_list = [
                            ["flux"],
                            ["visamp"],
                            ["visphi"],
                            ["visamp_to_vis2"],
                            ["t3"],
                        ]
                    ouput_file_path = (
                        inputdir_coh
                        + "/TARGET"
                        + nochop_txt
                        + "_FINALCAL_INT"
                        + k2n_tag
                        + ".fits"
                    )
                    if check_files_exist(input_fitsfile_list) == True:
                        st = oifits_patchwork(
                            input_fitsfile_list,
                            ouput_file_path,
                            oi_types_list=oi_types_list,
                            headerval=headerval,
                        )
                    else:
                        print("Missing files, oifits_patchwork() cannot be executed.")

                    if (
                        os.path.exists(inputdir_coh + "/TARGET_CAL_INT_0005.fits")
                        and os.path.exists(inputdir_coh + "/TARGET_CAL_INT_0006.fits")
                    ) or (
                        os.path.exists(inputdir_incoh + "/TARGET_CAL_INT_0005.fits")
                        and os.path.exists(inputdir_incoh + "/TARGET_CAL_INT_0006.fits")
                    ):

                        coh_fluxcal = (
                            inputdir_coh + "/TARGET_CHOPPED_AVGFLUXCAL_INT.fits"
                        )
                        incoh_fluxcal = (
                            inputdir_incoh + "/TARGET_CHOPPED_AVGFLUXCAL_INT.fits"
                        )
                        incoh_cal = inputdir_incoh + "/TARGET_CHOPPED_AVGCAL_INT.fits"
                        coh_bcdcal = inputdir_coh + "/TARGET_CHOPPED_BCDCAL_INT.fits"
                        incoh_bcdcal = (
                            inputdir_incoh + "/TARGET_CHOPPED_BCDCAL_INT.fits"
                        )
                        if os.path.exists(coh_bcdcal):
                            input_fitsfile_list = [
                                incoh_fluxcal,
                                coh_fluxcal,
                                coh_bcdcal,
                                incoh_cal,
                                incoh_bcdcal,
                            ]
                        else:
                            input_fitsfile_list = [
                                incoh_fluxcal,
                                coh_fluxcal,
                                incoh_bcdcal,
                                incoh_cal,
                                incoh_bcdcal,
                            ]
                        oi_types_list = [
                            ["flux"],
                            ["visamp"],
                            ["visphi"],
                            ["vis2"],
                            ["t3"],
                        ]
                        ouput_file_path = (
                            inputdir_coh + "/TARGET_CHOPPED_FINALCAL_INT.fits"
                        )
                        st = oifits_patchwork(
                            input_fitsfile_list,
                            ouput_file_path,
                            oi_types_list=oi_types_list,
                            headerval=headerval,
                        )

            # check if final calibrated files exist
            if cohs[0] == "ews":
                inputdir = RESDIR_MAIN + "/ews" + "/" + night + "/" + tpldir_sci + "/"
            else:
                inputdir = (
                    RESDIR_MAIN + "/coherent" + "/" + night + "/" + tpldir_sci + "/"
                )

            st_sci_mergedcal = np.zeros(n_exp, dtype=bool)
            for i in range(n_exp):
                if os.path.exists(
                    inputdir + "/TARGET_MERGEDCAL_INT_%04d.fits" % (i + 1)
                ) or os.path.exists(
                    inputdir + "/TARGET_MERGEDCAL_INT_K2N_%04d.fits" % (i + 1)
                ):
                    st_sci_mergedcal[i] = 1
                stfile.write("%1d" % st_sci_mergedcal[i])
            if band == "N":
                stfile.write("  ")
            stfile.write(" ")

            if os.path.exists(
                inputdir + "/TARGET" + nochop_txt + "_FINALCAL_INT.fits"
            ) or os.path.exists(inputdir + "/TARGET_FINALCAL_INT_K2N.fits"):
                stfile.write("  %1d " % 1)
            else:
                stfile.write("  %1d " % 0)
            if band == "L":
                if os.path.exists(inputdir + "/TARGET_CHOPPED_FINALCAL_INT.fits"):
                    stfile.write(" %1d" % 1)
                else:
                    stfile.write(" %1d" % 0)
            """
            if do_make_viscorrflux:
                flux_file_list = [
                'EX_Lup_2022-03-22T08_09_33_N_TARGET_FINALCAL_INT.fits',
                'EX_Lup_2021-03-27T06_56_05_N_TARGET_FINALCAL_INT.fits',
                ]
                for i in range(len(corr_file_list)):
                    calc_vis_from_corrflux(basedir+corr_file_list[i],
                        basedir+flux_file_list[i],
                        basedir+corr_file_list[i].replace('TARGET_FINALCAL_INT.fits','TARGET_FINALCAL_INT_VISCORRFLUX.fits'))

                show_allred_mosaic(basedir, outputdir=basedir,fn_pattern='VISCORRFLUX',fit_model = True)
            """
            fn_patterns = []
            if do_make_final_plots == True:
                fn_patterns = [
                    "TARGET" + nochop_txt + "_FINALCAL_INT",
                    "TARGET_CHOPPED_FINALCAL_INT",
                    "TARGET_MERGEDCAL_INT",
                ]
            if do_make_all_plots == True:
                fn_patterns = [
                    "TARGET" + nochop_txt + "_FINALCAL_INT",
                    "TARGET_CHOPPED_FINALCAL_INT",
                    "TARGET_MERGEDCAL_INT",
                    "TARGET_CAL_INT",
                    "TARGET_AVGCAL_INT",
                    "TARGET_CHOPPED_AVGCAL_INT",
                    "TARGET_FLUXCAL_INT",
                    "TARGET_AVGFLUXCAL_INT",
                    "TARGET_CHOPPED_AVGFLUXCAL_INT",
                    "TARGET_BCDCAL_INT",
                    "TARGET_CHOPPED_BCDCAL_INT",
                ]
            if len(fn_patterns) > 0:
                print("Make plots")
                for coh in cohs:
                    # make plots from the calibrated data
                    inputdir = (
                        RESDIR_MAIN + "/" + coh + "/" + night + "/" + tpldir_sci + "/"
                    )
                    # if len(fn_patterns) < 4:
                    for fn_pattern in fn_patterns:
                        if coh == "ews":
                            if "N" in band:
                                show_allred.show_allred_mosaic(
                                    inputdir,
                                    outputdir=inputdir + "/plots/",
                                    fn_pattern=fn_pattern,
                                    wl_lim=(7.5, 13.4),
                                    fit_model=True,
                                    annotate=False,
                                    redo_overwrite=redo_plots,
                                )
                            else:
                                show_allred.show_allred_mosaic(
                                    inputdir,
                                    outputdir=inputdir + "/plots/",
                                    fn_pattern=fn_pattern,
                                    fit_model=True,
                                    annotate=False,
                                    redo_overwrite=redo_plots,
                                )

                        else:
                            show_allred.show_allred_mosaic(
                                inputdir,
                                outputdir=inputdir + "/plots/",
                                fn_pattern=fn_pattern,
                                fit_model=True,
                                annotate=False,
                                redo_overwrite=redo_plots,
                            )  # ,sel_wl=10.5,bandwidth=0.8)
                    #     #     # show_allred.show_allred(inputdir, outputdir=inputdir + '/plots/', fn_pattern=fn_pattern,nbProc=6,fit_model=True)
                    #     #             # sel_wls=[3.53,3.58],bandwidths=[0.02,0.04])
                    #     #             # sel_wls=[3.2, 3.6, 4.0,10.7],bandwidths=[0.2, 0.2, 0.2,0.3])  # fn_pattern='TARGET*'
                    # pool = mp.Pool(processes=nbProc)
                    # pool.map(show_allred_mosaic_partial,fn_patterns)
                    # pool.terminate()
                    # pool.join()
            if cohs[0] == "ews":
                inputdir_coh = (
                    RESDIR_MAIN + "/ews" + "/" + night + "/" + tpldir_sci + "/"
                )
            else:
                inputdir_coh = (
                    RESDIR_MAIN + "/coherent" + "/" + night + "/" + tpldir_sci + "/"
                )
            dest_dir = targetdir + "/" + target_name_sci

            array = ""
            # search for calibrated files:
            cal_int_files = glob.glob(inputdir_coh + "TARGET_MERGEDCAL_INT_0*.fits")
            cal_int_files += glob.glob(
                inputdir_coh + "TARGET_MERGEDCAL_INT_K2N_0*.fits"
            )
            if len(cal_int_files) > 0:
                try:
                    hdul = fits.open(cal_int_files[0], ignore_missing_end=True)
                    hdr = hdul[0].header
                    for i in range(1, 5):
                        key = "HIERARCH ESO ISS CONF STATION%d" % i
                        if key in hdr:
                            array += hdr[key].strip()
                    hdul.close()
                except OSError as e:
                    print(os.path.basename(item))
                    print(e)

            dest_filetag = (
                target_name_sci
                + "_"
                + tpl_start_sci.replace(":", "_")
                + "_"
                + array
                + "_"
                + band
            )
            if do_copy == True:
                print("Copy final oifits files to the destination directory.")
                # copy final oifits files to the destination directory
                if not os.path.exists(dest_dir):
                    os.makedirs(dest_dir)

                # search for calibrated files:
                cal_int_files = glob.glob(inputdir_coh + "TARGET_*FINALCAL_INT.fits")
                cal_int_files += glob.glob(
                    inputdir_coh + "TARGET_*FINALCAL_INT_K2N.fits"
                )
                for cal_int_file in cal_int_files:
                    src_path = cal_int_file
                    dst_path = (
                        dest_dir
                        + "/"
                        + dest_filetag
                        + "_"
                        + os.path.basename(cal_int_file)
                    )
                    copyfile(src_path, dst_path)

                if not os.path.exists(dest_dir + "/misc/"):
                    os.makedirs(dest_dir + "/misc/")
                cal_int_files = glob.glob(inputdir_coh + "TARGET_MERGEDCAL_INT_0*.fits")
                cal_int_files += glob.glob(
                    inputdir_coh + "TARGET_MERGEDCAL_INT_K2N_0*.fits"
                )
                for cal_int_file in cal_int_files:
                    src_path = cal_int_file
                    dst_path = (
                        dest_dir
                        + "/misc/"
                        + dest_filetag
                        + "_"
                        + os.path.basename(cal_int_file)
                    )
                    copyfile(src_path, dst_path)

            # copy plots:
            if do_copy_plots:
                print("Copy final mosaic figure to the destination directory.")
                # copy final plots to the destination directory
                if not os.path.exists(dest_dir):
                    os.makedirs(dest_dir)
                cal_int_files = glob.glob(
                    inputdir_coh + "/plots/TARGET_*FINALCAL_INT_mosaic.png"
                )
                cal_int_files += glob.glob(
                    inputdir_coh + "/plots/TARGET_*FINALCAL_INT_K2N_mosaic.png"
                )
                for cal_int_file in cal_int_files:
                    src_path = cal_int_file
                    dst_path = (
                        dest_dir
                        + "/"
                        + dest_filetag
                        + "_"
                        + os.path.basename(cal_int_file)
                    )
                    copyfile(src_path, dst_path)

                if not os.path.exists(dest_dir + "/misc/"):
                    os.makedirs(dest_dir + "/misc/")
                cal_int_files = glob.glob(
                    inputdir_coh + "/plots/TARGET_MERGEDCAL_INT_0*_mosaic.png"
                )
                cal_int_files += glob.glob(
                    inputdir_coh + "/plots/TARGET_MERGEDCAL_INT_K2N_0*_mosaic.png"
                )
                for cal_int_file in cal_int_files:
                    src_path = cal_int_file
                    dst_path = (
                        dest_dir
                        + "/misc/"
                        + dest_filetag
                        + "_"
                        + os.path.basename(cal_int_file)
                    )
                    copyfile(src_path, dst_path)
            stfile.write("\n")
stfile.close()
print("EXTERMINATE!")
