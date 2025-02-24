# plot the contents of the fits files produced by the DRS
# created: Jozsef Varga, 2019
# varga@strw.leidenuniv.nl
#
# Example usage of the plot function (show_allred):
# from show_allred import show_allred
# datadir=r'/path/to/drs/products/'
# outputdir = datadir + 'plots/'
# show_allred(inputdir,outputdir=inputdir+'/plots/',nbProc=6)
#
# Full argument list (with default values):
# datadir: input data directory where the fits files are. The script searches for fits files, and then tries to plot them
# outputdir='plots': location for the plots (full path)
# fn_pattern = '': search pattern for filenames if you want that only specific files should be plotted
# verbose=False: whether to print many warning messages or not
# save_png=True: whether to save png figures
# save_eps=False: whether to save eps figures
# sel_wls=[np.nan]: list of wavelenghths for plotting flux/visibility against baseline
# bandwidths=[np.nan]: bandwidths corresponding to the list of wavelengths
# file_type='': if it is 'oifits', then the files in the datadir are treated like oifits files
# pro_catg='': if it is given, then only the files with the matching 'HIERARCH ESO PRO CATG' will be plotted
# nbProc=1: number of paralell processes
# annotate=True: whether to show extra info (seeing, tau0, airmass, etc.) on the plots or not
# wl_lim=(np.nan,np.nan): wavelength axis limits (um)
# B_lim=(np.nan,np.nan): baseline axis limits (m)

########################################################################

import functools
import glob
import itertools
import math
import multiprocessing as mp
import os
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits
from scipy.optimize import curve_fit

warnings.filterwarnings("ignore")


# plot style configuration
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
matplotlib.rcParams.update({"font.size": 10})


# this function sets the plot style and annotations
def plot_config(
    xlabel,
    ylabel,
    title,
    axes,
    dic,
    xlim=(None, None),
    ylim=(None, None),
    plot_legend=True,
    legend_loc="upper right",
    legend_ncol=1,
    annotate=True,
    annotate_va="top",
    annotate_fontsize=10,
    annotate_xy=(0, 1),
):
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    for ax in axes.flatten():
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        try:
            ax.set_xlim(xlim)
        except ValueError as e:
            print(e)
        try:
            ax.set_ylim(ylim)
        except ValueError as e:
            print(e)
        ax.xaxis.set_ticks_position("both")
        ax.yaxis.set_ticks_position("both")
        ax.xaxis.set_tick_params(direction="in", which="both")
        ax.yaxis.set_tick_params(direction="in", which="both")
        if plot_legend == True:
            ax.legend(
                loc=legend_loc,
                fontsize=8,
                fancybox=True,
                framealpha=0.5,
                ncol=legend_ncol,
            )
        if annotate == True:
            ax.set_title(title)
            try:
                if dic["CAL_NAME"] is not None:
                    cal_info = (
                        "\n\n"
                        + r"$\mathrm{Cal.\ name:\ %s}$" % (dic["CAL_NAME"])
                        + "\n"
                        + r"$\mathrm{Cal.\ Seeing} = %.2f''$" % (dic["CAL_SEEING"])
                        + "\n"
                        + r"$\mathrm{Cal.\ }\tau_0=%.2f$ ms"
                        % (dic["CAL_TAU0"] * 1000.0)
                        + "\n"
                        + r"$\mathrm{Cal.\ PWV} = %.2f$ mm" % (dic["CAL_PWV"])
                        + "\n"
                        + r"$\mathrm{Cal.\ Airm} = %.2f$" % (dic["CAL_AIRM"])
                        + "\n"
                        + r"$\mathrm{Cal.\ TPL\ start:}$ %s" % (dic["CAL_TPL_START"])
                        + "\n"
                        + r"$\mathrm{Cal.\ diameter} = %.2f$ mas" % (dic["CAL_DIAM"])
                        + "\n"
                        + r"$\mathrm{Cal.\ database:}$ %s" % (dic["CAL_DB"])
                    )
                else:
                    cal_info = ""

                ax.annotate(
                    r"$\alpha = %.4f^\circ$" % (dic["RA"])
                    + "\n"
                    + r"$\delta = %.4f^\circ$" % (dic["DEC"])
                    + "\n"
                    + r"$\mathrm{MJD\ obs} = %.5f$" % (dic["MJD-OBS"])
                    + "\n"
                    +
                    # r'$\mathrm{Wind sp} = %.1f$ m/s' % (dic['WINDSP']) + "\n" +
                    # r'$\mathrm{Wind dir} = %.1f^\circ$' % (dic['WINDDIR']) + "\n" +
                    r"$\mathrm{Seeing} = %.2f''$" % (dic["SEEING"])
                    + "\n"
                    + r"$\tau_0=%.2f$ ms" % (dic["TAU0"] * 1000.0)
                    + "\n"
                    +
                    # r'$T = %.1f^\circ$C' % (dic['TEMP']) + "\n" +
                    r"$\mathrm{PWV} = %.2f$ mm" % (dic["PWV"])
                    + "\n"
                    + r"$\mathrm{Airm} = %.2f$" % (dic["AIRM"])
                    + "\n"
                    + r"$\mathrm{DIT} = %.2f$ ms" % (dic["DIT"] * 1000.0)
                    + "\n"
                    + r"$\mathrm{DRS\ version:}$ %s" % (dic["PIPE_ID"])
                    + "\n"
                    + r"$\mathrm{FT\ sensor\ (status):}$ %s (%s)"
                    % (dic["FT_SENSOR"], dic["FT_STATUS"])
                    + "\n"
                    + r"$\mathrm{FT\ ratio} = %.2f$"
                    % (
                        np.mean(
                            [
                                dic["FT_RATIO1"],
                                dic["FT_RATIO2"],
                                dic["FT_RATIO3"],
                                dic["FT_RATIO4"],
                            ]
                        )
                    )
                    + cal_info,  # + "\n"
                    # r'$t_\mathrm{rel}=%.2f$ s' % (dic['TREL']),
                    xy=annotate_xy,
                    xytext=(12, -12),
                    va=annotate_va,
                    xycoords="axes fraction",
                    textcoords="offset points",
                    fontsize=annotate_fontsize,
                )
            except TypeError as e:
                print("Some header keywords missing. ")
                print(e)
        plt.tight_layout()


def check_oifits_valid(input_fitsfile):
    valid = False
    try:
        hdul = fits.open(input_fitsfile, memmap=True, ignore_missing_end=True)
        if (
            "OI_WAVELENGTH" in hdul
            and "OI_ARRAY" in hdul
            and "OI_TARGET" in hdul
            and ("OI_VIS2" in hdul or "OI_VIS" in hdul)
        ):
            valid = True
    except FileNotFoundError as e:
        print("File not found: " + fits_file)
        print(e)

    return valid


# This is a generic function to load a fits file produced by the DRS
# It supports various formats: OIFITS files, and also some DRS-specific formats
# Arguments:
# fits_file: path to the input fits file
# verbose (True or False): verbosity of output messages
# Output:
# dic: a dictionary containing the fits structures
# ext: which extension to read? OIFITS files can contain multiple extensions
def open_fits(fits_file, verbose=False, ext=0):
    dic = {}
    try:
        hdu = fits.open(fits_file, memmap=True, ignore_missing_end=True)
        hdr = hdu[0].header
        dic["IMG"] = hdu[0].data

        try:
            wl = hdu["OI_WAVELENGTH"].data["EFF_WAVE"]

            oi_wl_idx_lst = []
            # count the number of OI_WAVELENGTH tables
            for i in range(len(hdu)):
                if hdu[i].name == "OI_WAVELENGTH":
                    oi_wl_idx_lst.append(i)

            if ext >= len(oi_wl_idx_lst):
                print(
                    "Error: OI_WAVELENGTH table Ext No "
                    + "%d" % ext
                    + " does not exist."
                )
            sel_wl_idx = oi_wl_idx_lst[ext]
            wl = hdu[sel_wl_idx].data["EFF_WAVE"]
            # print(oi_wl_idx_lst)
            dic = {"WLEN": wl}
        except KeyError as e:
            if verbose:
                print("WARNING: No OI_WAVELENGTH table! ")
                print(e)

        dic["HDR"] = hdr
        dic["file"] = fits_file

        hdrkeys = [
            ["RA"],
            ["DEC"],
            ["MJD-OBS"],
            ["DATE-OBS"],
            ["HIERARCH ESO ISS AMBI TAU0 START"],
            ["HIERARCH ESO ISS AMBI FWHM START"],
            ["HIERARCH ESO ISS AMBI IWV START"],
            ["HIERARCH ESO ISS AIRM START"],
            ["HIERARCH ESO ISS AMBI TEMP"],
            ["HIERARCH ESO ISS AMBI WINDDIR"],
            ["HIERARCH ESO ISS AMBI WINDSP"],
            ["HIERARCH ESO INS BCD1 ID"],
            ["HIERARCH ESO INS BCD2 ID"],
            ["HIERARCH ESO DET NAME"],
            ["HIERARCH ESO INS DIL NAME", "HIERARCH ESO DET DISP NAME"],
            ["HIERARCH ESO INS DIN NAME", "HIERARCH ESO DET DISP NAME"],
            ["HIERARCH ESO DET SEQ1 DIT"],
            ["HIERARCH ESO PRO CATG"],
            ["HIERARCH ESO SEQ DIL WL0"],
            ["HIERARCH ESO ISS CONF STATION1"],
            ["HIERARCH ESO ISS CONF STATION2"],
            ["HIERARCH ESO ISS CONF STATION3"],
            ["HIERARCH ESO ISS CONF STATION4"],
            ["HIERARCH ESO PRO CAL DB NAME"],
            ["HIERARCH ESO PRO CAL FWHM"],
            ["HIERARCH ESO PRO CAL TAU0"],
            ["HIERARCH ESO PRO CAL IWV"],
            ["HIERARCH ESO PRO CAL AIRM"],
            ["HIERARCH ESO PRO CAL TPL START"],
            ["HIERARCH ESO PRO CAL DB DIAM"],
            ["HIERARCH ESO PRO CAL DB DBNAME"],
            ["HIERARCH ESO TPL START"],
            ["INSTRUME"],
            ["HIERARCH ESO OBS PROG ID"],
            ["HIERARCH ESO PRO REC1 PIPE ID"],
            ["HIERARCH ESO DEL FT SENSOR"],
            ["HIERARCH ESO DEL FT STATUS"],
            ["HIERARCH ESO DET TRACK MEAN1"],
            ["HIERARCH ESO DET TRACK MEAN2"],
            ["HIERARCH ESO DET TRACK MEAN3"],
            ["HIERARCH ESO DET TRACK MEAN4"],
            ["HIERARCH ESO ISS CHOP ST"],
        ]

        dickeys = [
            "RA",
            "DEC",
            "MJD-OBS",
            "DATE-OBS",
            "TAU0",
            "SEEING",
            "PWV",
            "AIRM",
            "TEMP",
            "WINDDIR",
            "WINDSP",
            "BCD1",
            "BCD2",
            "DETNAME",
            "DISPNAME_L",
            "DISPNAME_N",
            "DIT",
            "PRO_CATG",
            "WL_CENTRAL",
            "STA1",
            "STA2",
            "STA3",
            "STA4",
            "CAL_NAME",
            "CAL_SEEING",
            "CAL_TAU0",
            "CAL_PWV",
            "CAL_AIRM",
            "CAL_TPL_START",
            "CAL_DIAM",
            "CAL_DB",
            "TPL_START",
            "INSTRUMENT",
            "PROG_ID",
            "PIPE_ID",
            "FT_SENSOR",
            "FT_STATUS",
            "FT_RATIO1",
            "FT_RATIO2",
            "FT_RATIO3",
            "FT_RATIO4",
            "CHOP_ST",
        ]

        dic_default_values = [
            np.nan,
            np.nan,
            np.nan,
            "",
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            "",
            "",
            "",
            "",
            "",
            np.nan,
            "",
            np.nan,
            "",
            "",
            "",
            "",
            "",
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            "",
            np.nan,
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            "",
        ]

        for j in range(len(hdrkeys)):
            try:
                dic[dickeys[j]] = hdr[hdrkeys[j][0]]
                # print(dickeys[j],hdrkeys[j][0], hdr[hdrkeys[j][0]])
            except KeyError as e:
                if verbose:
                    print("WARNING: Header keyword " + hdrkeys[j][0] + " missing. ")
                    print(e)
                try:
                    dic[dickeys[j]] = hdr[hdrkeys[j][1]]
                except KeyError as e:
                    if verbose:
                        print("WARNING: Header keyword " + hdrkeys[j][1] + " missing. ")
                        print(e)
                    dic[dickeys[j]] = dic_default_values[j]
                except IndexError as e:
                    if verbose:
                        print("WARNING: No alternative header keyword. ")
                        print(e)
                    dic[dickeys[j]] = dic_default_values[j]

        dic["DISPNAME"] = ""
        if dic["DETNAME"]:
            if "L" in dic["DETNAME"]:
                dic["DISPNAME"] = dic["DISPNAME_L"]
            if "N" in dic["DETNAME"]:
                dic["DISPNAME"] = dic["DISPNAME_N"]
        if "PIONIER" in dic["INSTRUMENT"]:
            dic["DISPNAME"] = ""

        try:
            target_name = hdu["OI_TARGET"].data["TARGET"][0]
        except KeyError as e:
            if verbose:
                print("WARNING: No OI_TARGET table! ")
                print(e)
            target_name = ""
        if not target_name:
            try:
                target_name = hdr["HIERARCH ESO OBS TARG NAME"]
            except KeyError as e:
                if verbose:
                    print("Target name not found. ")
                    print(e)
                target_name = ""

        dic["TARGET"] = target_name

        # Fix eventual bad target identification
        # dic['TARGET'] = resolve_target(dic)

        try:
            target_category = hdu["OI_TARGET"].data["CATEGORY"][0]  # "CAL" or "SCI"
        except KeyError as e:
            if verbose:
                print("Target category not found. ")
                print(e)
            target_category = "CAL"
        dic["CATEGORY"] = target_category
        try:
            det_name = hdr["HIERARCH ESO DET CHIP NAME"]
        except KeyError as e:
            if verbose:
                print("Detector name not found. ")
                print(e)
            det_name = ""
        if det_name == "AQUARIUS":
            band = "N"
        elif det_name == "HAWAII-2RG":
            band = "LM"
        else:
            band = ""
        dic["BAND"] = band

        # non-OIFITS structures:
        # other structures:
        try:
            dic["IMG_DET"] = {}
            dic["IMG_DET"]["REGION"] = hdu["IMAGING_DETECTOR"].data["REGION"]
            dic["IMG_DET"]["REGNAME"] = hdu["IMAGING_DETECTOR"].data["REGNAME"]
        except KeyError as e:
            if verbose:
                print("WARNING: No IMAGING_DETECTOR table.")
                print(e)
        try:
            dic["IMG_DATA"] = {}
            dic["IMG_DATA"]["NREGION"] = hdu["IMAGING_DATA"].header["NREGION"]
            for i in range(1, dic["IMG_DATA"]["NREGION"] + 1):
                dic["IMG_DATA"]["DATA" + "%d" % (i)] = hdu["IMAGING_DATA"].data[
                    "DATA" + "%d" % (i)
                ]
            dic["IMG_DATA"]["MJD"] = hdu["IMAGING_DATA"].data["TIME"]
            dic["IMG_DATA"]["OPD"] = hdu["IMAGING_DATA"].data["OPD"]
            dic["IMG_DATA"]["LOCALOPD"] = hdu["IMAGING_DATA"].data["LOCALOPD"]
            dic["IMG_DATA"]["TARTYP"] = hdu["IMAGING_DATA"].data["TARTYP"]
        except KeyError as e:
            if verbose:
                print("WARNING: No IMAGING_DATA table.")
                print(e)
        except TypeError as e:
            print(e)
        try:
            dic["TEL_NAME"] = hdu["ARRAY_GEOMETRY"].data["TEL_NAME"]
            dic["STA_NAME"] = hdu["ARRAY_GEOMETRY"].data["STA_NAME"]
            dic["STA_INDEX"] = hdu["ARRAY_GEOMETRY"].data["STA_INDEX"]
        except KeyError as e:
            dic["TEL_NAME"] = {}
            dic["STA_NAME"] = {}
            dic["STA_INDEX"] = {}
            if verbose:
                print("WARNING: No ARRAY_GEOMETRY table.")
                print(e)
        # matis_eop.fits
        try:
            dic["EOP"] = {}
            dic["EOP"]["MJD"] = hdu[1].data["MJD"]
            dic["EOP"]["PMX"] = hdu[1].data["PMY"]
            dic["EOP"]["PMY"] = hdu[1].data["PMY"]
            dic["EOP"]["DUT"] = hdu[1].data["DUT"]
        except KeyError as e:
            if verbose:
                print("KeyError while attempting to read EOP table")
                print(e)
        except IndexError as e:
            if verbose:
                print("IndexError while attempting to read EOP table")
                print(e)
        except TypeError as e:
            print(e)

        try:
            dic["OPD"] = {}
            dic["OPD"]["MJD"] = hdu["OPD"].data["MJD"]
            dic["OPD"]["OPD"] = hdu["OPD"].data["OPD"]
            dic["OPD"]["OPDERR"] = hdu["OPD"].data["OPDERR"]
            dic["OPD"]["TEMPOFF"] = hdu["OPD"].data["TEMPOFF"]
            dic["OPD"]["TEMPOFFERR"] = hdu["OPD"].data["TEMPOFFERR"]
            dic["OPD"]["HUMOFF"] = hdu["OPD"].data["HUMOFF"]
            dic["OPD"]["HUMOFFERR"] = hdu["OPD"].data["HUMOFFERR"]
            dic["OPD"]["STA_INDEX"] = hdu["OPD"].data["STA_INDEX"]
        except KeyError as e:
            if verbose:
                print("WARNING: No OPD table.")
                print(e)

        # OIFITS structures
        try:
            dic["TEL_NAME"] = hdu["OI_ARRAY"].data["TEL_NAME"]
            dic["STA_NAME"] = hdu["OI_ARRAY"].data["STA_NAME"]
            dic["STA_INDEX"] = hdu["OI_ARRAY"].data["STA_INDEX"]
            dic["STAXYZ"] = hdu["OI_ARRAY"].data["STAXYZ"]
        except KeyError as e:
            dic["TEL_NAME"] = {}
            dic["STA_NAME"] = {}
            dic["STA_INDEX"] = {}
            dic["STAXYZ"] = {}
            if verbose:
                print("WARNING: No OI_ARRAY table.")
                print(e)

        oi_idx_lst = []
        dic["VIS"] = {}
        try:
            # count the number of OI_VIS tables
            for i in range(len(hdu)):
                if hdu[i].name == "OI_VIS":
                    if hdu[i].data["VISAMP"][0].size == len(wl):
                        oi_idx_lst.append(i)
            if len(oi_idx_lst) > 0:
                oi_idx_vis = oi_idx_lst[0]
                dic["VIS"]["VISAMP"] = hdu[oi_idx_vis].data["VISAMP"]
                dic["VIS"]["VISAMPERR"] = hdu[oi_idx_vis].data["VISAMPERR"]
                dic["VIS"]["VISPHI"] = hdu[oi_idx_vis].data["VISPHI"]
                dic["VIS"]["VISPHIERR"] = hdu[oi_idx_vis].data["VISPHIERR"]
                # try:
                #     dic['VIS']['CFXAMP'] = hdu[oi_idx_vis].data['CFXAMP']
                #     dic['VIS']['CFXAMPERR'] = hdu[oi_idx_vis].data['CFXAMPERR']
                # except KeyError as e:
                #     pass
                #     if verbose:
                #         print("WARNING: No correlated fluxes in this OI_VIS table!")
                #         print(e)

                dic["VIS"]["U"] = hdu[oi_idx_vis].data["UCOORD"]
                dic["VIS"]["V"] = hdu[oi_idx_vis].data["VCOORD"]
                dic["VIS"]["MJD"] = hdu[oi_idx_vis].data["MJD"]
                dic["VIS"]["STA_INDEX"] = hdu[oi_idx_vis].data["STA_INDEX"]
                dic["VIS"]["FLAG"] = hdu[oi_idx_vis].data["FLAG"]
            else:
                pass
                # ... raise an error - to be implemented
        except KeyError as e:
            if verbose:
                print("WARNING: No OI_VIS table!")
                print(e)

        if len(oi_idx_lst) > 0:
            try:
                dic["VIS"]["AMPTYP"] = hdu[oi_idx_vis].header["AMPTYP"]
            except KeyError as e:
                if verbose:
                    print("WARNING: No OI_VIS AMPTYP specified!")
                    print(e)
                dic["VIS"]["AMPTYP"] = ""
        else:
            dic["VIS"]["AMPTYP"] = ""

        oi_idx_lst = []
        try:
            # count the number of OI_VIS2 tables
            for i in range(len(hdu)):
                if hdu[i].name == "OI_VIS2":
                    if hdu[i].data["VIS2DATA"][0].size == len(wl):
                        oi_idx_lst.append(i)
            if len(oi_idx_lst) > 0:
                oi_idx_vis2 = oi_idx_lst[0]
                dic["VIS2"] = {}
                dic["VIS2"]["VIS2"] = hdu[oi_idx_vis2].data["VIS2DATA"]
                dic["VIS2"]["VIS2ERR"] = hdu[oi_idx_vis2].data["VIS2ERR"]
                dic["VIS2"]["U"] = hdu[oi_idx_vis2].data["UCOORD"]
                dic["VIS2"]["V"] = hdu[oi_idx_vis2].data["VCOORD"]
                dic["VIS2"]["MJD"] = hdu[oi_idx_vis2].data["MJD"]
                dic["VIS2"]["STA_INDEX"] = hdu[oi_idx_vis2].data["STA_INDEX"]
                dic["VIS2"]["FLAG"] = hdu[oi_idx_vis2].data["FLAG"]
        except KeyError as e:
            if verbose:
                print("WARNING: No OI_VIS2 table!")
                print(e)

        oi_idx_lst = []
        try:
            # count the number of TF2 tables
            for i in range(len(hdu)):
                if hdu[i].name == "TF2":
                    if hdu[i].data["TF2"][0].size == len(wl):
                        oi_idx_lst.append(i)
            if len(oi_idx_lst) > 0:
                oi_idx_tf2 = oi_idx_lst[0]
                dic["TF2"] = {}
                dic["TF2"]["TF2"] = hdu[oi_idx_tf2].data["TF2"]
                dic["TF2"]["TF2ERR"] = hdu[oi_idx_tf2].data["TF2ERR"]
                # dic['TF2']['U']       = hdu['OI_TF2'].data['UCOORD']
                # dic['TF2']['V']       = hdu['OI_TF2'].data['VCOORD']
                dic["TF2"]["MJD"] = hdu[oi_idx_tf2].data["MJD"]
                dic["TF2"]["STA_INDEX"] = hdu[oi_idx_tf2].data["STA_INDEX"]
                dic["TF2"]["FLAG"] = hdu[oi_idx_tf2].data["FLAG"]
        except KeyError as e:
            dic["TF2"] = {}
            dic["TF2"]["TF2"] = []
            dic["TF2"]["TF2ERR"] = []
            dic["TF2"]["MJD"] = []
            dic["TF2"]["STA_INDEX"] = []
            dic["TF2"]["FLAG"] = []
            if verbose:
                print("WARNING: No OI_TF2 table!")
                print(e)

        oi_idx_lst = []
        try:
            # count the number of OI_T3 tables
            for i in range(len(hdu)):
                if hdu[i].name == "OI_T3":
                    if hdu[i].data["T3PHI"][0].size == len(wl):
                        oi_idx_lst.append(i)
            if len(oi_idx_lst) > 0:
                oi_idx_t3 = oi_idx_lst[0]
                dic["T3"] = {}
                dic["T3"]["T3AMP"] = hdu[oi_idx_t3].data["T3AMP"]
                dic["T3"]["T3AMPERR"] = hdu[oi_idx_t3].data["T3AMPERR"]
                dic["T3"]["T3PHI"] = hdu[oi_idx_t3].data["T3PHI"]
                dic["T3"]["T3PHIERR"] = hdu[oi_idx_t3].data["T3PHIERR"]
                dic["T3"]["U1"] = hdu[oi_idx_t3].data["U1COORD"]
                dic["T3"]["V1"] = hdu[oi_idx_t3].data["V1COORD"]
                dic["T3"]["U2"] = hdu[oi_idx_t3].data["U2COORD"]
                dic["T3"]["V2"] = hdu[oi_idx_t3].data["V2COORD"]
                dic["T3"]["MJD"] = hdu[oi_idx_t3].data["MJD"]
                dic["T3"]["STA_INDEX"] = hdu[oi_idx_t3].data["STA_INDEX"]
                dic["T3"]["FLAG"] = hdu[oi_idx_t3].data["FLAG"]
        except KeyError as e:
            if verbose:
                print("WARNING: No OI_T3 table!")
                print(e)

        oi_idx_lst = []
        try:
            # count the number of OI_FLUX tables
            for i in range(len(hdu)):
                if hdu[i].name == "OI_FLUX":
                    if hdu[i].data["FLUXDATA"][0].size == len(wl):
                        oi_idx_lst.append(i)
            if len(oi_idx_lst) > 0:
                oi_idx_flux = oi_idx_lst[0]
                dic["FLUX"] = {}
                dic["FLUX"]["FLUX"] = hdu[oi_idx_flux].data["FLUXDATA"]
                dic["FLUX"]["FLUXERR"] = hdu[oi_idx_flux].data["FLUXERR"]
                dic["FLUX"]["MJD"] = hdu[oi_idx_flux].data["MJD"]
                dic["FLUX"]["STA_INDEX"] = hdu[oi_idx_flux].data["STA_INDEX"]
                dic["FLUX"]["FLAG"] = hdu[oi_idx_flux].data["FLAG"]
        except KeyError as e:
            if verbose:
                print("WARNING: No OI_FLUX table!")
                print(e)
    except FileNotFoundError as e:
        print("File not found: " + fits_file)
        print(e)
    return dic


def flip_wl(dic):
    dic_new = dic
    dic_new["WLEN"] = np.flip(dic["WLEN"])
    tables = ["VIS", "VIS2", "TF2", "T3", "FLUX"]
    columns_list = [
        ["VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR", "FLAG"],
        ["VIS2", "VIS2ERR", "FLAG"],
        ["TF2", "TF2ERR", "FLAG"],
        ["T3AMP", "T3AMPERR", "T3PHI", "T3PHIERR", "FLAG"],
        ["FLUX", "FLUXERR", "FLAG"],
    ]
    for table, columns in zip(tables, columns_list):
        if table in dic:
            for column in columns:
                if column in dic[table]:
                    for i in range(len(dic[table][column])):
                        dic_new[table][column][i] = np.flip(dic[table][column][i])
    return dic_new


# from menEWS: menEWS_utilities.py
def calculate_uv_points(inp, ha):
    #
    # calculates uv-point corresponding to inp (see "get_header_info"),
    # for hour angle(s) ha
    #

    ra, dec, BE, BN, BL, base = inp
    lat = -24.62587 * np.pi / 180.0  # Paranal lattitude in radians

    u = BE * np.cos(ha) - BN * np.sin(lat) * np.sin(ha) + BL * np.cos(lat) * np.sin(ha)
    v = (
        BE * np.sin(dec) * np.sin(ha)
        + BN * (np.sin(lat) * np.sin(dec) * np.cos(ha) + np.cos(lat) * np.cos(dec))
        - BL * (np.cos(lat) * np.sin(dec) * np.cos(ha) - np.sin(lat) * np.cos(dec))
    )
    return u, v


# B [m]
# wavelength [um]
# fwhm [mas]
def gaussian_fn(pbl, amp, fwhm, wavelength):
    x = math.pi * math.pi * 1e6 * fwhm * pbl / (1000.0 * 3600.0 * wavelength * 180.0)
    return amp * np.exp(-((x) ** 2) / (4 * np.log(2)))


# pbl [m]
# pbla [rad]
# fwhm_major [mas]
# axis ratio = b/a
# theta: gaussian major axis rotation angle (ccw) [rad]
# wavelength [um]
def twoD_gaussian_fn_ravel(xdata_tuple, amp, fwhm_major, axis_ratio, theta, wavelength):
    (pbl, pbla) = xdata_tuple
    sigma_X = (
        np.sqrt(2 * np.log(2))
        * (1000.0 * 3600.0 * wavelength * 180.0)
        / (math.pi * math.pi * 1e6 * fwhm_major)
    )
    sigma_Y = (
        np.sqrt(2 * np.log(2))
        * (1000.0 * 3600.0 * wavelength * 180.0)
        / (math.pi * math.pi * 1e6 * fwhm_major * axis_ratio)
    )
    a = np.cos(theta) ** 2 / (2 * sigma_X**2) + np.sin(theta) ** 2 / (2 * sigma_Y**2)
    b = +np.sin(2 * theta) / (4 * sigma_X**2) - np.sin(2 * theta) / (4 * sigma_Y**2)
    c = np.sin(theta) ** 2 / (2 * sigma_X**2) + np.cos(theta) ** 2 / (2 * sigma_Y**2)

    x = pbl * np.cos(pbla)
    y = pbl * np.sin(pbla)
    g = amp * np.exp(-(a * x**2 + 2 * b * x * y + c * y**2))
    return g.ravel()


# pbl [m]
# pbla [rad]
# fwhm_major [mas]
# axis ratio = b/a
# theta: gaussian major axis rotation angle (ccw) [rad]
# wavelength [um]
def twoD_gaussian_fn(xdata_tuple, amp, fwhm_major, axis_ratio, theta, wavelength):
    (pbl, pbla) = xdata_tuple
    sigma_X = (
        np.sqrt(2 * np.log(2))
        * (1000.0 * 3600.0 * wavelength * 180.0)
        / (math.pi * math.pi * 1e6 * fwhm_major)
    )
    sigma_Y = (
        np.sqrt(2 * np.log(2))
        * (1000.0 * 3600.0 * wavelength * 180.0)
        / (math.pi * math.pi * 1e6 * fwhm_major * axis_ratio)
    )
    a = np.cos(theta) ** 2 / (2 * sigma_X**2) + np.sin(theta) ** 2 / (2 * sigma_Y**2)
    b = +np.sin(2 * theta) / (4 * sigma_X**2) - np.sin(2 * theta) / (4 * sigma_Y**2)
    c = np.sin(theta) ** 2 / (2 * sigma_X**2) + np.cos(theta) ** 2 / (2 * sigma_Y**2)

    x = pbl * np.cos(pbla)
    y = pbl * np.sin(pbla)
    g = amp * np.exp(-(a * x**2 + 2 * b * x * y + c * y**2))
    return g


# mode = 'vis' or 'corr'
# B [m]
# wavelength [um]
# dist [pc]
def fit_gaussian_model(pbl, y, wavelength, dist=np.nan, mode="vis"):
    n = len(pbl)
    mean = sum(pbl * y) / n
    fwhm_init = 5.0  # [mas] #sum(y*(pbl-mean)**2)/n
    # print(np.nanmax(y),fwhm_init)
    # if mode == 'corr':
    try:
        popt, pcov = curve_fit(
            lambda pbl, amp, fwhm: gaussian_fn(pbl, amp, fwhm, wavelength),
            pbl,
            y,
            p0=[np.nanmax(y), fwhm_init],
            bounds=(0, +np.inf),
        )
        amp_fit = popt[0]
        fwhm_fit_mas = popt[1]
    except ValueError as e:
        # fit failed: a cause: array must not contain infs or NaNs
        print("Interferometric fit failed.")
        print(e)
        amp_fit = np.nan
        fwhm_fit_mas = np.nan
    except RuntimeError as e:
        print("Interferometric fit failed.")
        print(e)
        amp_fit = np.nan
        fwhm_fit_mas = np.nan
    # if mode == 'vis':
    #     popt,pcov = curve_fit(lambda pbl, fwhm: gaussian_fn(pbl, 1.0, fwhm, wavelength),
    #         pbl,y,p0=[fwhm_init],bounds=(0,+np.inf))
    #     amp_fit = 1.0
    #     fwhm_fit_mas = popt[0]
    fwhm_fit_au = fwhm_fit_mas / 1000.0 * dist

    return [fwhm_fit_mas, fwhm_fit_au, amp_fit]


# init_param: amp, fwhm_mas, axis_ratio, PA(rad)
def fit_2D_gaussian_model(
    pbl,
    pbla,
    y,
    wavelength,
    dist=np.nan,
    mode="vis",
    init_param=[np.nan, np.nan, np.nan, np.nan],
):
    n = len(pbl)
    mean = sum(pbl * y) / n
    if math.isnan(init_param[0]):
        init_param[0] = np.nanmax(y)
    if math.isnan(init_param[1]):
        init_param[1] = 2.0  # [mas] #sum(y*(pbl-mean)**2)/n
    if math.isnan(init_param[2]):
        init_param[2] = 0.7
    if math.isnan(init_param[3]):
        init_param[3] = 0.0
    xdata_tuple = (pbl, pbla)
    # if mode == 'corr':
    try:
        popt, pcov = curve_fit(
            lambda xdata_tuple, amp, fwhm_major, axis_ratio, theta: twoD_gaussian_fn_ravel(
                (pbl, pbla), amp, fwhm_major, axis_ratio, theta, wavelength
            ),
            (pbl, pbla),
            y,
            p0=init_param,
            bounds=((0.0, 0.0, 0.0, 0.0), (+np.inf, +np.inf, 1.0, math.pi)),
        )
    except:
        print("Interferometric fit (2D Gaussian) failed.")
    amp_fit = popt[0]
    fwhm_fit_mas = popt[1]
    axis_ratio_fit = popt[2]
    theta_fit = popt[3]  # radian
    # if mode == 'vis':
    # popt,pcov = curve_fit(
    #     lambda xdata_tuple, fwhm_major,axis_ratio,theta: twoD_gaussian_fn_ravel((pbl,pbla),1.0,fwhm_major,axis_ratio,theta,wavelength),
    #     (pbl,pbla),y,p0=init_param,bounds=((0.0,0.0,0.0),(+np.inf,1.0,math.pi)))
    #     amp_fit = 1.0
    #     fwhm_fit_mas = popt[0]
    # axis_ratio_fit = popt[1]
    # theta_fit = popt[2] #radian

    fwhm_fit_au = fwhm_fit_mas / 1000.0 * dist

    return [fwhm_fit_mas, fwhm_fit_au, amp_fit, axis_ratio_fit, theta_fit]


# from menEWS: menEWS_plot.py
def make_uv_tracks(
    uv,
    inp,
    flag,
    ax,
    bases=[],
    symbol="x",
    color="",
    print_station_names=True,
    sel_wl=1.0,
    plot_Mlambda=False,
):
    #
    # from coordinate + ha (range), calculate uv tracks
    #

    ra, dec, BE, BN, BL, base = inp
    lat = -24.62587 * np.pi / 180.0  # Paranal lattitude in radians
    mlim = 2.0  # airmass limit for tracks

    # uv-point
    u, v = uv
    if plot_Mlambda == True:
        u = u / sel_wl
        v = v / sel_wl

    # color depending on flag
    if not color:
        if np.all(flag) == "True":  # VJ 2018
            color = "r"
        else:
            color = "g"
    # color = (not flag)*'g' + flag*'r' #TypeError: an integer is required

    # if not yet done, plot curves
    if base not in bases:
        hamax = np.arccos(
            abs((1.0 / mlim - np.sin(lat) * np.sin(dec)) / (np.cos(lat) * np.cos(dec)))
        )
        harng = np.linspace(-hamax, hamax, 1000)

        ul, vl = calculate_uv_points(inp, harng)
        if plot_Mlambda == True:
            ul /= sel_wl
            vl /= sel_wl
        ax.plot(ul, vl, "-", color="grey", alpha=0.5)
        ax.plot(-ul, -vl, "-", color="grey", alpha=0.5)
        ax.plot([0.0], [0.0], "+k", markersize=5, markeredgewidth=2, alpha=0.5)
        if print_station_names:
            ax.text(-u - 7, -v - 3, base, color="0", alpha=0.8)
        bases.append(base)

    ax.plot(u, v, symbol, color=color, markersize=10, markeredgewidth=3)
    ax.plot(-u, -v, symbol, color=color, markersize=10, markeredgewidth=3)

    return bases


# from menEWS: menEWS_plot.py
def make_uv_plot(
    dic,
    ax,
    verbose=False,
    annotate=False,
    B_lim=(np.nan, np.nan),
    figsize=(5, 5),
    color="",
    print_station_names=True,
    sel_wl=1.0,
    plot_Mlambda=False,
):
    if plot_Mlambda == False:
        sel_wl = 1.0
    try:
        u = dic["VIS2"]["U"]
        v = dic["VIS2"]["V"]
        flag = dic["VIS2"]["FLAG"]
        sta_index = dic["VIS2"]["STA_INDEX"]
        mjd = dic["VIS2"]["MJD"]
    except KeyError as e:
        if verbose:
            print(e)
        u = [0.0]
        v = [0.0]
        flag = [False]
        sta_index = []
        mjd = [0.0]

    uvs = []
    inps = []
    flags = []
    umax = []
    vmax = []
    for j in range(len(u)):
        uvs.append([u[j], v[j]])
        try:
            BE, BN, BL = (
                dic["STAXYZ"][sta_index[j, 0] == dic["STA_INDEX"]][0]
                - dic["STAXYZ"][sta_index[j, 1] == dic["STA_INDEX"]][0]
            )
            sta_label = (
                dic["STA_NAME"][sta_index[j, 0] == dic["STA_INDEX"]][0]
                + "-"
                + dic["STA_NAME"][sta_index[j, 1] == dic["STA_INDEX"]][0]
            )
        except IndexError as e:
            print("make_uv_plot STA_INDEX error.")
            print(e)
            BE, BN, BL = [np.nan, np.nan, np.nan]
            sta_label = ""
        except TypeError as e:
            print("make_uv_plot STA_INDEX error.")
            print(e)
            BE, BN, BL = [np.nan, np.nan, np.nan]
            sta_label = ""
        inps.append(
            [
                dic["RA"] * np.pi / 180.0,
                dic["DEC"] * np.pi / 180.0,
                BE,
                BN,
                BL,
                sta_label,
            ]
        )
        flags.append(flag[j])
    bases = []
    umax = np.nanmax(np.abs(u))
    vmax = np.nanmax(np.abs(v))
    if not (dic["MJD-OBS"]):
        dic["MJD-OBS"] = np.amin(mjd[0])
    try:
        rel_time = (mjd - dic["MJD-OBS"]) * 24.0 * 3600.0  # (s)
        dic["TREL"] = rel_time[0]

        for k, uv in enumerate(uvs):
            bases = make_uv_tracks(
                uv,
                inps[k],
                flags[k],
                ax,
                bases,
                color=color,
                print_station_names=print_station_names,
                sel_wl=sel_wl,
                plot_Mlambda=plot_Mlambda,
            )

        if plot_Mlambda == False:
            xlabel = "$u$ (m)"
            ylabel = "$v$ (m)"
        else:
            xlabel = "$u$ ($M\lambda$)"
            ylabel = "$v$ ($M\lambda$)"
        ax.set_xlim((130, -130))
        ax.set_ylim((-130, 130))
        plotmax = 1.3 * np.amax([umax, vmax])

        plot_title = (
            dic["TARGET"]
            + "\n"
            + "date: "
            + dic["DATE-OBS"]
            + "\n"
            + "TPL start: "
            + dic["TPL_START"]
            + "\n"
            + dic["CATEGORY"]
            + " "
            + dic["BAND"]
            + " "
            + dic["DISPNAME"]
        )  # + ' ' + dic['BCD1'] + '-' + dic['BCD2']
        if math.isnan(B_lim[0]):
            xlim = (+plotmax / sel_wl, -plotmax / sel_wl)
            ylim = (-plotmax / sel_wl, +plotmax / sel_wl)
        else:
            xlim = (+B_lim[1] / sel_wl, -B_lim[1] / sel_wl)
            ylim = (-B_lim[1] / sel_wl, +B_lim[1] / sel_wl)
        # if plot_Mlambda == True:
        plot_config(
            xlabel,
            ylabel,
            plot_title,
            ax,
            dic,
            ylim=ylim,
            xlim=xlim,
            plot_legend=False,
            annotate=annotate,
        )
    except TypeError as e:
        if verbose:
            print("Unable to plot " + "uv")
        if verbose:
            print(e)
        return 1

    return 0


# oi_type = 'vis' (plots OI_VIS2 sqrt(VIS2)), 'vis2' (plots OI_VIS2 sqrt(VIS2)), 'visamp' (plots OI_VIS VISAMP)
def make_plot_with_baseline(
    dic,
    ax,
    tag,
    oi_type="vis",
    ylabel="",
    verbose=False,
    annotate=True,
    annotate_baselines=False,
    B_lim=(np.nan, np.nan),
    sel_wl=np.nan,
    bandwidth=np.nan,
    figsize=(5, 5),
    fit_model=False,
    input_fitsfile="",
    legend_loc="lower left",
    markercolor=None,
    leg_label=None,
):
    xlabel = "$B_\mathrm{p}$ (m)"
    try:
        try:
            x = dic["WLEN"] * 1e6  # (um)
            if oi_type == "vis":
                y = np.sqrt(dic["VIS2"]["VIS2"])
                yerr = np.abs(
                    0.5 * dic["VIS2"]["VIS2ERR"] / np.sqrt(dic["VIS2"]["VIS2"])
                )
                if ylabel == "":
                    ylabel = "Visibility"
            if oi_type == "vis2":
                y = dic["VIS2"]["VIS2"]
                yerr = dic["VIS2"]["VIS2ERR"]
                if ylabel == "":
                    ylabel = "Squared visibility"
            if oi_type == "vis" or oi_type == "vis2":
                u = dic["VIS2"]["U"]
                v = dic["VIS2"]["V"]
                mjd = dic["VIS2"]["MJD"]
                sta_index = dic["VIS2"]["STA_INDEX"]
            if oi_type == "visamp":
                y = dic["VIS"]["VISAMP"]
                yerr = dic["VIS"]["VISAMPERR"]
                u = dic["VIS"]["U"]
                v = dic["VIS"]["V"]
                mjd = dic["VIS"]["MJD"]
                sta_index = dic["VIS"]["STA_INDEX"]
                if ylabel == "":
                    ylabel = "VISAMP"
        except KeyError as e:
            if verbose:
                print("Unable to plot " + tag)
            if verbose:
                print(e)
            return 2

        if not (dic["MJD-OBS"]):
            dic["MJD-OBS"] = np.amin(mjd)

        if sel_wl:  # if it is not None
            if math.isnan(sel_wl):
                if "L" in dic["DETNAME"]:
                    sel_wl = dic["WL_CENTRAL"]
                    bandwidth = 0.2
                if "N" in dic["DETNAME"]:
                    sel_wl = 10.7
                    bandwidth = 0.5
        rel_time = (mjd - dic["MJD-OBS"]) * 24.0 * 3600.0  # (s)
        dic["TREL"] = rel_time[0]
        N = len(mjd)
        try:
            pbl = np.sqrt(u**2 + v**2)
            pbla = np.arctan2(u, v) * 180.0 / np.pi
            idx = np.where(pbla < 0.0)
            pbla[idx] = pbla[idx] + 180.0
        except TypeError as e:
            if verbose:
                print(e)
                print("pbl TypeError")
            pbl = [np.nan] * N
            pbla = [np.nan] * N
        ymin = [0.0]
        ymax = []  # [0.6]

        idx = np.logical_and(x > sel_wl - bandwidth / 2.0, x < sel_wl + bandwidth / 2.0)
        y_new = []
        yerr_new = []
        sta_labels = []
        for j in range(len(y)):
            y_new.append(np.nanmedian(y[j][idx]))
            yerr_new.append(np.nanmean(yerr[j][idx]))
            sta_labels.append(
                dic["STA_NAME"][sta_index[j, 0] == dic["STA_INDEX"]][0]
                + "-"
                + dic["STA_NAME"][sta_index[j, 1] == dic["STA_INDEX"]][0]
            )
        if not np.all(np.isnan(y_new)):
            plt.set_cmap("hsv")
            # print(pbl)
            # print(y_new)
            # print(yerr_new)
            # print(len(pbl),len(y_new),len(yerr_new))
            if markercolor is None:
                mcolor = pbla
            else:
                mcolor = markercolor
            if leg_label is None:
                llabel = "Data"
            else:
                llabel = leg_label
            ax.errorbar(
                pbl,
                y_new,
                yerr=yerr_new,
                fmt=".",
                ecolor="gray",
                alpha=0.6,
                zorder=0,
                capsize=3,
                label="_nolegend_",
            )
            sc_plot = ax.scatter(
                pbl,
                y_new,
                c=mcolor,
                vmin=0.0,
                vmax=180.0,
                marker="o",
                edgecolor="black",
                zorder=5,
                label=llabel,
            )  #'_nolegend_')
            if fit_model == True:
                if tag == "visibility_vs_pbl":
                    model_result = fit_gaussian_model(pbl, y_new, sel_wl, mode="vis")
                    # model_result_2D = fit_2D_gaussian_model(pbl,pbla*np.pi/180.0,y_new,sel_wl,mode='vis') #,init_param=[np.nan,1.9,0.7,0.7])
                else:
                    if tag == "corrflux_vs_pbl":
                        model_result = fit_gaussian_model(
                            pbl, y_new, sel_wl, mode="corr"
                        )
                        # model_result_2D = fit_2D_gaussian_model(pbl,pbla*np.pi/180.0,y_new,sel_wl,mode='corr') #,init_param=[np.nan,1.9,0.7,0.7])
                    else:
                        model_result = fit_gaussian_model(
                            pbl, y_new, sel_wl, mode="corr"
                        )
            else:
                model_result = fit_gaussian_model(pbl, y_new, sel_wl, mode="corr")
            if annotate == True or (annotate == False and annotate_baselines == True):
                for j in range(len(y)):
                    ax.text(
                        pbl[j] - 1.0,
                        y_new[j],
                        sta_labels[j],
                        color="gray",
                        alpha=0.6,
                        ha="right",
                    )
            pmin = np.nanmin(y_new)  # np.nanpercentile(y_new, 10.0)
            pmax = np.nanmax(y_new)  # np.nanpercentile(y_new, 90.0)
            pmid = (pmin + pmax) / 2.0
            pdiff = pmax - pmin
            lim_fact = 0.2
            ymin.append(pmin - pdiff * lim_fact)
            ymax.append(pmax + pdiff * lim_fact)
            # if ('AVGCAL_INT' in input_fitsfile) or ('AVGFLUXCAL_INT' in input_fitsfile) or \
            #            ('BCDCAL_INT' in input_fitsfile) or ('FINALCAL_INT' in input_fitsfile):
            #    plot_title = dic['TARGET'] + "\n" + "date: " + dic['DATE-OBS'] + "\n" + dic['CATEGORY'] + ' ' + \
            #        dic['BAND'] + ' ' + dic['DISPNAME'] + ' ' + ', ' + '%.2f $\mu$m'%sel_wl
            # else:
            plot_title = (
                dic["TARGET"]
                + "\n"
                + "date: "
                + dic["DATE-OBS"]
                + "\n"
                + "TPL start: "
                + dic["TPL_START"]
                + "\n"
                + dic["CATEGORY"]
                + " "
                + dic["BAND"]
                + " "
                + dic["DISPNAME"]
                + " "
                + dic["BCD1"]
                + "-"
                + dic["BCD2"]
                + ", "
                + "%.2f $\mu$m" % sel_wl
            )

            left, right = ax.set_xlim()
            if markercolor is None:
                cbar = plt.colorbar(sc_plot, orientation="vertical", ax=ax)
                cbar.ax.set_ylabel("Baseline PA ($^\circ$)")

            if math.isnan(B_lim[0]):
                xlim = (0.0, right)
            else:
                xlim = B_lim
            plot_leg = True  # False
            if fit_model == True:
                # if tag == 'visibility_vs_pbl' or tag == 'corrflux_vs_pbl':
                pbl_model = np.arange(0.0, xlim[1], 1.0)
                y_model = gaussian_fn(
                    pbl_model, model_result[2], model_result[0], sel_wl
                )
                ax.plot(
                    pbl_model,
                    y_model,
                    "-r",
                    label="1D Gaussian fit\n FWHM = %.2f mas" % model_result[0],
                )
                plot_leg = True
                # y_model2 = twoD_gaussian_fn((pbl,pbla*np.pi/180.0),model_result_2D[2],model_result_2D[0],model_result_2D[3],model_result_2D[4],sel_wl)
                # ax.scatter(pbl, y_model2,c=pbla,s=60,vmin=0.0,vmax=180.0,marker='X', edgecolor='black',zorder=10,
                #     label='2D Gaussian fit \n $\mathrm{FWHM}_a = %.2f$ mas \n $b/a = %.2f$ \n PA = %.1f$^\circ$'
                #     % (model_result_2D[0],model_result_2D[3],model_result_2D[4]*180.0/math.pi))
            plot_config(
                xlabel,
                ylabel,
                plot_title,
                ax,
                dic,
                ylim=(np.nanmin(ymin), np.nanmax(ymax)),
                xlim=xlim,
                plot_legend=plot_leg,
                annotate=annotate,
                legend_loc=legend_loc,
            )
    except TypeError as e:
        if verbose:
            print("Unable to plot " + tag)
        if verbose:
            print(e)
        return 1

    return 0


# tag =  'flux', 'visibility', 'corrflux','allflux' 'visamp', 'visphi', 't3phi', 't3amp'
# oi_type = 'vis','vis2','visamp','visphi','t3amp','t3phi','flux','tf','tf2'
# legend_style='full': prints station names, baseline lengths and PAs; or 'short': prints only station names, 'no': '_nolegend_'
# plot_errorbars: 'full', 'some', 'no', 'fill'
# ax: if a single value: all plots will be plotted on one axis, if a list of axes: lines will be plotted one by one
# xaxis: 'wl' (default) or 'sp_freq'
# color_by_wl: only applied if xaxis is 'sp_freq'
def make_plot_with_wavelength(
    dic,
    ax,
    tag,
    oi_type="vis",
    ylabel="",
    verbose=False,
    annotate=True,
    wl_lim=(np.nan, np.nan),
    figsize=(5, 5),
    input_fitsfile="",
    legend_loc="upper right",
    plot_legend=True,
    legend_ncol=1,
    legend_style="full",
    n_max_legend_entries=6,
    convolution_kernel_width=0.0,
    normalize=0,
    linestyle="-",
    linecolors=None,
    custom_legend_label="",
    alpha=1.0,
    plot_errorbars="no",
    err_show_ratio=0.04,
    pbls=None,
    ylim_percentiles=[1.0, 97.0],
    nodata=False,
    xaxis="wl",
    plot_cbar=False,
    color_by_wl=False,
    plot_bl_indices=None,
):
    if xaxis == "wl":
        xlabel = "Wavelength ($\mu$m)"
    if xaxis == "sp_freq":
        xlabel = "$B_\mathrm{p}$ (M$\lambda$)"
        if oi_type == "t3amp" or oi_type == "t3phi":
            xlabel = "$B_\mathrm{p,max}$ (M$\lambda$)"
    try:
        ymin = []
        ymax = []
        # print('1')
        try:
            wl = dic["WLEN"] * 1e6
            x = wl
            if oi_type == "vis":
                y = np.sqrt(dic["VIS2"]["VIS2"])
                yerr = np.abs(
                    0.5 * dic["VIS2"]["VIS2ERR"] / np.sqrt(dic["VIS2"]["VIS2"])
                )
                if ylabel == "":
                    ylabel = "Visibility"
            if oi_type == "vis2":
                y = dic["VIS2"]["VIS2"]
                yerr = dic["VIS2"]["VIS2ERR"]
                if ylabel == "":
                    ylabel = "Squared visibility"
            if oi_type == "vis" or oi_type == "vis2":
                u = dic["VIS2"]["U"]
                v = dic["VIS2"]["V"]
                mjd = dic["VIS2"]["MJD"]
                sta_index = dic["VIS2"]["STA_INDEX"]
            if oi_type == "visamp":
                y = dic["VIS"]["VISAMP"]
                yerr = dic["VIS"]["VISAMPERR"]
                if ylabel == "":
                    if tag == "corrflux" or tag == "allflux":
                        ylabel = "Correlated flux (Jy)"
                    else:
                        ylabel = "VISAMP"
                if tag == "allflux":
                    y2 = dic["FLUX"]["FLUX"]
                    y2err = dic["FLUX"]["FLUXERR"]
                    y3 = np.concatenate((y, y2))
                    y = y3
                    y3err = np.concatenate((yerr, y2err))
                    yerr = y3err
            if oi_type == "visphi":
                y = dic["VIS"]["VISPHI"]
                yerr = dic["VIS"]["VISPHIERR"]
                if ylabel == "":
                    ylabel = "Differential phase ($^\circ$)"
            if oi_type == "visamp" or oi_type == "visphi":
                u = dic["VIS"]["U"]
                v = dic["VIS"]["V"]
                mjd = dic["VIS"]["MJD"]
                sta_index = dic["VIS"]["STA_INDEX"]
                if tag == "allflux":
                    y2 = dic["FLUX"]["FLUX"]
                    u2 = np.zeros_like(y2[:, 0])
                    v2 = np.zeros_like(y2[:, 0])
                    mjd2 = dic["FLUX"]["MJD"]
                    sta_index2 = dic["FLUX"]["STA_INDEX"]
                    u = np.concatenate((u, u2))
                    v = np.concatenate((v, v2))
                    mjd = np.concatenate((mjd, mjd2))
                    sta_index = np.concatenate(
                        (sta_index, np.transpose(np.tile(sta_index2, (2, 1))))
                    )
            if oi_type == "t3amp":
                y = dic["T3"]["T3AMP"]
                yerr = dic["T3"]["T3AMPERR"]
                if ylabel == "":
                    ylabel = "T3AMP"
            if oi_type == "t3phi":
                y = dic["T3"]["T3PHI"]
                yerr = dic["T3"]["T3PHIERR"]
                if ylabel == "":
                    ylabel = "Closure phase ($^\circ$)"
            if oi_type == "t3amp" or oi_type == "t3phi":
                u = [dic["T3"]["U1"], dic["T3"]["U2"]]
                v = [dic["T3"]["V1"], dic["T3"]["V2"]]
                mjd = dic["T3"]["MJD"]
                sta_index = dic["T3"]["STA_INDEX"]
            if oi_type == "flux":
                y = dic["FLUX"]["FLUX"]
                yerr = dic["FLUX"]["FLUXERR"]
                u = []
                v = []
                mjd = dic["FLUX"]["MJD"]
                sta_index = dic["FLUX"]["STA_INDEX"]
                if ylabel == "":
                    ylabel = "Flux (Jy)"
            if oi_type == "tf":
                y = np.sqrt(dic["TF2"]["TF2"])
                if ylabel == "":
                    ylabel = "TF"
            if oi_type == "tf2":
                y = dic["TF2"]["TF2"]
                if ylabel == "":
                    ylabel = "TF2"
            if oi_type == "tf" or oi_type == "tf2":
                u = []
                v = []
                mjd = dic["TF2"]["MJD"]
                sta_index = dic["TF2"]["STA_INDEX"]
        except KeyError as e:
            if verbose:
                print("Unable to plot " + tag)
            if verbose:
                print(e)
            return 2

        if not (dic["MJD-OBS"]):
            dic["MJD-OBS"] = np.amin(mjd)

        if "mjd" in locals():
            rel_time = (mjd - dic["MJD-OBS"]) * 24.0 * 3600.0  # (s)
            dic["TREL"] = rel_time[0]
            N = len(mjd)
        sort_idx = np.array([])
        if tag != "flux" and tag != "t3amp" and tag != "t3phi":
            try:
                pbl = np.sqrt(u**2 + v**2)
                pbla = np.arctan2(u, v) * 180.0 / np.pi
                idx = np.where(pbla < 0.0)
                pbla[idx] = pbla[idx] + 180.0
                sort_idx = np.argsort(pbl)
            except TypeError as e:
                if verbose:
                    print(e)
                    print("pbl TypeError")
                pbl = [np.nan] * N
                pbla = [np.nan] * N
        elif oi_type == "t3phi" or oi_type == "t3amp":
            pbl1 = np.sqrt(u[0] ** 2 + v[0] ** 2)
            pbl2 = np.sqrt(u[1] ** 2 + v[1] ** 2)
            pbl3 = np.sqrt((u[0] + u[1]) ** 2 + (v[0] + v[1]) ** 2)
            pbl = np.max([pbl1, pbl2, pbl3], axis=0)
            pbla = [np.nan] * N
        # if oi_type == 'vis2' or oi_type == 'vis' or oi_type == 'visamp' \
        #     or oi_type == 'tf2' or oi_type == 't3amp':
        #     ymin = [0.0]
        # ymax = [1.0]
        # ymin = [0.0]

        if sort_idx.size:
            for_range = sort_idx
            if pbls is not None:
                pbls.append(pbl[sort_idx])
        else:
            for_range = range(len(y))
            if pbls is not None:
                if oi_type == "t3phi":
                    pbl1 = np.sqrt(u[0] ** 2 + v[0] ** 2)
                    pbl2 = np.sqrt(u[1] ** 2 + v[1] ** 2)
                    pbl3 = np.sqrt((u[0] + u[1]) ** 2 + (v[0] + v[1]) ** 2)
                    pbls.append([pbl1, pbl2, pbl3])
                else:
                    pbls.append(pbl)

        if linecolors is None:
            if isinstance(ax, np.ndarray):
                colors_ = ["blue"]
            else:
                colors_ = [
                    "#1f77b4",
                    "#ff7f0e",
                    "#2ca02c",
                    "#d62728",
                    "#9467bd",
                    "#8c564b",
                    "#e377c2",
                    "#7f7f7f",
                    "#bcbd22",
                    "#17becf",
                ]
        else:
            colors_ = linecolors
        for j, color in zip(for_range, itertools.cycle(colors_)):
            if plot_bl_indices is None or j in plot_bl_indices:
                if oi_type == "flux":
                    try:
                        sta_label = (
                            dic["STA_NAME"][sta_index[j] == dic["STA_INDEX"]][0]
                            + ", "
                            + dic["TEL_NAME"][sta_index[j] == dic["STA_INDEX"]][0]
                        )
                        lab = r"%s" % (sta_label)
                    except IndexError as e:
                        print(e)
                        lab = r""
                try:
                    if oi_type == "t3amp" or oi_type == "t3phi":
                        sta_label = (
                            dic["STA_NAME"][sta_index[j, 0] == dic["STA_INDEX"]][0]
                            + "-"
                            + dic["STA_NAME"][sta_index[j, 1] == dic["STA_INDEX"]][0]
                            + "-"
                            + dic["STA_NAME"][sta_index[j, 2] == dic["STA_INDEX"]][0]
                        )
                        lab = r"%s" % (sta_label)
                    if oi_type != "flux" and oi_type != "t3amp" and oi_type != "t3phi":
                        sta_label = (
                            dic["STA_NAME"][sta_index[j, 0] == dic["STA_INDEX"]][0]
                            + "-"
                            + dic["STA_NAME"][sta_index[j, 1] == dic["STA_INDEX"]][0]
                        )
                        if legend_style == "full":
                            lab = r"%s $B_\mathrm{p}=%.1f$ m, $\phi = %.1f^{\circ}$" % (
                                sta_label,
                                pbl[j],
                                pbla[j],
                            )
                        if legend_style == "short":
                            # lab = r'%s $%.1f$ m' % (sta_label, pbl[j])
                            lab = r"%s $B_\mathrm{p}=%.1f$ m" % (sta_label, pbl[j])
                        if legend_style == "no":
                            lab = "_nolegend_"
                except IndexError as e:
                    print("Cannot make legend entry.")
                    print(e)
                    lab = r""
                if tag == "allflux":
                    if pbl[j] < 0.0001:
                        lab = "Single-dish spectrum"
                if j == n_max_legend_entries:
                    lab = "..."
                if j > n_max_legend_entries:
                    lab = "_nolegend_"
                if custom_legend_label != "":
                    lab = custom_legend_label + " " + lab
                Ny = len(y[j])

                if xaxis == "sp_freq":
                    x = pbl[j] / wl
                if len(y[j]) > 1:
                    if nodata == False:
                        if convolution_kernel_width > 0.0:
                            # Create kernel
                            g = Gaussian1DKernel(stddev=convolution_kernel_width)
                            # Convolve data
                            y[j] = convolve(y[j], g)
                            yerr[j] = convolve(yerr[j], g)
                        if normalize == 1:
                            y[j] = y[j] / np.max(y[j])
                        if isinstance(ax, np.ndarray):
                            (base_line,) = ax.flatten()[j].plot(
                                x, y[j], label=lab, linestyle=linestyle, alpha=alpha
                            )
                        else:
                            if xaxis == "sp_freq":
                                if not (math.isnan(wl_lim[0])) and not (
                                    math.isnan(wl_lim[1])
                                ):
                                    idx = np.logical_and(
                                        wl > (wl_lim[0]), wl <= (wl_lim[1])
                                    )
                                    xp = x[idx]
                                    yp = y[j][idx]
                                    yerrp = yerr[j][idx]
                                    wlp = wl[idx]
                                    wl_min = wl_lim[0]
                                    wl_max = wl_lim[1]
                                else:
                                    xp = x
                                    yp = y[j]
                                    yerrp = yerr[j]
                                    wlp = wl
                                    wl_min = np.min(wl)
                                    wl_max = np.max(wl)
                                if (
                                    plot_errorbars == "full"
                                    or plot_errorbars == "some"
                                    or plot_errorbars == "fill"
                                ):

                                    # XX = np.concatenate([xp, np.flip(xp)])
                                    # YY = np.concatenate([yp-yerrp, np.flip(yp+yerrp)])
                                    # #print(XX.shape,YY.shape,yp.shape,yerrp.shape)
                                    # ax.fill(XX, YY, alpha=0.33,color='grey')
                                    err_idx = np.where(wlp > 0.0)
                                    XX = np.concatenate(
                                        [xp[err_idx], np.flip(xp[err_idx])]
                                    )
                                    # print(XX)
                                    YY = np.nan_to_num(
                                        np.concatenate(
                                            [
                                                yp[err_idx] - yerrp[err_idx],
                                                np.flip(yp[err_idx] + yerrp[err_idx]),
                                            ]
                                        )
                                    )
                                    # print(YY)
                                    if isinstance(ax, np.ndarray):
                                        ax.flatten()[j].fill(
                                            XX, YY, alpha=0.1, color="grey"
                                        )
                                    else:
                                        ax.fill(XX, YY, alpha=0.1, color="grey")

                                if len(yp) < 10:
                                    ax.plot(xp, yp, alpha=0.5, color="black")
                                    plt.set_cmap("rainbow")
                                    sc_plot = ax.scatter(
                                        xp,
                                        yp,
                                        s=60,
                                        c=wlp,
                                        vmin=wl_min,
                                        vmax=wl_max,
                                        marker="o",
                                        edgecolor="black",
                                        zorder=5,
                                    )

                                else:
                                    if color_by_wl:
                                        points = np.array([xp, yp]).T.reshape(-1, 1, 2)
                                        segments = np.concatenate(
                                            [points[:-1], points[1:]], axis=1
                                        )
                                        # Create a continuous norm to map from data points to colors
                                        norm = plt.Normalize(wl_min, wl_max)
                                        from matplotlib.collections import (
                                            LineCollection,
                                        )
                                        from matplotlib.colors import (
                                            BoundaryNorm,
                                            ListedColormap,
                                        )

                                        lc = LineCollection(
                                            segments, cmap="rainbow", norm=norm
                                        )
                                        # Set the values used for colormapping
                                        lc.set_array(wlp)
                                        lc.set_linewidth(2)
                                        line = ax.add_collection(lc)
                                    else:
                                        (base_line,) = ax.plot(
                                            xp,
                                            yp,
                                            color=color,
                                            label=lab,
                                            linestyle=linestyle,
                                            alpha=alpha,
                                        )
                                # plt.colorbar(line, ax=ax)
                            else:
                                (base_line,) = ax.plot(
                                    x,
                                    y[j],
                                    color=color,
                                    label=lab,
                                    linestyle=linestyle,
                                    alpha=alpha,
                                )
                                # if tag == 'allflux':
                                #    ax.plot(x, y2[0], color='black',label='Single-dish spectrum',linestyle=linestyle,alpha=alpha)

                else:
                    if nodata == False:
                        if isinstance(ax, np.ndarray):
                            (base_line,) = ax.flatten()[j].plot(
                                x,
                                y[j],
                                "o",
                                label=lab,
                                linestyle=linestyle,
                                alpha=alpha,
                            )
                            ax.flatten()[j].errorbar(
                                x,
                                y[j],
                                yerr=yerr[j],
                                fmt="none",
                                ecolor=base_line.get_color(),
                                capsize=3,
                                label="_nolegend_",
                            )
                        else:
                            (base_line,) = ax.plot(
                                x,
                                y[j],
                                "o",
                                label=lab,
                                linestyle=linestyle,
                                alpha=alpha,
                            )
                            ax.errorbar(
                                x,
                                y[j],
                                yerr=yerr[j],
                                fmt="none",
                                ecolor=color,
                                capsize=3,
                                label="_nolegend_",
                            )
                # if j==3: plt.show()
                if (
                    plot_errorbars == "full"
                    or plot_errorbars == "some"
                    or plot_errorbars == "fill"
                ):
                    if xaxis == "wl":
                        if plot_errorbars == "fill":
                            if nodata == False:
                                err_idx = np.where(wl > 0.0)
                                XX = np.concatenate([x[err_idx], np.flip(x[err_idx])])
                                # print(XX)
                                YY = np.nan_to_num(
                                    np.concatenate(
                                        [
                                            y[j][err_idx] - yerr[j][err_idx],
                                            np.flip(y[j][err_idx] + yerr[j][err_idx]),
                                        ]
                                    )
                                )
                                # print(YY)
                                if isinstance(ax, np.ndarray):
                                    ax.flatten()[j].fill(
                                        XX, YY, alpha=0.1, color=base_line.get_color()
                                    )
                                else:
                                    ax.fill(
                                        XX, YY, alpha=0.1, color=base_line.get_color()
                                    )
                        else:
                            if plot_errorbars == "some":
                                err_idx = np.sort(
                                    np.random.randint(0, Ny, int(err_show_ratio * Ny))
                                )
                                cs = 3
                                al = 0.9
                            elif plot_errorbars == "full":
                                err_idx = np.where(wl > 0.0)
                                cs = 0
                                al = 0.25
                            if nodata == False:
                                if isinstance(ax, np.ndarray):
                                    ax.flatten()[j].errorbar(
                                        x[err_idx],
                                        y[j][err_idx],
                                        yerr=yerr[j][err_idx],
                                        fmt="none",
                                        ecolor=base_line.get_color(),
                                        capsize=cs,
                                        alpha=al,
                                        label="_nolegend_",
                                    )
                                else:
                                    ax.errorbar(
                                        x[err_idx],
                                        y[j][err_idx],
                                        yerr=yerr[j][err_idx],
                                        fmt="none",
                                        ecolor=color,
                                        capsize=cs,
                                        alpha=al,
                                        label="_nolegend_",
                                    )  # ecolor=base_line.get_color()
                        # print(j)
                # wavelengths  considered calculation of y plot limits
                # M_idx = np.logical_not(np.logical_or(np.logical_and(x > 1.6,x < 1.8),np.logical_and(x > 3.0,x < 4.0),np.logical_and(x > 8.0,x < 12.5)))
                if "MATISSE" in dic["INSTRUMENT"]:
                    if math.isnan(wl_lim[0]) or math.isnan(wl_lim[1]):
                        M_idx = np.logical_not(
                            np.logical_and(wl > 1.6, wl < 1.8)
                            + np.logical_and(wl > 3.0, wl < 4.0)
                            + np.logical_and(wl > 8.0, wl < 12.5)
                        )
                    else:
                        M_idx = np.logical_not(
                            np.logical_and(wl > 1.6, wl < 1.8)
                            + np.logical_and(wl > 3.0, wl < 4.0)
                            + np.logical_and(wl > 8.0, wl < 12.5)
                            + np.logical_and(wl > wl_lim[0], wl < wl_lim[1])
                        )

                y_new = y[j]
                # print(M_idx)
                # print(x)
                # print(y[j])
                if "MATISSE" in dic["INSTRUMENT"]:
                    y_new[M_idx] = np.nan
                # print(y_new)
                if "PIONIER" not in dic["INSTRUMENT"]:
                    if oi_type == "flux":
                        pmin = np.nanpercentile(
                            y_new[int(0.2 * Ny) : int(0.8 * Ny)], ylim_percentiles[0]
                        )
                        pmax = np.nanpercentile(
                            y_new[int(0.2 * Ny) : int(0.8 * Ny)], ylim_percentiles[1]
                        )
                        pmid = (pmin + pmax) / 2.0
                        pdiff = pmax - pmin
                        lim_fact = 1.0
                        ymin.append(pmin - pdiff * lim_fact)
                        ymax.append(pmax + pdiff * lim_fact)
                    else:
                        # print(np.nanmin(y_new),np.nanmax(y_new))
                        pmin = np.nanpercentile(
                            y_new[int(0.2 * Ny) : int(0.8 * Ny)], ylim_percentiles[0]
                        )
                        pmax = np.nanpercentile(
                            y_new[int(0.2 * Ny) : int(0.8 * Ny)], ylim_percentiles[1]
                        )
                        pmid = (pmin + pmax) / 2.0
                        pdiff = pmax - pmin
                        lim_fact = 2.0
                        ymin.append(pmid - pdiff / 2.0 * lim_fact)
                        ymax.append(pmid + pdiff / 2.0 * lim_fact)
                        # print('p',pmin,pmax)
                else:
                    pmin = np.nanmin(y_new)
                    pmax = np.nanmax(y_new)
                    pmid = (pmin + pmax) / 2.0
                    pdiff = pmax - pmin
                    lim_fact = 1.0
                    ymin.append(pmin - pdiff * lim_fact)
                    ymax.append(pmax + pdiff * lim_fact)
        if xaxis == "sp_freq":
            if color_by_wl:
                if plot_cbar == True:
                    if len(yp) < 10:
                        cbar = plt.colorbar(sc_plot, orientation="vertical", ax=ax)
                    else:
                        cbar = plt.colorbar(line, ax=ax)
                    cbar.ax.set_ylabel("Wavelength ($\mu$m)")

        # print(dic['TARGET'],dic['DATE-OBS'],dic['CATEGORY'],dic['BAND'],dic['DISPNAME'],dic['BCD1'],dic['BCD2'])
        # if ('AVGCAL_INT' in input_fitsfile) or ('AVGFLUXCAL_INT' in input_fitsfile) or \
        #                ('BCDCAL_INT' in input_fitsfile) or ('FINALCAL_INT' in input_fitsfile):
        #    plot_title = dic['TARGET'] + "\n" + "date: " + dic['DATE-OBS'] + "\n" + dic['CATEGORY'] + ' ' + \
        #                    dic['BAND'] + ' ' + dic['DISPNAME'] #+ ' ' + dic['BCD1'] + '-' + dic['BCD2']
        # else:
        # print(ymin)
        # print(ymax)
        # plt.show()

        plot_title = (
            dic["TARGET"]
            + "\n"
            + "date: "
            + dic["DATE-OBS"]
            + "\n"
            + "TPL start: "
            + dic["TPL_START"]
            + "\n"
            + dic["CATEGORY"]
            + " "
            + dic["BAND"]
            + " "
            + dic["DISPNAME"]
            + " "
            + dic["BCD1"]
            + "-"
            + dic["BCD2"]
        )
        if xaxis == "wl":
            wl_lim_temp = wl_lim
            if wl_lim_temp[0]:
                if math.isnan(wl_lim_temp[0]):
                    wl_lim_temp = (None, None)
        if xaxis == "sp_freq":
            wl_lim_temp = (np.min(pbl) / wl_lim[1], np.max(pbl) / wl_lim[0])
        # print('nm',np.nanmin(ymin), np.nanmax(ymax))
        # print(ymin,ymax)
        plot_config(
            xlabel,
            ylabel,
            plot_title,
            ax,
            dic,
            ylim=(np.nanmin(ymin), np.nanmax(ymax)),
            xlim=wl_lim_temp,
            annotate=annotate,
            # ylim=(-80.0,30.0),xlim=wl_lim_temp,annotate=annotate,
            legend_loc=legend_loc,
            legend_ncol=legend_ncol,
            plot_legend=plot_legend,
        )
        # print('2')
    except TypeError as e:  # TypeError as e:
        print("Unable to plot " + tag)
        print(e)
        return 1
    # print('3')
    return 0


###############################################################################################################
# show_allred
###############################################################################################################
# the main plotting function: makes plots from the results of the DRS reduction
# arguments:
# 	datadir: the directory containing the fits files produced by the DRS
# 	outputdir (optional): directory where the plots are placed
# 	verbose (optional): turn on verbose messages
def show_allred(
    datadir,
    outputdir="plots",
    fn_pattern="",
    verbose=False,
    save_png=True,
    save_eps=False,
    sel_wls=[np.nan],
    bandwidths=[np.nan],
    file_type="",
    pro_catg="",
    nbProc=1,
    annotate=True,
    wl_lim=(np.nan, np.nan),
    B_lim=(np.nan, np.nan),
    fit_model=False,
    figsize=(5, 5),
):
    # check if output directory exists
    # if not, create it
    fitsfiles = sorted(glob.glob(datadir + "/*" + fn_pattern + "*fits"))
    if fitsfiles != []:
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
    # for input_fitsfile in fitsfiles:
    #     make_plots(input_fitsfile, outputdir=outputdir, fn_pattern = fn_pattern, verbose=verbose, save_png=save_png,
    #     save_eps=save_eps,sel_wls=sel_wls,bandwidths=bandwidths,file_type=file_type,pro_catg=pro_catg)
    pool = mp.Pool(processes=nbProc)
    make_plots_fn = functools.partial(
        make_plots,
        outputdir=outputdir,
        fn_pattern=fn_pattern,
        verbose=verbose,
        save_png=save_png,
        save_eps=save_eps,
        sel_wls=sel_wls,
        bandwidths=bandwidths,
        file_type=file_type,
        pro_catg=pro_catg,
        annotate=annotate,
        wl_lim=wl_lim,
        B_lim=B_lim,
        fit_model=fit_model,
        figsize=figsize,
    )
    pool.map(make_plots_fn, fitsfiles)
    pool.terminate()
    pool.join()


###############################################################################################################
# show_allred_mosaic
###############################################################################################################
# makes mosaic plots from the results of the DRS reduction
# arguments:
# 	datadir: the directory containing the fits files produced by the DRS
# 	outputdir (optional): directory where the plots are placed
# 	verbose (optional): turn on verbose messages
#   plot_errorbars_wlplot: 'full', 'fill', 'some', 'no'
def show_allred_mosaic(
    datadir,
    outputdir="plots",
    fn_pattern="",
    oi_types_list=[
        ["uv", "vis_pbl", "visamp_pbl", "flux", "vis", "visamp", "visphi", "t3phi"]
    ],
    verbose=False,
    save_png=True,
    save_eps=False,
    sel_wl=np.nan,
    bandwidth=np.nan,
    annotate=False,
    wl_lim=(np.nan, np.nan),  # file_type='',pro_catg=''
    B_lim=(np.nan, np.nan),
    fit_model=False,
    figsize=(15, 15),
    ext=0,
    redo_overwrite=True,
    plot_errorbars_wlplot="fill",
    err_show_ratio=0.04,
):
    # check if output directory exists
    # if not, create it
    fitsfiles = sorted(glob.glob(datadir + "/*" + fn_pattern + "*fits"))
    if fitsfiles != []:
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
    for input_fitsfile in fitsfiles:
        if ext > 0:
            ext_tag = "_%d" % ext
        else:
            ext_tag = ""
        output_file_path = (
            outputdir
            + "/"
            + os.path.splitext(os.path.basename(input_fitsfile))[0]
            + ext_tag
            + "_mosaic"
        )
        if redo_overwrite:
            make_plot = True
        else:
            # check if plot file already exists
            if os.path.exists(output_file_path + ".png"):
                make_plot = False
                print(
                    os.path.basename(output_file_path)
                    + " already exists. Skip making plot."
                )
            else:
                make_plot = True
        if make_plot:
            make_mosaic_plot(
                [input_fitsfile],
                output_file_path,
                oi_types_list=oi_types_list,
                verbose=verbose,
                save_png=save_png,
                save_eps=save_eps,
                sel_wl=sel_wl,
                bandwidth=bandwidth,
                annotate=annotate,
                legend_loc="upper right",
                legend_loc_pbl="lower left",
                wl_lim=wl_lim,
                B_lim=B_lim,
                fit_model=fit_model,
                figsize=(15, 15),
                ext=ext,
                plot_errorbars_wlplot=plot_errorbars_wlplot,
                err_show_ratio=err_show_ratio,
            )


###############################################################################################################
# make_plots
###############################################################################################################
# default wavelengths: in L band: the central wavelength from the observation, bandwidth = 0.2 um
#                      in N band: 10.7 um, bandwidth = 0.2 um
def make_plots(
    input_fitsfile,
    outputdir="plots",
    fn_pattern="",
    verbose=False,
    save_png=True,
    save_eps=False,
    sel_wls=[np.nan],
    bandwidths=[np.nan],
    file_type="",
    pro_catg="",
    annotate=True,
    wl_lim=(np.nan, np.nan),
    B_lim=(np.nan, np.nan),
    fit_model=False,
    figsize=(5, 5),
):
    print(input_fitsfile)
    # open fits file
    dic = {}
    if "OBJ_CORR_FLUX" not in input_fitsfile and "CALIB_CAL" not in input_fitsfile:
        dic = open_fits(input_fitsfile, verbose)
    if dic:
        if not (pro_catg == ""):
            if not (dic["PRO_CATG"] == pro_catg):
                return

        #####################################
        # make uv plot:
        #####################################
        if (
            ("RAW_INT" in input_fitsfile)
            or ("CAL_INT" in input_fitsfile)
            or ("oifits" in input_fitsfile)
            or ("OIFITS" in input_fitsfile or file_type == "oifits")
        ):
            if verbose:
                print("Make uv-plot.")
            fig, ((ax1)) = plt.subplots(
                1, 1, sharey=False, sharex=False, figsize=figsize
            )
            return_val = make_uv_plot(
                dic,
                ax1,
                verbose=verbose,
                annotate=annotate,
                B_lim=B_lim,
                figsize=figsize,
            )
            if return_val == 0:
                outputfig = (
                    outputdir
                    + "_".join(os.path.basename(input_fitsfile).split(".")[:-1])
                    + "_"
                    + "uv"
                )
                # outputfig = outputdir + '/uv'
                if save_png:
                    plt.savefig(outputfig + ".png", dpi=200)
                if save_eps:
                    plt.savefig(outputfig + ".eps", format="eps", dpi=300)
            plt.close(fig)

        #####################################
        # plot as a function of baseline
        #####################################
        if (
            ("RAW_INT" in input_fitsfile)
            or ("CAL_INT" in input_fitsfile)
            or ("oifits" in input_fitsfile)
            or ("OIFITS" in input_fitsfile or file_type == "oifits")
        ):
            if verbose:
                print("Plots as a function of baseline.")

            if "CAL_INT" not in input_fitsfile:
                ylabels = ["Raw squared visibility", "Raw visibility"]
            else:
                ylabels = ["Squared visibility", "Visibility"]
            tags = ["squaredvisibility_vs_pbl", "visibility_vs_pbl"]
            oi_types = ["vis2", "vis"]

            if "correlated" in dic["VIS"]["AMPTYP"]:
                if "CAL_INT" not in input_fitsfile:
                    ylabels.append("Raw correlated flux")
                else:
                    ylabels.append("Correlated flux (Jy)")
                tags.append("corrflux_vs_pbl")
                oi_types.append("visamp")

            for k in range(len(tags)):
                for l in range(len(sel_wls)):
                    if sel_wls[l]:  # if it is not None
                        if math.isnan(sel_wls[l]):
                            if "L" in dic["DETNAME"]:
                                sel_wls[l] = dic["WL_CENTRAL"]
                                bandwidths[l] = 0.2
                            if "N" in dic["DETNAME"]:
                                sel_wls[l] = 10.7
                                bandwidths[l] = 0.5
                    if (sel_wls[l] % 1) == 0:
                        lambda_tag = "%dum" % sel_wls[l]
                    else:
                        int_part = np.floor(sel_wls[l])
                        fract_part = sel_wls[l] - 1.0 * int_part
                        lambda_tag = "%dum%d" % (int_part, np.round(fract_part * 10.0))
                    fig, ((ax1)) = plt.subplots(
                        1, 1, sharey=False, sharex=False, figsize=figsize
                    )
                    return_val = make_plot_with_baseline(
                        dic,
                        ax1,
                        tags[k],
                        oi_type=oi_types[k],
                        ylabel=ylabels[k],
                        annotate=True,
                        B_lim=B_lim,
                        verbose=False,
                        sel_wl=sel_wls[l],
                        bandwidth=bandwidths[l],
                        figsize=figsize,
                        fit_model=fit_model,
                    )
                    if return_val == 0:
                        outputfig = (
                            outputdir
                            + "_".join(os.path.basename(input_fitsfile).split(".")[:-1])
                            + "_"
                            + tags[k]
                            + "_"
                            + lambda_tag
                        )
                        if save_png:
                            plt.savefig(outputfig + ".png", dpi=200)
                        if save_eps:
                            plt.savefig(outputfig + ".eps", format="eps", dpi=300)
                    plt.close(fig)

        #####################################
        # plots as a function of wavelength
        #####################################
        # if 'CALIB_CAL' in input_fitsfile:
        # it contains many image frames, not practical to plot
        # if 'OBJ_CORR_FLUX' in input_fitsfile:
        # it contains many image frames, not practical to plot
        # if 'PHOT_BEAMS' in input_fitsfile:
        # it contains many image frames, not practical to plot
        tags = []
        if verbose:
            print("Plots as a function of wavelength.")
        if (
            ("RAW_INT" in input_fitsfile)
            or ("CAL_INT" in input_fitsfile)
            or ("oifits" in input_fitsfile)
            or ("OIFITS" in input_fitsfile or file_type == "oifits")
        ):
            # check if the CFXAMP is present
            # try:
            #     a = dic['VIS']['CFXAMP']
            # except KeyError as e:
            #     if verbose: print(e)
            #     dic['VIS']['CFXAMP'] = 0.0*dic['VIS']['VISAMP']
            oi_types = ["vis2", "vis", "visamp", "visphi", "flux", "t3phi"]
            if "CAL_INT" not in input_fitsfile:
                ylabels = [
                    "Raw squared visibility",
                    "Raw visibility",
                    "VISAMP",
                    "VISPHI ($^\circ$)",
                    "Raw total flux",
                    "Raw closure phase ($^\circ$)",
                ]
                if "correlated" in dic["VIS"]["AMPTYP"]:
                    ylabels[2] = "Raw correlated flux"
            else:
                ylabels = [
                    "Squared visibility",
                    "Visibility",
                    "VISAMP",
                    "VISPHI ($^\circ$)",
                    "Total flux (Jy)",
                    "Closure phase ($^\circ$)",
                ]
                if "correlated" in dic["VIS"]["AMPTYP"]:
                    ylabels[2] = "Correlated flux (Jy)"
            tags = [
                "squaredvisibility",
                "visibility",
                "visamp",
                "visphi",
                "flux",
                "t3phi",
            ]
            if "correlated" in dic["VIS"]["AMPTYP"]:
                tags[2] = "corrflux"
            if "CALIB_RAW_INT" in input_fitsfile:
                # append TF2 table
                oi_types.append("tf2")
                ylabels.append("TF2")
                tags.append("tf2")
        if "RAW_CPHASE" in input_fitsfile:
            oi_types = ["t3amp", "t3phi"]
            ylabels = ["T3AMP", "Raw closure phase ($^\circ$)"]
            tags = ["t3amp", "t3phi"]
        if "RAW_DPHASE" in input_fitsfile:
            oi_types = ["visamp", "visphi", "visamp"]
            ylabels = ["VISAMP", "VISPHI ($^\circ$)", "CFXAMP"]
            tags = ["visamp", "visphi", "cfxamp"]
        if "RAW_SPECTRUM" in input_fitsfile:
            oi_types = ["flux"]
            ylabels = ["Raw total flux"]
            tags = ["flux"]
        if "RAW_TF2" in input_fitsfile:
            oi_types = ["tf2"]
            ylabels = ["TF2"]
            tags = ["tf2"]
        if "RAW_VIS2" in input_fitsfile:
            oi_types = ["vis2"]
            ylabels = ["VIS2"]
            tags = ["squaredvisibility"]

        for k in range(len(tags)):
            fig, ((ax1)) = plt.subplots(
                1, 1, sharey=False, sharex=False, figsize=figsize
            )
            # print(k)
            # print(len(tags))
            # print(tags[k])
            # print(oi_types[k])
            # print(ylabels[k])
            return_val = make_plot_with_wavelength(
                dic,
                ax1,
                tags[k],
                oi_type=oi_types[k],
                ylabel=ylabels[k],
                verbose=verbose,
                annotate=annotate,
                wl_lim=wl_lim,
                figsize=figsize,
                input_fitsfile=input_fitsfile,
            )
            if return_val == 0:
                outputfig = (
                    outputdir
                    + "_".join(os.path.basename(input_fitsfile).split(".")[:-1])
                    + "_"
                    + tags[k]
                )
                if save_png:
                    plt.savefig(outputfig + ".png", dpi=200)
                if save_eps:
                    plt.savefig(outputfig + ".eps", format="eps", dpi=300)
            plt.close(fig)

        #####################################
        # plot images:
        #####################################
        imgs = []
        if "BSimag" in input_fitsfile or "BSreal" in input_fitsfile:
            imgs = [dic["IMG"]]
            if verbose:
                print("Plot images.")
        if "DSPtarget" in input_fitsfile or "fringePeak" in input_fitsfile:
            imgs = [np.arcsinh(dic["IMG"])]
            if verbose:
                print("Plot images.")
        for k in range(len(imgs)):
            img = imgs[k]
            xlabel = "x"
            ylabel = "y"
            fig, ((ax1)) = plt.subplots(
                1, 1, sharey=False, sharex=False, figsize=figsize
            )
            # ax1.plot(0, np.nan, "-", label=lab ,alpha=0,color='white')
            plt.imshow(img, origin="lower", aspect="auto")  # extent=extent
            # if math.isnan(xlim[0]):
            xlim = (None, None)
            plot_config(xlabel, ylabel, "", ax1, dic, xlim=xlim)
            outputfig = outputdir + "_".join(
                os.path.basename(input_fitsfile).split(".")[:-1]
            )
            if save_png:
                fig.savefig(outputfig + ".png", dpi=200)
            if save_eps:
                fig.savefig(outputfig + ".eps", format="eps", dpi=300)
            plt.close(fig)
        # if 'nrjImag' in input_fitsfile:
        # if 'nrjReal' in input_fitsfile:

        # plots as a function of time
        k = 0
        xs = []
        ys = []
        if "matis_eop" in input_fitsfile:
            xs = [dic["EOP"]["MJD"]] * 3
            ys = [[dic["EOP"]["PMX"]], [dic["EOP"]["PMY"]], [dic["EOP"]["DUT"]]]
            xlabels = ["MJD (d)"] * 3
            ylabels = ["PMX (arcsec)", "PMY (arcsec)", "DUT (s)"]
            tags = ["pmx", "pmy", "dut"]
        if "OI_OPDWVPO" in input_fitsfile:
            xs = [dic["OPD"]["MJD"]] * 3
            ys = [
                np.transpose(dic["OPD"]["OPD"]),
                np.transpose(dic["OPD"]["TEMPOFF"]),
                np.transpose(dic["OPD"]["HUMOFF"]),
            ]
            xlabels = ["MJD (d)"] * 3
            ylabels = ["OPD ($\mu$m)", "TEMPOFF", "HUMOFF"]
            tags = ["opd", "tempoff", "humoff"]
            sta_indices = [dic["OPD"]["STA_INDEX"]] * 3
        for k in range(len(xs)):
            x = xs[k]
            y = ys[k]
            if tags[k] == "opd" or tags[k] == "tempoff" or tags[k] == "humoff":
                sta_index = sta_indices[k][0]
            fig, ((ax1)) = plt.subplots(
                1, 1, sharey=False, sharex=False, figsize=figsize
            )
            for j in range(len(y)):
                if tags[k] == "opd" or tags[k] == "tempoff" or tags[k] == "humoff":
                    # print(sta_index)
                    sta_label = (
                        dic["STA_NAME"][sta_index[j * 2] == dic["STA_INDEX"]][0]
                        + "-"
                        + dic["STA_NAME"][sta_index[j * 2 + 1] == dic["STA_INDEX"]][0]
                    )
                    lab = r"%s" % (sta_label)
                    try:
                        title = (
                            dic["TARGET"]
                            + "\n"
                            + "date: "
                            + dic["DATE-OBS"]
                            + "\n"
                            + "TPL start: "
                            + dic["TPL_START"]
                            + "\n"
                            + dic["CATEGORY"]
                            + " "
                            + dic["BAND"]
                            + " "
                            + dic["DISPNAME"]
                            + " "
                            + dic["BCD1"]
                            + "-"
                            + dic["BCD2"]
                        )
                    except TypeError as e:
                        if verbose:
                            print("Unable to make the plot title. ")
                            if verbose:
                                print(e)
                        title = ""
                else:
                    lab = "_nolegend_"
                    title = ""
                plt.plot(x, y[j], label=lab)
            # if math.isnan(xlim[0]):
            xlim = (None, None)
            plot_config(
                xlabels[k], ylabels[k], title, ax1, dic, xlim=xlim, annotate=annotate
            )
            outputfig = (
                outputdir
                + "_".join(os.path.basename(input_fitsfile).split(".")[:-1])
                + "_"
                + tags[k]
            )
            if save_png:
                fig.savefig(outputfig + ".png", dpi=200)
            if save_eps:
                fig.savefig(outputfig + ".eps", format="eps", dpi=300)
            plt.close(fig)
    plt.close("all")


###############################################################################################################
# make_mosaic_plot
###############################################################################################################
# default wavelengths: in L band: the central wavelength from the observation, bandwidth = 0.2 um
#                      in N band: 10.7 um, bandwidth = 0.2 um
# input_fitsfile_list = [incoherent_CAL_INT,coherent_FLUXCAL_int]
# oi_types_list = [ ['uv','vis','t3phi'], ['visamp'] ]
# use only for oifits files
def make_mosaic_plot(
    input_fitsfile_list,
    output_file_path,
    oi_types_list=[
        ["uv", "vis_pbl", "visamp_pbl", "flux", "vis", "visamp", "visphi", "t3phi"]
    ],
    verbose=False,
    save_png=True,
    save_eps=False,
    sel_wl=np.nan,
    bandwidth=np.nan,
    annotate=False,
    legend_loc="upper right",
    legend_loc_pbl="lower left",
    wl_lim=(np.nan, np.nan),
    B_lim=(np.nan, np.nan),
    fit_model=False,
    figsize=(15, 15),
    ext=0,
    plot_errorbars_wlplot="fill",
    err_show_ratio=0.04,
):
    if verbose:
        print("Make mosaic plot")
    outfig, axes = plt.subplots(3, 3, sharey=False, sharex=False, figsize=figsize)

    dic_layout = {
        "uv": [0, 1],
        "flux": [0, 2],
        "t3phi": [1, 0],
        "vis": [1, 1],
        "visamp": [1, 2],
        "visphi": [2, 0],
        "vis_pbl": [2, 1],
        "visamp_pbl": [2, 2],
    }

    figsize = (figsize[0] / 3.0, figsize[1] / 3.0)
    for i in range(len(input_fitsfile_list)):
        input_fitsfile = input_fitsfile_list[i]
        print(input_fitsfile)

        # open fits file
        # check if oifits file is valid
        oifits_valid = check_oifits_valid(input_fitsfile)
        if not oifits_valid:
            print(input_fitsfile + " not a valid OIFITS file.")
            continue
        dic = open_fits(input_fitsfile, verbose, ext=ext)
        if dic:
            # print(dic['DISPNAME'])
            for oi_type in oi_types_list[i]:
                if oi_type == "uv":
                    return_val = make_uv_plot(
                        dic,
                        axes[dic_layout["uv"][0], dic_layout["uv"][1]],
                        verbose=verbose,
                        annotate=annotate,
                        B_lim=B_lim,
                        figsize=figsize,
                    )

                    # if ('AVGCAL_INT' in input_fitsfile) or ('AVGFLUXCAL_INT' in input_fitsfile) or \
                    #    ('BCDCAL_INT' in input_fitsfile) or ('FINALCAL_INT' in input_fitsfile):
                    #    plot_title = dic['TARGET'] + "\n" + \
                    #        "date: " + dic['DATE-OBS'] + "\n" + \
                    #        dic['CATEGORY'] + ' ' + dic['BAND'] + ' ' + dic['DISPNAME'] + '\n' + \
                    #        dic['STA1'] + '-' + dic['STA2'] + '-' + dic['STA3'] + '-' + dic['STA4']
                    # else:
                    # print('d'+dic['DISPNAME'] +'/')
                    chop_str = ""
                    if "L" in dic["BAND"] or "M" in dic["BAND"]:
                        if dic["CHOP_ST"] == "T":
                            chop_str = "chopped"
                        if dic["CHOP_ST"] == "F":
                            chop_str = "nochop"
                    plot_title = (
                        dic["TARGET"]
                        + "\n"
                        + "date: "
                        + dic["DATE-OBS"]
                        + "\n"
                        + "TPL start: "
                        + dic["TPL_START"]
                        + "\n"
                        + dic["CATEGORY"]
                        + " "
                        + dic["BAND"]
                        + " "
                        + dic["DISPNAME"]
                        + " "
                        + dic["BCD1"]
                        + "-"
                        + dic["BCD2"]
                        + " "
                        + chop_str
                        + "\n"
                        + dic["STA1"]
                        + "-"
                        + dic["STA2"]
                        + "-"
                        + dic["STA3"]
                        + "-"
                        + dic["STA4"]
                    )
                    axes[0, 0].annotate(
                        plot_title,
                        xy=(0, 1),
                        xytext=(12, -12),
                        va="top",
                        xycoords="axes fraction",
                        textcoords="offset points",
                        fontsize=16,
                    )
                    plot_config(
                        "",
                        "",
                        "",
                        axes[0, 0],
                        dic,
                        plot_legend=False,
                        annotate=True,
                        annotate_va="bottom",
                        annotate_fontsize=10,
                        annotate_xy=(0, 0),
                    )
                    axes[0, 0].axis("off")
                else:
                    if oi_type == "vis_pbl" or oi_type == "visamp_pbl":
                        if oi_type == "vis_pbl":
                            if ("CAL_INT" not in input_fitsfile) or (
                                "RAW_INT" in dic["PRO_CATG"]
                            ):
                                ylabel = "Raw visibility"
                            else:
                                ylabel = "Visibility"
                            tag = "visibility_vs_pbl"
                        if oi_type == "visamp_pbl":
                            if "correlated" in dic["VIS"]["AMPTYP"]:
                                if ("CAL_INT" not in input_fitsfile) or (
                                    "RAW_INT" in dic["PRO_CATG"]
                                ):
                                    ylabel = "Raw correlated flux"
                                else:
                                    ylabel = "Correlated flux (Jy)"
                                tag = "corrflux_vs_pbl"
                            else:
                                if ("CAL_INT" not in input_fitsfile) or (
                                    "RAW_INT" in dic["PRO_CATG"]
                                ):
                                    ylabel = "Raw VISAMP"
                                else:
                                    ylabel = "VISAMP"
                                tag = "visamp_vs_pbl"
                        if sel_wl:  # if it is not None
                            if math.isnan(sel_wl):
                                print(dic["DISPNAME_L"])
                                print(dic["DETNAME"])
                                if "L" in dic["DETNAME"]:
                                    if "LOW" in dic["DISPNAME"]:
                                        sel_wl = 3.4
                                    else:
                                        # print (np.nanmin(dic['WLEN']),np.nanmax(dic['WLEN']))
                                        if (
                                            np.nanmin(dic["WLEN"]) < 3.2e-6
                                            and np.nanmax(dic["WLEN"]) > 3.6e-6
                                        ):
                                            sel_wl = 3.4
                                        elif (
                                            np.nanmin(dic["WLEN"]) < 3.58e-6
                                            and np.nanmax(dic["WLEN"]) > 3.9e-6
                                        ):
                                            sel_wl = 3.7
                                        else:
                                            sel_wl = dic["WL_CENTRAL"]
                                    bandwidth = 0.2
                                if "N" in dic["DETNAME"]:
                                    sel_wl = 8.8
                                    bandwidth = 0.8
                        else:
                            if "L" in dic["DETNAME"]:
                                if "LOW" in dic["DISPNAME"]:
                                    sel_wl = 3.4
                                else:
                                    if (
                                        np.nanmin(dic["WLEN"]) < 3.2e-6
                                        and np.nanmax(dic["WLEN"]) > 3.6e-6
                                    ):
                                        sel_wl = 3.4
                                    elif (
                                        np.nanmin(dic["WLEN"]) < 3.58e-6
                                        and np.nanmax(dic["WLEN"]) > 3.9e-6
                                    ):
                                        sel_wl = 3.7
                                    else:
                                        sel_wl = dic["WL_CENTRAL"]
                                bandwidth = 0.2
                            if "N" in dic["DETNAME"]:
                                sel_wl = 8.8
                                bandwidth = 0.8
                        # print('vis visamp')
                        return_val = make_plot_with_baseline(
                            dic,
                            axes[dic_layout[oi_type][0], dic_layout[oi_type][1]],
                            tag,
                            oi_type=oi_type.replace("_pbl", ""),
                            ylabel=ylabel,
                            annotate=False,
                            annotate_baselines=True,
                            legend_loc=legend_loc_pbl,
                            B_lim=B_lim,
                            verbose=False,
                            sel_wl=sel_wl,
                            bandwidth=bandwidth,
                            figsize=figsize,
                            fit_model=fit_model,
                            input_fitsfile=input_fitsfile,
                        )
                        axes[dic_layout[oi_type][0], dic_layout[oi_type][1]].annotate(
                            r"$\lambda = %.2f\ \mu$m" % sel_wl,
                            xy=(0, -0.01),
                            xytext=(12, -12),
                            va="top",
                            xycoords="axes fraction",
                            textcoords="offset points",
                            fontsize=12,
                        )
                        # print(sel_wl)
                    else:
                        if oi_type == "flux":
                            if ("CAL_INT" not in input_fitsfile) or (
                                "RAW_INT" in dic["PRO_CATG"]
                            ):
                                ylabel = "Raw total flux"
                            else:
                                ylabel = "Total flux (Jy)"
                            tag = "flux"
                        if oi_type == "vis":
                            if ("CAL_INT" not in input_fitsfile) or (
                                "RAW_INT" in dic["PRO_CATG"]
                            ):
                                ylabel = "Raw visibility"
                            else:
                                ylabel = "Visibility"
                            tag = "visibility"
                        if oi_type == "visamp":
                            if "correlated" in dic["VIS"]["AMPTYP"]:
                                if ("CAL_INT" not in input_fitsfile) or (
                                    "RAW_INT" in dic["PRO_CATG"]
                                ):
                                    ylabel = "Raw correlated flux"
                                else:
                                    ylabel = "Correlated flux (Jy)"
                                tag = "corrflux"
                            else:
                                if ("CAL_INT" not in input_fitsfile) or (
                                    "RAW_INT" in dic["PRO_CATG"]
                                ):
                                    ylabel = "Raw VISAMP"
                                else:
                                    ylabel = "VISAMP"
                                tag = "visamp"
                        if oi_type == "visphi":
                            if ("CAL_INT" not in input_fitsfile) or (
                                "RAW_INT" in dic["PRO_CATG"]
                            ):
                                ylabel = "Raw differential phase ($^\circ$)"
                            else:
                                ylabel = "Differential phase ($^\circ$)"
                            tag = "visphi"
                        if oi_type == "t3phi":
                            if ("CAL_INT" not in input_fitsfile) or (
                                "RAW_INT" in dic["PRO_CATG"]
                            ):
                                ylabel = "Raw closure phase ($^\circ$)"
                            else:
                                ylabel = "Closure phase ($^\circ$)"
                            tag = "t3phi"
                        if oi_type == "t3amp":
                            if ("CAL_INT" not in input_fitsfile) or (
                                "RAW_INT" in dic["PRO_CATG"]
                            ):
                                ylabel = "Raw T3AMP"
                            else:
                                ylabel = "T3AMP"
                            tag = "t3amp"

                        # print('vmi')
                        return_val = make_plot_with_wavelength(
                            dic,
                            axes[dic_layout[oi_type][0], dic_layout[oi_type][1]],
                            tag,
                            oi_type=oi_type,
                            ylabel=ylabel,
                            verbose=verbose,
                            legend_loc=legend_loc,
                            annotate=annotate,
                            wl_lim=wl_lim,
                            figsize=figsize,
                            input_fitsfile=input_fitsfile,
                            plot_errorbars=plot_errorbars_wlplot,
                            err_show_ratio=err_show_ratio,
                        )
    # print('vmi2')
    if save_png:
        plt.savefig(output_file_path + ".png", dpi=200)
    if save_eps:
        plt.savefig(output_file_path + ".eps", format="eps", dpi=300)
    plt.close(outfig)
    plt.close("all")
    # print('vmi3')


# input_fitsfile_types: list of strings: 'c': take only correlated flux, 'f': take only total flux
def show_corr_total_flux(
    input_fitsfiles,
    input_fitsfile_types,
    outputdir="plots",
    fn_pattern="",
    verbose=False,
    save_png=True,
    save_eps=False,
    sel_wl=3.6,
    bandwidth=0.2,
    file_type="",
    pro_catg="",
    annotate=False,
    wl_lim=(np.nan, np.nan),
    B_lim=(np.nan, np.nan),
):
    # check if output directory exists
    # if not, create it
    if input_fitsfiles != []:
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)

    #####################################
    # plot as a function of baseline
    #####################################
    if verbose:
        print("Plots as a function of baseline.")
    k = 0
    xs = []
    ys = []
    yerrs = []
    us = []
    vs = []
    pbls = []
    pblas = []
    mjds = []
    sta_indices = []
    sta_labels = []
    xlabels = []
    ylabels = []
    data_types = []
    tags = []
    xlabels.append("Projected baseline (m)")  #'$B_\mathrm{p}$ (m)')
    if "CAL_INT" not in input_fitsfiles[0]:
        ylabels.append("Raw correlated flux")
    else:
        ylabels.append("Correlated flux (Jy)")
    tags.append("corr_total_flux_vs_pbl")

    for i in range(len(input_fitsfiles)):
        input_fitsfile = input_fitsfiles[i]
        input_fitsfile_type = input_fitsfile_types[i]
        print(input_fitsfile)
        # open fits file
        dic = {}
        if "OBJ_CORR_FLUX" not in input_fitsfile and "CALIB_CAL" not in input_fitsfile:
            dic = open_fits(input_fitsfile, verbose)
        if dic:
            if not (pro_catg == ""):
                if not (dic["PRO_CATG"] == pro_catg):
                    return

            if input_fitsfile_type == "c":
                # if ('RAW_INT' in input_fitsfile) or ('CAL_INT' in input_fitsfile) or ('oifits' in input_fitsfile) or (
                #         'OIFITS' in input_fitsfile or file_type == 'oifits'):

                # xs = [dic['WLEN'] * 1e6] * 2  # (um)
                # try:
                #     ys = [dic['VIS2']['VIS2'], np.sqrt(dic['VIS2']['VIS2']) ]
                #     yerrs = [dic['VIS2']['VIS2ERR'], np.abs(0.5*dic['VIS2']['VIS2ERR']/np.sqrt(dic['VIS2']['VIS2']))]
                #     us = [dic['VIS2']['U'], dic['VIS2']['U']]
                #     vs = [dic['VIS2']['V'], dic['VIS2']['V']]
                #     mjds = [dic['VIS2']['MJD'], dic['VIS2']['MJD']]
                #     sta_indices = [dic['VIS2']['STA_INDEX'], dic['VIS2']['STA_INDEX']]
                # except KeyError as e:
                #     if verbose: print(e)
                #     ys = [0.0]
                #     yerrs = [0.0]
                #     us = [0.0]
                #     vs = [0.0]
                #     mjds = [0.0]
                #     sta_indices = []
                # xlabels = ['$B_\mathrm{p}$ (m)'] * 2
                # if ('CAL_INT' not in input_fitsfile):
                #     ylabels = ['Raw squared visibility', 'Raw visibility']
                # else:
                #     ylabels = ['Squared visibility', 'Visibility']
                # tags = ['squaredvisibility_vs_pbl', 'visibility_vs_pbl']

                # if 'correlated' in dic['VIS']['AMPTYP']:
                for j in range(len(dic["VIS"]["VISAMP"])):
                    xs.append(dic["WLEN"] * 1e6)
                    ys.append(dic["VIS"]["VISAMP"][j])
                    yerrs.append(dic["VIS"]["VISAMPERR"][j])
                    us.append(dic["VIS"]["U"][j])
                    vs.append(dic["VIS"]["V"][j])
                    mjds.append(dic["VIS"]["MJD"][j])
                    sta_indices.append(dic["VIS"]["STA_INDEX"][j])
                    pbls.append(np.sqrt(us[-1] ** 2 + vs[-1] ** 2))
                    pbla = np.arctan2(us[-1], vs[-1]) * 180.0 / np.pi
                    if pbla < 0.0:
                        pbla = pbla + 180.0
                    pblas.append(pbla)
                    data_types.append(input_fitsfile_type)
                    sta_labels.append(
                        dic["STA_NAME"][
                            dic["VIS"]["STA_INDEX"][j, 0] == dic["STA_INDEX"]
                        ][0]
                        + "-"
                        + dic["STA_NAME"][
                            dic["VIS"]["STA_INDEX"][j, 1] == dic["STA_INDEX"]
                        ][0]
                    )

            if input_fitsfile_type == "f":
                for j in range(len(dic["FLUX"]["FLUX"])):
                    xs.append(dic["WLEN"] * 1e6)  # (um)
                    ys.append(dic["FLUX"]["FLUX"][j])
                    yerrs.append(dic["FLUX"]["FLUXERR"][j])
                    us.append(0.0)
                    vs.append(0.0)
                    mjds.append(dic["FLUX"]["MJD"][j])
                    sta_indices.append(dic["FLUX"]["STA_INDEX"][j])
                    pbls.append(0.0)
                    pblas.append(0.0)
                    data_types.append(input_fitsfile_type)
                    sta_labels.append("")
                    # sta_labels.append(dic['STA_NAME'][dic['FLUX']['STA_INDEX'][j] == dic['STA_INDEX']][0] + '-' + \
                    # dic['STA_NAME'][dic['FLUX']['STA_INDEX'][j] == dic['STA_INDEX']][0])

    ymin = []
    ymax = []
    xmin = []
    xmax = []
    y_new = []
    yerr_new = []
    if (sel_wl % 1) == 0:
        lambda_tag = "%dum" % sel_wl
    else:
        int_part = np.floor(sel_wl)
        fract_part = sel_wl - 1.0 * int_part
        lambda_tag = "%dum%d" % (int_part, np.round(fract_part * 10.0))

    for k in range(len(xs)):
        x = xs[k]
        y = ys[k]
        pbl = pbls[k]
        pbla = pblas[k]
        yerr = yerrs[k]
        mjd = mjds[k]
        sta_index = sta_indices[k]
        sta_label = sta_labels[k]

        # if 'squaredvisibility' in tags[k] or 'visibility' in tags[k]  or 'visamp' in tags[k] or 'tf2' in tags[
        #     k] or 't3amp' in  tags[k]:
        #     ymin = [0.0]
        # ymax = [1.0]
        ymin = [0.0]
        idx = np.logical_and(x > sel_wl - bandwidth / 2.0, x < sel_wl + bandwidth / 2.0)

        y_new.append(np.nanmedian(y[idx]))
        yerr_new.append(np.nanmean(yerr[idx]))

        pmin = np.nanpercentile(y_new, 10.0)
        pmax = np.nanpercentile(y_new, 95.0)
        pmid = (pmin + pmax) / 2.0
        pdiff = pmax - pmin
        lim_fact = 1.75
        ymin.append(pmid - pdiff / 2.0 * lim_fact)
        ymax.append(1.05 * np.amax(y_new))  # (pmid + pdiff / 2.0 * lim_fact)
        xmax.append(np.amax(pbl))

    if not np.all(np.isnan(y_new)):
        fig, ((ax1)) = plt.subplots(
            1, 1, sharey=False, sharex=False, figsize=(3.5, 3.5)
        )
        plt.set_cmap("hsv")

        plt.errorbar(
            pbls,
            y_new,
            yerr=yerr_new,
            fmt=".",
            ecolor="gray",
            alpha=0.6,
            zorder=0,
            capsize=3,
        )
        if annotate == True:
            plt.scatter(
                pbls,
                y_new,
                c=pblas,
                vmin=0.0,
                vmax=180.0,
                marker="o",
                edgecolor="black",
                zorder=5,
            )  # c=pbla
        else:
            plt.scatter(
                pbls,
                y_new,
                c="blue",
                vmin=0.0,
                vmax=180.0,
                marker="o",
                edgecolor="black",
                zorder=5,
            )  # c=pbla

        if annotate == True:
            for k in range(len(xs)):
                plt.text(
                    pbls[k] - 1.0,
                    y_new[k],
                    sta_labels[k],
                    color="gray",
                    alpha=0.6,
                    ha="right",
                )

        plot_title = dic["TARGET"] + " @%.2f $\mu$m" % sel_wl
        # + "\n" + "date: " + dic['DATE-OBS'] + "\n" + dic['CATEGORY'] + ' ' + \
        #    dic['BAND'] + ' ' + dic['DISPNAME'] + ' ' + dic['BCD1'] + '-' + dic['BCD2'] + ', ' + '%.2f $\mu$m'%sel_wl
        # left, right = plt.xlim()
        if annotate == True:
            cbar = plt.colorbar(orientation="vertical")
            cbar.ax.set_ylabel("Baseline PA ($^\circ$)")
        if math.isnan(B_lim[0]):
            xlim = (-1.5, np.nanmax(xmax) * 1.1)
        else:
            xlim = B_lim
        plot_config(
            xlabels[0],
            ylabels[0],
            plot_title,
            ax1,
            dic,
            ylim=(np.nanmin(ymin), np.nanmax(ymax)),
            xlim=xlim,
            plot_legend=False,
            annotate=False,
        )
        outputfig = (
            outputdir
            + "_".join(os.path.basename(input_fitsfile).split(".")[:-1])
            + "_"
            + tags[0]
            + "_"
            + lambda_tag
        )
        if save_png:
            fig.savefig(outputfig + ".png", dpi=200)
        if save_eps:
            fig.savefig(outputfig + ".eps", format="eps", dpi=300)
        plt.close(fig)

    #####################################
    # plot as a function of wavelength: total and correlated spectra
    #####################################
    if verbose:
        print("Plots as a function of wavelength.")
    k = 0
    xs = []
    ys = []
    yerrs = []
    us = []
    vs = []
    pbls = []
    pblas = []
    mjds = []
    sta_indices = []
    sta_labels = []
    xlabels = []
    ylabels = []
    data_types = []
    tags = []
    xlabels.append("Wavelength ($\mu$m)")  #'$B_\mathrm{p}$ (m)')
    if "CAL_INT" not in input_fitsfiles[0]:
        ylabels.append("Raw correlated flux")
    else:
        ylabels.append("Flux (Jy)")
    tags.append("corr_total_flux")

    for i in range(len(input_fitsfiles)):
        input_fitsfile = input_fitsfiles[i]
        input_fitsfile_type = input_fitsfile_types[i]
        print(input_fitsfile)
        # open fits file
        dic = {}
        if "OBJ_CORR_FLUX" not in input_fitsfile and "CALIB_CAL" not in input_fitsfile:
            dic = open_fits(input_fitsfile, verbose)
        if dic:
            if not (pro_catg == ""):
                if not (dic["PRO_CATG"] == pro_catg):
                    return

            if input_fitsfile_type == "c":
                # if 'correlated' in dic['VIS']['AMPTYP']:
                xs.append(dic["WLEN"] * 1e6)
                ys.append(dic["VIS"]["VISAMP"])
                yerrs.append(dic["VIS"]["VISAMPERR"])
                us.append(dic["VIS"]["U"])
                vs.append(dic["VIS"]["V"])
                mjds.append(dic["VIS"]["MJD"])
                sta_indices.append(dic["VIS"]["STA_INDEX"])
                pbls.append(np.sqrt(us[-1] ** 2 + vs[-1] ** 2))
                pblas.append(np.arctan2(us[-1], vs[-1]) * 180.0 / np.pi)
                idx = np.where(pblas[-1] < 0.0)
                pblas[-1][idx] = pblas[-1][idx] + 180.0
                data_types.append(input_fitsfile_type)
                sta_label = []
                for j in range(len(dic["VIS"]["STA_INDEX"])):
                    sta_label.append(
                        dic["STA_NAME"][
                            dic["VIS"]["STA_INDEX"][j, 0] == dic["STA_INDEX"]
                        ][0]
                        + "-"
                        + dic["STA_NAME"][
                            dic["VIS"]["STA_INDEX"][j, 1] == dic["STA_INDEX"]
                        ][0]
                    )
                sta_labels.append([sta_label])

            if input_fitsfile_type == "f":
                xs.append(dic["WLEN"] * 1e6)  # (um)
                ys.append(dic["FLUX"]["FLUX"])
                yerrs.append(dic["FLUX"]["FLUXERR"])
                us.append([0.0])
                vs.append([0.0])
                mjds.append(dic["FLUX"]["MJD"])
                sta_indices.append(dic["FLUX"]["STA_INDEX"])
                pbls.append([0.0])
                pblas.append([0.0])
                data_types.append(input_fitsfile_type)
                sta_label = []
                for j in range(len(dic["FLUX"]["STA_INDEX"])):
                    sta_label.append(
                        dic["STA_NAME"][
                            dic["FLUX"]["STA_INDEX"][j] == dic["STA_INDEX"]
                        ][0]
                        + "-"
                        + dic["STA_NAME"][
                            dic["FLUX"]["STA_INDEX"][j] == dic["STA_INDEX"]
                        ][0]
                    )
                sta_labels.append([sta_label])

    fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(3.5, 3.5))
    ymin = []
    ymax = []
    xmin = []
    xmax = []
    for k in range(len(xs)):
        x = xs[k]
        y = ys[k]
        pbl = pbls[k]
        pbla = pblas[k]
        yerr = yerrs[k]
        mjd = mjds[k]
        sta_index = sta_indices[k]
        # sta_label = sta_labels[k]
        N = len(mjd)

        # if 'squaredvisibility' in tags[k] or 'visibility' in tags[k]  or 'visamp' in tags[k] or 'tf2' in tags[
        #     k] or 't3amp' in  tags[k]:
        #     ymin = [0.0]
        # ymax = [1.0]
        ymin = [0.0]

        sta_labels = []
        for j in range(len(y)):
            y_new.append(np.nanmedian(y[j][idx]))
            yerr_new.append(np.nanmean(yerr[j][idx]))
            if data_types[k] == "c":
                # lab = r'%s' % (sta_labels[-1])
                # lab = r'%s $B_\mathrm{p}=%.1f$ m, $\phi = %.1f^{\circ}$' % (sta_label[j], pbl[j], pbla[j])
                lab = r"$B_\mathrm{p}=%.0f$ m" % (pbl[j])
            if data_types[k] == "f":
                lab = "Total spectrum"
            if data_types[k] == "c":
                plt.plot(x, y[j], label=lab)
            if data_types[k] == "f":
                plt.plot(x, y[j], "--", label=lab)
            # for j in range(len(y)):
            #     plt.text(pbl[j]-1.0, y_new[j], sta_labels[j], color='gray', alpha=0.6,ha='right')
            ynew2 = y[j]
            # exclude M band (4.1-4.6 um) for the calculation of y plot limits
            M_idx = np.logical_and(x > 4.1, x < 4.6)
            ynew2[M_idx] = np.nan
            pmin = np.nanpercentile(ynew2, 10.0)
            pmax = np.nanpercentile(ynew2, 95.0)
            pmid = (pmin + pmax) / 2.0
            pdiff = pmax - pmin
            lim_fact = 1.75
            ymin.append(pmid - pdiff / 2.0 * lim_fact)
            ymax.append(1.05 * np.nanmax(ynew2))  # (pmid + pdiff / 2.0 * lim_fact)
            # xmax.append(np.amax(pbl))
    plot_title = dic["TARGET"]
    # + "\n" + "date: " + dic['DATE-OBS'] + "\n" + dic['CATEGORY'] + ' ' + \
    #    dic['BAND'] + ' ' + dic['DISPNAME'] + ' ' + dic['BCD1'] + '-' + dic['BCD2'] + ', ' + '%.1f $\mu$m'%sel_wl
    # left, right = plt.xlim()
    # cbar = plt.colorbar(orientation='vertical')
    # cbar.ax.set_ylabel('Baseline PA ($^\circ$)')
    # print((np.nanmin(ymin), np.nanmax(ymax)))
    # print((0.0,right))
    wl_lim_temp = wl_lim
    if wl_lim_temp[0]:
        if math.isnan(wl_lim_temp[0]):
            wl_lim_temp = (np.amin(x), np.amax(x))
    plot_config(
        xlabels[0],
        ylabels[0],
        plot_title,
        ax1,
        dic,
        ylim=(np.nanmin(ymin), np.nanmax(ymax)),
        xlim=wl_lim_temp,
        plot_legend=True,
        annotate=False,
    )
    outputfig = (
        outputdir
        + "_".join(os.path.basename(input_fitsfile).split(".")[:-1])
        + "_"
        + tags[0]
    )
    if save_png:
        fig.savefig(outputfig + ".png", dpi=200)
    if save_eps:
        fig.savefig(outputfig + ".eps", format="eps", dpi=300)
    plt.close(fig)

    #####################################
    # plot as a function of wavelength: visibilities calculated from total and correlated spectra
    #####################################
    if verbose:
        print("Plots as a function of wavelength.")
    k = 0
    xs = []
    ys = []
    yerrs = []
    us = []
    vs = []
    pbls = []
    pblas = []
    mjds = []
    sta_indices = []
    # sta_labels = []
    xlabels = []
    ylabels = []
    data_types = []
    tags = []
    total_flux = []
    xlabels.append("Wavelength ($\mu$m)")  #'$B_\mathrm{p}$ (m)')
    if "CAL_INT" not in input_fitsfiles[0]:
        ylabels.append("Raw visibility")
    else:
        ylabels.append("Visibility")
    tags.append("coh_vis")

    for i in range(len(input_fitsfiles)):
        input_fitsfile = input_fitsfiles[i]
        input_fitsfile_type = input_fitsfile_types[i]
        print(input_fitsfile)
        # open fits file
        dic = {}
        if "OBJ_CORR_FLUX" not in input_fitsfile and "CALIB_CAL" not in input_fitsfile:
            dic = open_fits(input_fitsfile, verbose)
        if dic:
            if not (pro_catg == ""):
                if not (dic["PRO_CATG"] == pro_catg):
                    return

            if input_fitsfile_type == "c":
                # if 'correlated' in dic['VIS']['AMPTYP']:
                xs.append(dic["WLEN"] * 1e6)
                ys.append(dic["VIS"]["VISAMP"])
                yerrs.append(dic["VIS"]["VISAMPERR"])
                us.append(dic["VIS"]["U"])
                vs.append(dic["VIS"]["V"])
                mjds.append(dic["VIS"]["MJD"])
                sta_indices.append(dic["VIS"]["STA_INDEX"])
                pbls.append(np.sqrt(us[-1] ** 2 + vs[-1] ** 2))
                pblas.append(np.arctan2(us[-1], vs[-1]) * 180.0 / np.pi)
                idx = np.where(pblas[-1] < 0.0)
                pblas[-1][idx] = pblas[-1][idx] + 180.0
                data_types.append(input_fitsfile_type)
                sta_label = []
                for j in range(len(dic["VIS"]["STA_INDEX"])):
                    sta_label.append(
                        dic["STA_NAME"][
                            dic["VIS"]["STA_INDEX"][j, 0] == dic["STA_INDEX"]
                        ][0]
                        + "-"
                        + dic["STA_NAME"][
                            dic["VIS"]["STA_INDEX"][j, 1] == dic["STA_INDEX"]
                        ][0]
                    )
                sta_labels.append([sta_label])

            if input_fitsfile_type == "f":
                xs.append(dic["WLEN"] * 1e6)  # (um)
                ys.append(dic["FLUX"]["FLUX"])
                yerrs.append(dic["FLUX"]["FLUXERR"])
                us.append([0.0])
                vs.append([0.0])
                mjds.append(dic["FLUX"]["MJD"])
                sta_indices.append(dic["FLUX"]["STA_INDEX"])
                pbls.append([0.0])
                pblas.append([0.0])
                data_types.append(input_fitsfile_type)
                sta_label = []
                for j in range(len(dic["FLUX"]["STA_INDEX"])):
                    sta_label.append(
                        dic["STA_NAME"][
                            dic["FLUX"]["STA_INDEX"][j] == dic["STA_INDEX"]
                        ][0]
                        + "-"
                        + dic["STA_NAME"][
                            dic["FLUX"]["STA_INDEX"][j] == dic["STA_INDEX"]
                        ][0]
                    )
                sta_labels.append([sta_label])
                total_flux.append(dic["FLUX"]["FLUX"][0])

    total_flux = np.array(total_flux)
    total_flux = np.nanmean(total_flux, axis=0)
    for i in range(len(ys)):
        if data_types[i] == "c":
            ys[i] = ys[i] / total_flux
            yerrs[i] = yerrs[i] / total_flux

    fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(3.5, 3.5))
    ymin = []
    ymax = []
    xmin = []
    xmax = []
    for k in range(len(xs)):
        x = xs[k]
        y = ys[k]
        pbl = pbls[k]
        pbla = pblas[k]
        yerr = yerrs[k]
        mjd = mjds[k]
        sta_index = sta_indices[k]
        # sta_label = sta_labels[k]
        N = len(mjd)

        # if 'squaredvisibility' in tags[k] or 'visibility' in tags[k]  or 'visamp' in tags[k] or 'tf2' in tags[
        #     k] or 't3amp' in  tags[k]:
        #     ymin = [0.0]
        # ymax = [1.0]
        ymin = [0.0]

        sta_labels = []
        for j in range(len(y)):
            y_new.append(np.nanmedian(y[j][idx]))
            yerr_new.append(np.nanmean(yerr[j][idx]))
            if data_types[k] == "c":
                # lab = r'%s' % (sta_label[j])
                # lab = r'%s $B_\mathrm{p}=%.1f$ m, $\phi = %.1f^{\circ}$' % (sta_label[j], pbl[j], pbla[j])
                lab = r"$B_\mathrm{p}=%.0f$ m" % (pbl[j])
            if data_types[k] == "f":
                lab = "Total spectrum"
            if data_types[k] == "c":
                plt.plot(x, y[j], label=lab)
            # if data_types[k] == 'f':
            #     plt.plot(x, y[j],'--', label=lab)
            # for j in range(len(y)):
            #     plt.text(pbl[j]-1.0, y_new[j], sta_labels[j], color='gray', alpha=0.6,ha='right')
            pmin = np.nanpercentile(y[j], 10.0)
            pmax = np.nanpercentile(y[j], 95.0)
            pmid = (pmin + pmax) / 2.0
            pdiff = pmax - pmin
            lim_fact = 1.75
            ymin.append(pmid - pdiff / 2.0 * lim_fact)
            ymax.append(1.05 * np.amax(y_new))  # (pmid + pdiff / 2.0 * lim_fact)
            # xmax.append(np.amax(pbl))
    plot_title = dic["TARGET"]
    # + "\n" + "date: " + dic['DATE-OBS'] + "\n" + dic['CATEGORY'] + ' ' + \
    #    dic['BAND'] + ' ' + dic['DISPNAME'] + ' ' + dic['BCD1'] + '-' + dic['BCD2'] + ', ' + '%.1f $\mu$m'%sel_wl
    # left, right = plt.xlim()
    # cbar = plt.colorbar(orientation='vertical')
    # cbar.ax.set_ylabel('Baseline PA ($^\circ$)')
    # print((np.nanmin(ymin), np.nanmax(ymax)))
    # print((0.0,right))
    wl_lim_temp = wl_lim
    if wl_lim_temp[0]:
        if math.isnan(wl_lim_temp[0]):
            wl_lim_temp = (np.amin(x), np.amax(x))
    plot_config(
        xlabels[0],
        ylabels[0],
        plot_title,
        ax1,
        dic,
        ylim=(0.0, 1.0),
        xlim=wl_lim_temp,
        plot_legend=True,
        legend_loc="best",
        annotate=False,
    )  # ylim=(np.nanmin(ymin), np.nanmax(ymax))
    outputfig = (
        outputdir
        + "_".join(os.path.basename(input_fitsfile).split(".")[:-1])
        + "_"
        + tags[0]
    )
    if save_png:
        fig.savefig(outputfig + ".png", dpi=200)
    if save_eps:
        fig.savefig(outputfig + ".eps", format="eps", dpi=300)
    plt.close(fig)


def create_obs_lst(
    rawdir,
    write_output_file=False,
    outfile_path="./fits_cat.txt",
    concat_tpl_starts=True,
    filter_terms=[],
):
    keys = []
    keys.append(
        {
            "key": "HIERARCH ESO TPL START",
            "name": "tpl_start",
            "format": "%-19s",
            "nameformat": "%-19s",
        }
    )
    # keys.append({'key': 'HIERARCH ESO OBS TARG NAME','name':'obs_targ_name','format':'"%-26s"', 'nameformat':'%-26s'})
    # keys.append({'key': 'HIERARCH ESO OBS NAME','name':'obs_name','format':'"%-32s"', 'nameformat':'%-32s'})
    keys.append(
        {
            "key": "HIERARCH ESO OBS TARG NAME",
            "name": "obs_targ_name",
            "format": "%-26s",
            "nameformat": "%-26s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO OBS NAME",
            "name": "obs_name",
            "format": "%-32s",
            "nameformat": "%-32s",
        }
    )
    keys.append(
        {
            "key": "DATE-OBS",
            "name": "date_obs",
            "format": "%-.24s",
            "nameformat": "%-24s",
        }
    )
    keys.append({"key": "RA", "name": "ra", "format": "%9.5f", "nameformat": "%-9s"})
    keys.append({"key": "DEC", "name": "dec", "format": "%9.5f", "nameformat": "%-9s"})
    # keys.append({'key': 'HIERARCH ESO INS BCD1 ID',  'name':'BCD1',         'format':'%-4s' , 'nameformat':'%-4s' })
    # keys.append({'key': 'HIERARCH ESO INS BCD2 ID',  'name':'BCD2',         'format':'%-4s' , 'nameformat':'%-4s' })
    keys.append(
        {
            "key": "HIERARCH ESO DET NAME",
            "name": "det_name",
            "format": "%-10s",
            "nameformat": "%-10s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO INS DIL NAME",
            "name": "dil_name",
            "format": "%-10s",
            "nameformat": "%-10s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO INS DIN NAME",
            "name": "din_name",
            "format": "%-10s",
            "nameformat": "%-10s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO DET SEQ1 DIT",
            "name": "DIT",
            "format": "%5.3f",
            "nameformat": "%-5s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO SEQ DIL WL0",
            "name": "wl_c",
            "format": "%5.2f",
            "nameformat": "%-5s",
        }
    )
    # keys.append({'key': 'OBJECT',                    'name':'object',       'format':'%-16s', 'nameformat':'%-16s'})
    keys.append(
        {
            "key": "HIERARCH ESO ISS AMBI FWHM START",
            "name": "seeing",
            "format": "%6.2f",
            "nameformat": "%-6s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO ISS AMBI TAU0 START",
            "name": "tau0",
            "format": "%7.5f",
            "nameformat": "%-7s",
        }
    )
    keys.append(
        {
            "key": "TELESCOP",
            "name": "telescope",
            "format": "%15s",
            "nameformat": "%-15s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO SEQ TARG FLUX L",
            "name": "Lflux",
            "format": "%8.2f",
            "nameformat": "%-8s",
        }
    )
    keys.append(
        {
            "key": "HIERARCH ESO SEQ TARG FLUX N",
            "name": "Nflux",
            "format": "%8.2f",
            "nameformat": "%-8s",
        }
    )
    if concat_tpl_starts == False:
        keys.append(
            {
                "key": "HIERARCH ESO ISS CHOP ST",
                "name": "chop_st",
                "format": "%-7s",
                "nameformat": "%-7s",
            }
        )
        keys.append(
            {
                "key": "HIERARCH ESO INS BCD1 ID",
                "name": "BCD1",
                "format": "%-4s",
                "nameformat": "%-4s",
            }
        )
        keys.append(
            {
                "key": "HIERARCH ESO INS BCD2 ID",
                "name": "BCD2",
                "format": "%-4s",
                "nameformat": "%-4s",
            }
        )

    if write_output_file:
        try:
            outf = open(outfile_path, "w")
            outf.write("night      ")
            for item in keys:
                outf.write(item["nameformat"] % (item["name"]) + " ")
            outf.write("\n")
        except FileNotFoundError as e:
            print(e)
            write_output_file = False
        except PermissionError as e:
            print(e)
            write_output_file = False

    obs_lst = []
    night = os.path.basename(os.path.dirname(rawdir))
    flst = sorted(glob.glob(rawdir + "/*.fits"))
    tpl_start_prev = ""
    cal_txt = ""

    for i in range(len(flst)):
        # print(flst[i])
        hdu = fits.open(flst[i], memmap=True, ignore_missing_end=True)
        hdr = hdu[0].header
        hdu.close()
        if "HIERARCH ESO TPL START" in hdr:
            tpl_start = hdr["HIERARCH ESO TPL START"]
        else:
            continue
        dic_val = {}
        if "HIERARCH ESO DPR TYPE" in hdr:
            dpr_type = hdr["HIERARCH ESO DPR TYPE"]
        else:
            dpr_type = ""
        if "SKY" not in dpr_type:
            if tpl_start_prev != tpl_start or concat_tpl_starts == False:
                outstr = sprintf_header(hdr, keys, dic_val)
                if (
                    "Calibrations" not in dic_val["obs_name"]
                    and "Maintenance" not in dic_val["obs_name"]
                    and "OPEN" not in dic_val["dil_name"]
                    and os.path.basename(flst[i]).startswith("MATIS.")
                ):
                    if not filter_terms:
                        select_filter = True
                    else:
                        select_filter = False
                        for item in filter_terms:
                            try:
                                hdr_val = hdr[item["key"]]
                                if hdr_val == item["value"]:
                                    select_filter = True
                            except KeyError as e:
                                pass
                    if select_filter == True:
                        if write_output_file:
                            outf.write(night + " ")
                            outf.write(outstr)
                            outf.write("\n")
                        # print(flst[i])

                        # print("{'night':'"+night+"', ", end = '')
                        # print("'tpl_start':'"+dic_val['tpl_start']+"', ", end = '')
                        if "U1234" in dic_val["telescope"]:
                            tel = "UTs"
                        elif "A1234" in dic_val["telescope"]:
                            tel = "ATs"
                        else:
                            tel = ""
                        # print("'tel':'"+tel+"',", end = '')
                        # print("'diL':'"+dic_val['dil_name']+"',", end = '')
                        # print("'diN':'"+dic_val['din_name']+"'},", end = '')
                        # print("  #"+dic_val['obs_targ_name']+", ", end = '')
                        # print(dic_val['obs_name'])
                        obs_txt = (
                            "{'night':'"
                            + night
                            + "', "
                            + "'tpl_start':'"
                            + dic_val["tpl_start"]
                            + "', "
                            + "'tel':'"
                            + tel
                            + "',"
                            + "'diL':'"
                            + dic_val["dil_name"]
                            + "',"
                            + "'diN':'"
                            + dic_val["din_name"]
                            + "'},"
                            + "  #"
                            + dic_val["obs_targ_name"]
                            + ", "
                            + dic_val["obs_name"]
                        )
                        print(obs_txt)
                        # print(hdr['HIERARCH ESO DPR CATG'])
                        if "SCIENCE" in hdr["HIERARCH ESO DPR CATG"]:
                            cal_txt += "'sci':"
                        else:
                            cal_txt += "'cal':"
                        cal_txt += obs_txt
                        cal_txt += "\n"
                        if "SCIENCE" in hdr["HIERARCH ESO DPR CATG"]:
                            cal_txt += (
                                "'sci_name':'" + dic_val["obs_targ_name"] + "'}," + "\n"
                            )
                        obs_lst.append(
                            {
                                "night": night,
                                "tpl_start": dic_val["tpl_start"],
                                "tel": tel,
                                "diL": dic_val["dil_name"],
                                "diN": dic_val["din_name"],
                            }
                        )
            tpl_start_prev = tpl_start
    if write_output_file:
        outf.close()

    return [obs_lst, cal_txt]


def sprintf_header(hdr, keys, dic):
    outstr = ""
    for item in keys:
        try:
            val = hdr[item["key"]]
        except KeyError as e:
            # print(e)
            if "s" in item["format"]:
                val = ""
            else:
                val = np.nan
        dic[item["name"]] = val
        outstr += item["format"] % (val) + " "
    return outstr
