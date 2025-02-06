# plot the contents of the fits files produced by the DRS
# created: Jozsef Varga, 2019-02-27
# varga@strw.leidenuniv.nl
#
# Example usage DRS to EWS format converter (calib_raw_int_to_pk):
# outputdir = r'D:\jvarga\Data\MATISSE\commissioning\dec\2018-12-09T01_17_59\result/'
# tags = ['dil1_00','orig', 'dil0_01','dil0_10','dil0_50']
# for i in range(len(tags)):
#     inputdir = r'D:\jvarga\Data\MATISSE\commissioning\dec\2018-12-09T01_17_59\result/vis2/'+tags[i]+r'\Iter1\mat_raw_estimates.2018-12-09T01_17_59.HAWAII-2RG.rb/'
#     outputfile = outputdir + 'alfEri_drs_dec09_01_17_59_'+'vis2_'+tags[i]+'.pk'
#     calib_raw_int_to_pk(inputdir,outputfile,1,verbose=False)
#
# NOT YET FINISHED
########################################################################

import pickle
from shutil import copyfile

import matplotlib
import numpy as np
from astropy.io import fits
from mcdb import wutil as wu

from ...plotting.show_allred import open_fits

# plot style configuration
matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
matplotlib.rcParams.update({"font.size": 10})


# outhdul['OI_VIS2'].data['VIS2DATA'][j] = np.abs(dic['opddata']['vis'][j])**2
# outhdul['OI_VIS'].data['VISAMP'][j]    = np.abs(dic['opddata']['flux'][j])
# outhdul['OI_VIS'].data['VISPHI'][j]    = np.angle(dic['opddata']['flux'][j]) * 180.0 / np.pi
# outhdul['OI_T3'].data['T3PHI'][j]      = dic['opddata']['closurephase'][j]
# outhdul['OI_FLUX'].data['FLUXDATA'][j] = dic['photdata']['photave'][j+6]
# or outhdul['OI_FLUX'].data['FLUXDATA'][j] = 0.0*outhdul['OI_FLUX'].data['FLUX'][j]


# pkfile
# if string: path to pkfile path
# if list: reduced data object (list of dictionaries)
# version: 'alpha' wl array is dic['wave'], 'omicron': wl array is dic['opddata']['wave']
def pk_2_oifits(
    pkfile, outdir, oifits_template="template_oifits.fits", version="omicron"
):
    # read in ews pickle file
    if isinstance(pkfile, str):
        try:
            data = wu.mrestore(pkfile)
        except FileNotFoundError as e:
            print(e)
            return
        N_data = len(data)
    else:
        data = pkfile
        N_data = len(data)

    for i in range(N_data):
        dic = data[i]
        hdr = dic["header"]
        dpr_catg = hdr["HIERARCH ESO DPR CATG"]
        bcd1 = dic["bcd1name"]
        bcd2 = dic["bcd2name"]
        hdr["HIERARCH ESO INS BCD1 ID"] = bcd1
        hdr["HIERARCH ESO INS BCD2 ID"] = bcd2
        chop_status = hdr["HIERARCH ESO ISS CHOP ST"]
        band = dic["band"]
        # print('wl old new ',len(dic['wave']),len(dic['opddata']['wave']))
        if version == "alpha":
            wav = dic["wave"]
        if version == "omicron":
            wav = dic["opddata"]["wave"]
        if "N" in band:
            if np.min(wav) < 6.5:
                incl_idx = wav > 6.5
                wav = wav[incl_idx]
                if dic["opddata"]["vis"] is not None:
                    dic["opddata"]["vis"] = dic["opddata"]["vis"][:, incl_idx]
                dic["opddata"]["flux"] = dic["opddata"]["flux"][:, incl_idx]
                dic["opddata"]["closurephase"] = dic["opddata"]["closurephase"][
                    :, incl_idx
                ]
                if "photdata" in dic:
                    dic["photdata"]["photave"] = dic["photdata"]["photave"][:, incl_idx]

        # print(bcd1,bcd2)
        if bcd1 == "OUT":
            if bcd2 == "OUT":
                # OO
                if "L" in band:
                    if chop_status == "F":
                        id = 1
                    if chop_status == "T":
                        id = 5
                if "N" in band:
                    id = 1
            else:
                # OI
                id = 4
        else:
            if bcd2 == "OUT":
                # IO
                id = 3
            else:
                # II
                if "L" in band:
                    if chop_status == "F":
                        id = 2
                    if chop_status == "T":
                        id = 6
                if "N" in band:
                    id = 2
        if bcd1 == "":
            id = "AVG"
            if "L" in band:
                if chop_status == "T":
                    id = "CHOPPED_AVG"
        if isinstance(id, str):
            if dpr_catg == "SCIENCE":
                oifits_file = "TARGET_RAW_INT_%s.fits" % id
            else:
                oifits_file = "CALIB_RAW_INT_%s.fits" % id
        else:
            if dpr_catg == "SCIENCE":
                oifits_file = "TARGET_RAW_INT_%04d.fits" % id
            else:
                oifits_file = "CALIB_RAW_INT_%04d.fits" % id

        copyfile(oifits_template, outdir + "/" + oifits_file)
        outhdul = fits.open(outdir + "/" + oifits_file, mode="update")
        outhdul[0].header = hdr
        if dic["header"]["HIERARCH ESO DPR CATG"] == "CALIB":
            outhdul[0].header["HIERARCH ESO PRO SCIENCE"] = False
            outhdul[0].header["HIERARCH ESO PRO CATG"] = "CALIB_RAW_INT"
        else:
            outhdul[0].header["HIERARCH ESO PRO SCIENCE"] = True
            outhdul[0].header["HIERARCH ESO PRO CATG"] = "TARGET_RAW_INT"
        outhdul[0].header["HIERARCH ESO CFG BCD MODE"] = (
            dic["bcd1name"] + "-" + dic["bcd2name"]
        )

        # outhdul[0].header['HIERARCH ESO INS FIL NAME'] =
        # outhdul[0].header['HIERARCH ESO INS FIL ID'] =
        # outhdul[0].header['HIERARCH ESO INS FIL NO'] =

        # for i in hdr.keys():
        #     print(i)

        # OI_TARGET
        target_id = 1
        c1 = fits.Column(name="TARGET_ID", format="1I", array=[target_id])
        c2 = fits.Column(
            name="TARGET", format="%dA" % len(dic["targ"]), array=[dic["targ"]]
        )
        c3 = fits.Column(name="RAEP0", format="1D", unit="deg", array=[dic["ra"]])
        c4 = fits.Column(name="DECEP0", format="1D", unit="deg", array=[dic["dec"]])
        c5 = fits.Column(name="EQUINOX", format="1E", unit="yr", array=[0.0])
        c6 = fits.Column(name="RA_ERR", format="1D", unit="deg", array=[0.0])
        c7 = fits.Column(name="DEC_ERR", format="1D", unit="deg", array=[0.0])
        c8 = fits.Column(name="SYSVEL", format="1D", unit="m/s", array=[0.0])
        c9 = fits.Column(name="VELTYP", format="8A", array=["HELIOCEN"])
        c10 = fits.Column(name="VELDEF", format="8A", array=["OPTICAL"])
        c11 = fits.Column(name="PMRA", format="1D", unit="deg/yr", array=[0.0])
        c12 = fits.Column(name="PMDEC", format="1D", unit="deg/yr", array=[0.0])
        c13 = fits.Column(name="PMRA_ERR", format="1D", unit="deg/yr", array=[0.0])
        c14 = fits.Column(name="PMDEC_ERR", format="1D", unit="deg/yr", array=[0.0])
        c15 = fits.Column(name="PARALLAX", format="1E", unit="deg", array=[0.0])
        c16 = fits.Column(name="PARA_ERR", format="1E", unit="deg", array=[0.0])
        c17 = fits.Column(name="SPECTYP", format="7A", array=["UNKNOWN"])
        if dic["header"]["HIERARCH ESO DPR CATG"] == "CALIB":
            c18 = fits.Column(name="CATEGORY", format="3A", array=["CAL"])
        else:
            c18 = fits.Column(name="CATEGORY", format="3A", array=["SCI"])
        del outhdul["OI_TARGET"]
        oi_target_hdu = fits.BinTableHDU.from_columns(
            [
                c1,
                c2,
                c3,
                c4,
                c5,
                c6,
                c7,
                c8,
                c9,
                c10,
                c11,
                c12,
                c13,
                c14,
                c15,
                c16,
                c17,
                c18,
            ],
            name="OI_TARGET",
        )
        outhdul.append(oi_target_hdu)
        outhdul["OI_TARGET"].header["EXTNAME"] = "OI_TARGET"
        outhdul["OI_TARGET"].header["OI_REVN"] = 1

        # OI_ARRAY
        for j in range(4):
            outhdul["OI_ARRAY"].data["TEL_NAME"][j] = hdr[
                "HIERARCH ESO ISS CONF T%dNAME" % (j + 1)
            ]
            outhdul["OI_ARRAY"].data["STA_NAME"][j] = hdr[
                "HIERARCH ESO ISS CONF STATION%d" % (j + 1)
            ]
            outhdul["OI_ARRAY"].data["STA_INDEX"][j] = j + 1
            if outhdul["OI_ARRAY"].data["TEL_NAME"][j].startswith("U"):
                outhdul["OI_ARRAY"].data["DIAMETER"][j] = 8.2
                outhdul["OI_ARRAY"].data["FOV"][j] = 1.0  # arcsec, check!!!
            else:
                outhdul["OI_ARRAY"].data["DIAMETER"][j] = 1.8
                outhdul["OI_ARRAY"].data["FOV"][j] = 1.0  # arcsec

            outhdul["OI_ARRAY"].data["STAXYZ"][j] = np.array(
                [
                    hdr["HIERARCH ESO ISS CONF T%dX" % (j + 1)],
                    hdr["HIERARCH ESO ISS CONF T%dY" % (j + 1)],
                    hdr["HIERARCH ESO ISS CONF T%dZ" % (j + 1)],
                ]
            )
            outhdul["OI_ARRAY"].data["FOVTYPE"][j] = "RADIUS"

        # OI_WAVELENGTH
        # print(len(dic['wave']))

        outhdul["OI_WAVELENGTH"].data = np.resize(
            outhdul["OI_WAVELENGTH"].data, (len(wav),)
        )
        outhdul["OI_WAVELENGTH"].data["EFF_WAVE"] = wav * 1e-6  # m
        outhdul["OI_WAVELENGTH"].data["EFF_BAND"] = wav * 0.0 + 7.521082e-08
        Nwl = len(wav)

        sta_hdr = np.array(
            [
                hdr["HIERARCH ESO ISS CONF STATION1"],
                hdr["HIERARCH ESO ISS CONF STATION2"],
                hdr["HIERARCH ESO ISS CONF STATION3"],
                hdr["HIERARCH ESO ISS CONF STATION4"],
            ]
        )
        sta_nr = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        pbl_hdr = np.array(
            [
                hdr["HIERARCH ESO ISS PBL12 START"],
                hdr["HIERARCH ESO ISS PBL13 START"],
                hdr["HIERARCH ESO ISS PBL14 START"],
                hdr["HIERARCH ESO ISS PBL23 START"],
                hdr["HIERARCH ESO ISS PBL24 START"],
                hdr["HIERARCH ESO ISS PBL34 START"],
            ]
        )
        pbla_hdr = np.array(
            [
                hdr["HIERARCH ESO ISS PBLA12 START"],
                hdr["HIERARCH ESO ISS PBLA13 START"],
                hdr["HIERARCH ESO ISS PBLA14 START"],
                hdr["HIERARCH ESO ISS PBLA23 START"],
                hdr["HIERARCH ESO ISS PBLA24 START"],
                hdr["HIERARCH ESO ISS PBLA34 START"],
            ]
        )
        u_hdr = np.array(pbl_hdr * np.sin(pbla_hdr * np.pi / 180.0))
        v_hdr = np.array(pbl_hdr * np.cos(pbla_hdr * np.pi / 180.0))

        # OI_VIS2
        NV2 = len(dic["pbl"])
        del outhdul["OI_VIS2"]
        c1 = fits.Column(name="TARGET_ID", format="I")
        c2 = fits.Column(name="TIME", format="D", unit="s")
        c3 = fits.Column(name="MJD", format="D", unit="day")
        c4 = fits.Column(name="INT_TIME", format="D", unit="s")
        c5 = fits.Column(name="VIS2DATA", format="%dD" % Nwl)
        c6 = fits.Column(name="VIS2ERR", format="%dD" % Nwl)
        c7 = fits.Column(name="UCOORD", format="D", unit="m")
        c8 = fits.Column(name="VCOORD", format="D", unit="m")
        c9 = fits.Column(name="STA_INDEX", format="2I")
        c10 = fits.Column(name="FLAG", format="%dL" % Nwl)
        # coldefs=fits.ColDefs(outhdul['OI_VIS2'])
        # coldefs.change_attrib('VIS2DATA','format','%dD'%Nwl)
        # coldefs.change_attrib('VIS2ERR','format','%dD'%Nwl)
        # coldefs.change_attrib('FLAG','format','%dL'%Nwl)
        # oi_vis2_hdu = fits.BinTableHDU.from_columns(coldefs)
        oi_vis2_hdu = fits.BinTableHDU.from_columns(
            [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10],
            name="OI_VIS2",
            nrows=NV2,
            fill=True,
        )
        outhdul.append(oi_vis2_hdu)

        outhdul["OI_VIS2"].header["DATE-OBS"] = outhdul[0].header["DATE-OBS"]
        outhdul["OI_VIS2"].header["EXTNAME"] = "OI_VIS2"
        outhdul["OI_VIS2"].header["EXTVER"] = 1
        outhdul["OI_VIS2"].header["OI_REVN"] = 1
        outhdul["OI_VIS2"].header["ARRNAME"] = "VLTI"
        outhdul["OI_VIS2"].header["INSNAME"] = "MATISSE"
        outhdul["OI_VIS2"].data["TARGET_ID"] = np.array([target_id] * NV2)
        outhdul["OI_VIS2"].data["TIME"] = np.array([0.0] * NV2)
        outhdul["OI_VIS2"].data["MJD"] = np.array([dic["mjd-obs"]] * NV2)
        outhdul["OI_VIS2"].data["INT_TIME"] = np.array([hdr["EXPTIME"]] * NV2)
        outhdul["OI_VIS2"].data["UCOORD"] = np.array(
            [dic["uvcoord"][i][0] for i in range(NV2)]
        )
        outhdul["OI_VIS2"].data["VCOORD"] = np.array(
            [dic["uvcoord"][i][1] for i in range(NV2)]
        )
        # print(dic['pbl'])
        # print(dic['pbla'])
        # print(dic['pbl'] * np.sin(dic['pbla'] * np.pi / 180.0),dic['pbl'] * np.cos(dic['pbla'] * np.pi / 180.0))
        for j in range(NV2):
            # print(len(dic['opddata']['vis'][j]))
            # print(len(dic['opddata']['flux'][j]))
            try:
                if dic["opddata"]["vis"] is not None:
                    outhdul["OI_VIS2"].data["VIS2DATA"][j] = (
                        np.abs(dic["opddata"]["vis"][j]) ** 2
                    )
                    outhdul["OI_VIS2"].data["FLAG"][j] = np.array([False] * Nwl)
                else:
                    # print((outhdul['OI_VIS2'].columns))
                    # print(np.shape(outhdul['OI_VIS2'].data['VIS2DATA']))
                    outhdul["OI_VIS2"].data["VIS2DATA"][j] = (
                        np.abs(dic["opddata"]["flux"][j])
                    ) ** 2  # squared corr. flux
                    outhdul["OI_VIS2"].data["FLAG"][j] = np.array([True] * Nwl)
            except KeyError as e:
                print(e)
                outhdul["OI_VIS2"].data["VIS2DATA"][j] = (
                    0.0 * outhdul["OI_VIS2"].data["VIS2DATA"][j]
                )
                outhdul["OI_VIS2"].data["FLAG"][j] = np.array([True] * Nwl)
            outhdul["OI_VIS2"].data["VIS2ERR"][j] = (
                0.0 * outhdul["OI_VIS2"].data["VIS2ERR"][j]
            )
            duv = (outhdul["OI_VIS2"].data["UCOORD"][j] - u_hdr) ** 2 + (
                outhdul["OI_VIS2"].data["VCOORD"][j] - v_hdr
            ) ** 2
            idx = sta_nr[np.argmin(duv)]
            oi_arr_idx1 = np.where(
                sta_hdr[idx[0]] == outhdul["OI_ARRAY"].data["STA_NAME"]
            )[0][0]
            oi_arr_idx2 = np.where(
                sta_hdr[idx[1]] == outhdul["OI_ARRAY"].data["STA_NAME"]
            )[0][0]
            outhdul["OI_VIS2"].data["STA_INDEX"][j] = np.array(
                [
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx1],
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx2],
                ]
            )

        # OI_T3
        NT3 = len(dic["opddata"]["closurephase"])
        del outhdul["OI_T3"]
        c1 = fits.Column(name="TARGET_ID", format="1I")
        c2 = fits.Column(name="TIME", format="1D", unit="s")
        c3 = fits.Column(name="MJD", format="1D", unit="day")
        c4 = fits.Column(name="INT_TIME", format="1D", unit="s")
        c5 = fits.Column(name="T3AMP", format="%dD" % Nwl)
        c6 = fits.Column(name="T3AMPERR", format="%dD" % Nwl)
        c7 = fits.Column(name="T3PHI", format="%dD" % Nwl, unit="deg")
        c8 = fits.Column(name="T3PHIERR", format="%dD" % Nwl, unit="deg")
        c9 = fits.Column(name="U1COORD", format="1D", unit="m")
        c10 = fits.Column(name="V1COORD", format="1D", unit="m")
        c11 = fits.Column(name="U2COORD", format="1D", unit="m")
        c12 = fits.Column(name="V2COORD", format="1D", unit="m")
        c13 = fits.Column(name="STA_INDEX", format="3I")
        c14 = fits.Column(name="FLAG", format="%dL" % Nwl)
        oi_t3_hdu = fits.BinTableHDU.from_columns(
            [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14],
            name="OI_T3",
            nrows=NT3,
            fill=True,
        )
        outhdul.append(oi_t3_hdu)

        outhdul["OI_T3"].header["DATE-OBS"] = outhdul[0].header["DATE-OBS"]
        outhdul["OI_T3"].header["EXTNAME"] = "OI_T3"
        outhdul["OI_T3"].header["EXTVER"] = 1
        outhdul["OI_T3"].header["OI_REVN"] = 1
        outhdul["OI_T3"].header["ARRNAME"] = "VLTI"
        outhdul["OI_T3"].header["INSNAME"] = "MATISSE"
        outhdul["OI_T3"].data["TARGET_ID"] = np.array([target_id] * NT3)
        outhdul["OI_T3"].data["TIME"] = np.array([0.0] * NT3)
        outhdul["OI_T3"].data["MJD"] = np.array([dic["mjd-obs"]] * NT3)
        outhdul["OI_T3"].data["INT_TIME"] = np.array([hdr["EXPTIME"]] * NT3)
        # print(u_hdr)
        # print(v_hdr)
        # print(sta_hdr)
        # print(outhdul['OI_ARRAY'].data['STA_NAME'])
        # print(outhdul['OI_ARRAY'].data['STA_INDEX'])
        for j in range(NT3):
            outhdul["OI_T3"].data["T3AMP"][j] = 0.0 * dic["opddata"]["closurephase"][j]
            outhdul["OI_T3"].data["T3AMPERR"][j] = (
                0.0 * dic["opddata"]["closurephase"][j]
            )
            outhdul["OI_T3"].data["T3PHI"][j] = dic["opddata"]["closurephase"][j]
            outhdul["OI_T3"].data["T3PHIERR"][j] = (
                0.0 * dic["opddata"]["closurephase"][j]
            )
            outhdul["OI_T3"].data["U1COORD"][j] = dic["uvcoord1"][j][0]
            outhdul["OI_T3"].data["V1COORD"][j] = dic["uvcoord1"][j][1]
            outhdul["OI_T3"].data["U2COORD"][j] = dic["uvcoord2"][j][0]
            outhdul["OI_T3"].data["V2COORD"][j] = dic["uvcoord2"][j][1]
            duv = (outhdul["OI_T3"].data["U1COORD"][j] - u_hdr) ** 2 + (
                outhdul["OI_T3"].data["V1COORD"][j] - v_hdr
            ) ** 2
            idx = sta_nr[np.argmin(duv)]
            min_duv = np.sqrt(duv[np.argmin(duv)])
            match = False
            # print(duv[idx])
            if min_duv < 2.0:
                oi_arr_idx1 = np.where(
                    sta_hdr[idx[0]] == outhdul["OI_ARRAY"].data["STA_NAME"]
                )[0][0]
                oi_arr_idx2 = np.where(
                    sta_hdr[idx[1]] == outhdul["OI_ARRAY"].data["STA_NAME"]
                )[0][0]
                match = True
            else:
                oi_arr_idx1 = 0
                oi_arr_idx2 = 0
                match = False

            duv = (outhdul["OI_T3"].data["U2COORD"][j] - u_hdr) ** 2 + (
                outhdul["OI_T3"].data["V2COORD"][j] - v_hdr
            ) ** 2
            idx = sta_nr[np.argmin(duv)]
            min_duv = np.sqrt(duv[np.argmin(duv)])
            if min_duv < 2.0:
                oi_arr_idx3 = np.where(
                    sta_hdr[idx[0]] == outhdul["OI_ARRAY"].data["STA_NAME"]
                )[0][0]
                oi_arr_idx4 = np.where(
                    sta_hdr[idx[1]] == outhdul["OI_ARRAY"].data["STA_NAME"]
                )[0][0]
                match = match & True
            else:
                oi_arr_idx3 = 0
                oi_arr_idx4 = 0
                match = False

            arr = np.array(
                [
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx1],
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx2],
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx3],
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx4],
                ]
            )
            if match == True:
                idx = np.unique(arr, return_index=True)[1]
            else:
                print("Unable to identify closure phase baselines.")
                idx = np.array([0, 0, 0])

            # print(j)
            # print(dic['uvcoord1'][j],dic['uvcoord2'][j])
            # print([arr[idx] for idx in sorted(idx)])
            # print(np.array([arr[idx] for idx in sorted(idx)]))
            outhdul["OI_T3"].data["STA_INDEX"][j] = np.array(
                [arr[idx] for idx in sorted(idx)]
            )
            outhdul["OI_T3"].data["FLAG"][j] = np.array([False] * Nwl)

        # OI_VIS
        NV = len(dic["pbl"])
        del outhdul["OI_VIS"]
        c1 = fits.Column(name="TARGET_ID", format="1I")
        c2 = fits.Column(name="TIME", format="1D", unit="s")
        c3 = fits.Column(name="MJD", format="1D", unit="day")
        c4 = fits.Column(name="INT_TIME", format="1D", unit="s")
        c5 = fits.Column(name="VISAMP", format="%dD" % Nwl)
        c6 = fits.Column(name="VISAMPERR", format="%dD" % Nwl)
        c7 = fits.Column(name="VISPHI", format="%dD" % Nwl, unit="deg")
        c8 = fits.Column(name="VISPHIERR", format="%dD" % Nwl, unit="deg")
        c9 = fits.Column(name="UCOORD", format="1D", unit="m")
        c10 = fits.Column(name="VCOORD", format="1D", unit="m")
        c11 = fits.Column(name="STA_INDEX", format="2I")
        c12 = fits.Column(name="FLAG", format="%dL" % Nwl)
        oi_vis_hdu = fits.BinTableHDU.from_columns(
            [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12],
            name="OI_VIS",
            nrows=NV,
            fill=True,
        )
        outhdul.append(oi_vis_hdu)

        outhdul["OI_VIS"].header["EXTNAME"] = "OI_VIS"
        outhdul["OI_VIS"].header["EXTVER"] = 1
        outhdul["OI_VIS"].header["OI_REVN"] = 1
        outhdul["OI_VIS"].header["ARRNAME"] = "VLTI"
        outhdul["OI_VIS"].header["INSNAME"] = "MATISSE"
        outhdul["OI_VIS"].header["DATE-OBS"] = outhdul[0].header["DATE-OBS"]
        outhdul["OI_VIS"].header["AMPTYP"] = "correlated flux"
        outhdul["OI_VIS"].header["PHITYP"] = "differential"
        outhdul["OI_VIS"].data["TARGET_ID"] = np.array([target_id] * NV)
        outhdul["OI_VIS"].data["TIME"] = np.array([0.0] * NV)
        outhdul["OI_VIS"].data["MJD"] = np.array([dic["mjd-obs"]] * NV)
        outhdul["OI_VIS"].data["INT_TIME"] = np.array([hdr["EXPTIME"]] * NV)
        outhdul["OI_VIS"].data["UCOORD"] = np.array(
            [dic["uvcoord"][i][0] for i in range(NV)]
        )
        outhdul["OI_VIS"].data["VCOORD"] = np.array(
            [dic["uvcoord"][i][1] for i in range(NV)]
        )
        for j in range(NV):
            outhdul["OI_VIS"].data["VISAMP"][j] = np.abs(dic["opddata"]["flux"][j])
            outhdul["OI_VIS"].data["VISAMPERR"][j] = (
                0.0 * outhdul["OI_VIS"].data["VISAMPERR"][j]
            )
            outhdul["OI_VIS"].data["VISPHI"][j] = (
                np.angle(dic["opddata"]["flux"][j]) * 180.0 / np.pi
            )
            outhdul["OI_VIS"].data["VISPHIERR"][j] = (
                0.0 * outhdul["OI_VIS"].data["VISPHIERR"][j]
            )
            duv = (outhdul["OI_VIS"].data["UCOORD"][j] - u_hdr) ** 2 + (
                outhdul["OI_VIS"].data["VCOORD"][j] - v_hdr
            ) ** 2
            idx = sta_nr[np.argmin(duv)]
            oi_arr_idx1 = np.where(
                sta_hdr[idx[0]] == outhdul["OI_ARRAY"].data["STA_NAME"]
            )[0][0]
            oi_arr_idx2 = np.where(
                sta_hdr[idx[1]] == outhdul["OI_ARRAY"].data["STA_NAME"]
            )[0][0]
            outhdul["OI_VIS"].data["STA_INDEX"][j] = np.array(
                [
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx1],
                    outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx2],
                ]
            )
            outhdul["OI_VIS"].data["FLAG"][j] = np.array([False] * Nwl)

        # OI_FLUX
        NF = 4  #!!!!!!!
        del outhdul["OI_FLUX"]
        if "photdata" in dic:  # check
            c1 = fits.Column(name="TARGET_ID", format="1I")
            c2 = fits.Column(name="TIME", format="1D", unit="s")
            c3 = fits.Column(name="MJD", format="1D", unit="day")
            c4 = fits.Column(name="INT_TIME", format="1D", unit="s")
            c5 = fits.Column(name="FLUXDATA", format="%dD" % Nwl)
            c6 = fits.Column(name="FLUXERR", format="%dD" % Nwl)
            c7 = fits.Column(name="STA_INDEX", format="1I")
            c8 = fits.Column(name="FLAG", format="%dL" % Nwl)
            oi_flux_hdu = fits.BinTableHDU.from_columns(
                [c1, c2, c3, c4, c5, c6, c7, c8], name="OI_FLUX", nrows=NF, fill=True
            )
            outhdul.append(oi_flux_hdu)

            outhdul["OI_FLUX"].header["EXTNAME"] = "OI_FLUX"
            outhdul["OI_FLUX"].header["EXTVER"] = 1
            outhdul["OI_FLUX"].header["OI_REVN"] = 1
            outhdul["OI_FLUX"].header["DATE-OBS"] = outhdul[0].header["DATE-OBS"]
            outhdul["OI_FLUX"].header["ARRNAME"] = "VLTI"
            outhdul["OI_FLUX"].header["INSNAME"] = "MATISSE"
            outhdul["OI_FLUX"].header["FOV"] = 1.0
            outhdul["OI_FLUX"].header["FOVTYPE"] = "RADIUS"
            outhdul["OI_FLUX"].header["CALSTAT"] = "U"
            outhdul["OI_FLUX"].data["TARGET_ID"] = np.array([target_id] * NF)
            outhdul["OI_FLUX"].data["TIME"] = np.array([0.0] * NF)
            outhdul["OI_FLUX"].data["MJD"] = np.array([dic["mjd-obs"]] * NF)
            outhdul["OI_FLUX"].data["INT_TIME"] = np.array([hdr["EXPTIME"]] * NF)
            for j in range(NF):
                try:
                    # print(dic.keys())
                    # print(['data'])
                    # print(dic['photdata']['data'])
                    # print('meandata')
                    # print(dic['photdata']['meandata'])
                    # print('photdata')
                    # print(dic['photdata']['photdata'])
                    # print('photdata photdata',dic['photdata']['photdata'][0].shape)
                    # print('photdata photflag',len(dic['photdata']['photflag']))
                    # print(dic['photdata']['photflag'])
                    # print('opddata vis', dic['opddata']['vis'][0].shape)
                    # print('opddata flux',dic['opddata']['flux'][0].shape)
                    # print('photave',dic['photdata']['photave'][0].shape)
                    # print('data',dic['photdata']['data'][0].shape)
                    # print('meandata',dic['photdata']['meandata'][0].shape)
                    # print('photave')
                    # print(dic['photdata']['photave'])
                    # print('vis')
                    # print(dic['opddata']['vis'])
                    # print(dic['photdata']['photflag'])
                    # print(dic['photdata']['photflag4'])
                    outhdul["OI_FLUX"].data["FLUXDATA"][j] = dic["photdata"]["photave"][
                        j + 6
                    ]
                except KeyError as e:
                    print(e)
                    outhdul["OI_FLUX"].data["FLUXDATA"][j] = (
                        0.0 * outhdul["OI_FLUX"].data["FLUX"][j]
                    )
                outhdul["OI_FLUX"].data["FLUXERR"][j] = (
                    0.0 * outhdul["OI_FLUX"].data["FLUXERR"][j]
                )
                oi_arr_idx = np.where(
                    sta_hdr[j] == outhdul["OI_ARRAY"].data["STA_NAME"]
                )[0][0]
                outhdul["OI_FLUX"].data["STA_INDEX"][j] = outhdul["OI_ARRAY"].data[
                    "STA_INDEX"
                ][oi_arr_idx]
                outhdul["OI_FLUX"].data["FLAG"][j] = np.array([False] * Nwl)

        # TF2
        if dic["header"]["HIERARCH ESO DPR CATG"] == "CALIB":
            NV2 = len(dic["pbl"])
            c1 = fits.Column(name="TIME", format="1D", unit="s")
            c2 = fits.Column(name="MJD", format="1D", unit="day")
            c3 = fits.Column(name="INT_TIME", format="1D", unit="s")
            c4 = fits.Column(name="TF2", format="%dD" % Nwl)
            c5 = fits.Column(name="TF2ERR", format="%dD" % Nwl)
            c6 = fits.Column(name="TF", format="%dD" % Nwl)
            c7 = fits.Column(name="TFERR", format="%dD" % Nwl)
            c8 = fits.Column(name="STA_INDEX", format="2I")
            oi_tf2_hdu = fits.BinTableHDU.from_columns(
                [c1, c2, c3, c4, c5, c6, c7, c8], name="TF2", nrows=NV2, fill=True
            )
            outhdul.append(oi_tf2_hdu)

            outhdul["TF2"].header["DATE-OBS"] = outhdul[0].header["DATE-OBS"]
            outhdul["TF2"].header["EXTNAME"] = "TF2"
            outhdul["TF2"].header["EXTVER"] = 1
            outhdul["TF2"].header["OI_REVN"] = 1
            outhdul["TF2"].header["ARRNAME"] = "VLTI"
            outhdul["TF2"].header["INSNAME"] = "MATISSE"
            outhdul["TF2"].data["TIME"] = np.array([0.0] * NV2)
            outhdul["TF2"].data["MJD"] = np.array([dic["mjd-obs"]] * NV2)
            outhdul["TF2"].data["INT_TIME"] = np.array([hdr["EXPTIME"]] * NV2)
            # print(dic['pbl'])
            # print(dic['pbla'])
            # print(dic['pbl'] * np.sin(dic['pbla'] * np.pi / 180.0),dic['pbl'] * np.cos(dic['pbla'] * np.pi / 180.0))
            for j in range(NV2):
                #!!!!!!!!!!!!!!!! apply cal diameter !!!!!!!!!!!!!!!!!!!!!
                try:
                    if dic["opddata"]["vis"] is not None:
                        outhdul["TF2"].data["TF2"][j] = (
                            np.abs(dic["opddata"]["vis"][j]) ** 2
                        )
                        outhdul["TF2"].data["TF"][j] = np.abs(dic["opddata"]["vis"][j])
                    else:
                        outhdul["TF2"].data["TF2"][j] = (
                            np.abs(dic["opddata"]["flux"][j])
                        ) ** 2  # squared corr. flux
                        outhdul["TF2"].data["TF"][j] = np.abs(dic["opddata"]["flux"][j])
                except KeyError as e:
                    print(e)
                outhdul["TF2"].data["TF2ERR"][j] = (
                    0.0 * outhdul["TF2"].data["TF2ERR"][j]
                )
                outhdul["TF2"].data["TFERR"][j] = 0.0 * outhdul["TF2"].data["TFERR"][j]
                duv = (outhdul["OI_VIS2"].data["UCOORD"][j] - u_hdr) ** 2 + (
                    outhdul["OI_VIS2"].data["VCOORD"][j] - v_hdr
                ) ** 2
                idx = sta_nr[np.argmin(duv)]
                oi_arr_idx1 = np.where(
                    sta_hdr[idx[0]] == outhdul["OI_ARRAY"].data["STA_NAME"]
                )[0][0]
                oi_arr_idx2 = np.where(
                    sta_hdr[idx[1]] == outhdul["OI_ARRAY"].data["STA_NAME"]
                )[0][0]
                outhdul["TF2"].data["STA_INDEX"][j] = np.array(
                    [
                        outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx1],
                        outhdul["OI_ARRAY"].data["STA_INDEX"][oi_arr_idx2],
                    ]
                )
        outhdul.flush()  # changes are written back to original.fits
        outhdul.close()


# dict_keys(['phot', 'wave', 'cflux', 'flux', 'vis', 'mjd', 'meanspectrum', 'photimage', 'opd', 'pbl', 'pbla', 'header',
# 'bcd', 'flags', 'baseflag', 'fflag', 'tau', 'ra', 'dec', 'mjd-obs', 'targ', 'file', 'sky', 'closurePhase', 'uvcoord1', 'uvcoord2'])


# Convert the final results of the DRS to a pickle file
# The output pickle file has the same format as the result of the EWS reduction
# The conversion from the DRS quantities to the EWS quantities is problematic.
# Example: the EWS pickle file stores the complex visibility. However, the DRS
# outputs (CALIB_RAW_INT) contain the squared visibility (OI_VIS2), thus the phase information is lost.
# The OI_VIS table contains VISAMP and VISPHI, but it seems that these quantities are not
# the modulus and phase part of the complex visibility.
# Arguments:
# inputdir: directory containing the DRS data products
# outputfile: path of the output pickle file
# id: select which CALIB_RAW_INT file you want to convert (the id e.g., can correspond to a specific BCD combination,
# typical values: 1,2,3,4)
# verbose (True or False): verbosity of output messages
def calib_raw_int_to_pk(inputdir, outputfile, id, verbose=False):
    dic = open_fits(
        inputdir + "/CALIB_RAW_INT_" + "%04d" % (id) + ".fits", verbose=verbose
    )
    try:
        phot = dic["FLUX"]["FLUX"]
    except KeyError as e:
        print("Total flux table not found. ")
        print(e)
        phot = []
    wl = dic["WLEN"] * 1e6  # (um)
    flux = dic["VIS"]["CFXAMP"]  #!!! the OIFITS files do not contain CFXAMP
    modulus = dic["VIS"]["VISAMP"]
    phase = dic["VIS"]["VISPHI"] * np.pi / 180.0
    vis = modulus * np.cos(phase) + modulus * np.sin(phase) * np.complex(
        0, 1
    )  # np.sqrt(dic['VIS2']['VIS2']) #or ???? dic['VIS']['VISPHI']
    h = dic["HDR"]
    tau = dic["TAU0"]
    ra = dic["RA"]
    dec = dic["DEC"]
    mjdobs = dic["MJD-OBS"]
    targ = dic["TARGET"]
    file = dic["HDR"]["ARCFILE"]
    skyfile = ""
    u = dic["VIS2"]["U"]
    v = dic["VIS2"]["V"]
    pbl = np.sqrt(u**2 + v**2)

    dic = open_fits(inputdir + "OI_OPDWVPO_" + "%04d" % (id) + ".fits", verbose=verbose)
    try:
        mjd = dic["OPD"]["MJD"]
        opd = dic["OPD"]["OPD"]
    except KeyError:
        mjd = []
        opd = []
    cflux = {}
    cflux["mjd"] = mjd
    cflux["wave"] = wl
    cflux["delay"] = []
    cflux["fdata"] = []
    cflux["cdata"] = []
    cflux["pdata"] = []
    cflux["localopd"] = []
    cflux["opd"] = opd
    cflux["opdp"] = []
    cflux["phase"] = []
    cflux["opdf"] = []
    cflux["flux"] = flux
    cflux["opdimage"] = []

    out = []
    out.append(
        {
            "phot": phot,
            "wave": wl,
            "cflux": cflux,
            "flux": flux,
            "vis": vis,
            "mjd": mjd,
            "opd": cflux["opd"],
            "pbl": pbl,
            "header": h,
            "tau": tau,
            "ra": ra,
            "dec": dec,
            "mjd-obs": mjdobs,
            "targ": targ,
            "file": file,
            "sky": skyfile,
        }
    )
    rdata = out
    # print(rdata)

    sf = open(outputfile, "wb")
    p = pickle.Pickler(sf)
    p.dump(rdata)
    sf.close()
    del p


# RESDIR_L = '/Users/jvarga/Data/MATISSE/matisse_red/ews/2019-03-22/2019-03-23T08_41_19/Iter1/raw_estimates.2019-03-23T08_41_19.HAWAII-2RG.rb'
# TPL_START = '2019-03-23T08:41:19'
# pk_2_oifits(RESDIR_L+'/'+TPL_START.replace(':', '_')+'.tpl.pk', RESDIR_L,oifits_template='/Users/jvarga/Dokumentumok/MATISSE/pro/template_oifits.fits')
# inputdir = RESDIR_L
# show_allred(inputdir, outputdir=inputdir+'/plots/', verbose=False)
