# note: closing '/' in the paths are important! without it create_obs_lst does not work properly (os.path.basename won't return the night folder but the folder above it)
DATADIR = "/allegro6/matisse/rawdata/"  #'/allegro6/matisse/varga/rawdata/Z_CMa_img/' #'/allegro6/matisse/rawdata/'
# complete until 2025-01-06
nightdirs = [
    DATADIR + "/2024-10-05/",
    DATADIR + "/2024-10-18/",
    DATADIR + "/2024-10-19/",
    DATADIR + "/2024-10-20/",
    DATADIR + "/2024-10-26/",
    DATADIR + "/2024-10-27/",
    DATADIR + "/2024-11-09/",
    DATADIR + "/2024-11-10/",
    DATADIR + "/2024-11-26/",
    DATADIR + "/2024-11-27/",
    DATADIR + "/2024-11-28/",
    DATADIR + "/2024-11-29/",
    DATADIR + "/2024-11-30/",
    DATADIR + "/2024-12-02/",
    DATADIR + "/2024-12-08/",
    DATADIR + "/2024-12-15/",
    DATADIR + "/2024-12-16/",
    DATADIR + "/2024-12-18/",
    DATADIR + "/2024-12-29/",
    DATADIR + "/2025-01-02/",
    DATADIR + "/2025-01-06/",
    DATADIR + "/2025-01-15/",
    DATADIR + "/2025-01-16/",
    DATADIR + "/2025-01-19/",
]
