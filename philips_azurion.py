#Creation of the dataset for Philips Azurion
#.npz to .npz files - image normalization
#we deleted a frame around the image and resized it

import cv2
import numpy as np


def normalize_img(img):
    """
    normalize function to [0, 255]

    img - the initial range is [0, 4095]
    """

    img = img / 4095 * 255
    return img


#series pathes with a picularity
series_paths = [
    "***.npz",
    "***.npz",
    "***.npz",
    "***.npz",
]

for j, series_path in enumerate(series_paths):

    #load angiogrpahic videos from .npz
    hearth_npz_info = np.load(series_path)
    pixel_array = hearth_npz_info["pixel_array"]
    pixel_array = np.copy(normalize_img(pixel_array).astype(np.uint8))

    #new video
    new_pixel_array = np.zeros((len(pixel_array), 512, 512))

    #delete frames
    for i in range(len(pixel_array)):
        if j == 0:
            сur_img = np.copy(pixel_array[i][115:390, 105:410])
        elif j == 1:
            сur_img = np.copy(pixel_array[i][122:388, 112:403])
        elif j == 2:
            сur_img = np.copy(pixel_array[i][120:390, 98:420])
        elif j == 3:
            сur_img = np.copy(pixel_array[i][120:390, 115:405])
        else:
            raise ValueError

        #resize back to (512, 512)
        сur_img = cv2.resize(сur_img, (512, 512), interpolation=cv2.INTER_CUBIC)
        new_pixel_array[i] = np.copy(сur_img)

    # сохраняем npz в новом формате
    new_series_path = series_path.replace(
        "phillips_dataset", "nonframe_phillips_dataset"
    )
    np.savez(
        new_series_path,
        pixel_array=new_pixel_array,
        pixel_path=hearth_npz_info["pixel_path"],
        artery_type=hearth_npz_info["artery_type"],
        seriesid=hearth_npz_info["seriesid"],
        studyid=hearth_npz_info["studyid"],
        series_number=hearth_npz_info["series_number"],
        primary_angle=hearth_npz_info["primary_angle"],
        secondary_angle=hearth_npz_info["secondary_angle"],
        is_collaterals=hearth_npz_info["is_collaterals"],
        is_occlusion=hearth_npz_info["is_occlusion"],
        is_undefined_type=hearth_npz_info["is_undefined_type"],
        is_artifact=hearth_npz_info["is_artifact"],
    )
