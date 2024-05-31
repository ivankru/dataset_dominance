#Creation of the main dataset for Philips Allura Clarity
#DICOM convertion into .npz files
#additional tags extraction from .xlsx file

import os
import shutil

import numpy as np
import pandas as pd
import pydicom as dcm
from tqdm import tqdm


def create_dir_for_column(
    column_name, df, cur_dirname, path_new_data, sheet_name=None, column_type=None
):
    """
    Function to create a folder for a specific tag "normal", "occlusion", "bad_quality", "artifact" or "uncertainty"
    column_name - path
    df - Dataframe with additional information about patients
    cur_dirname - название папки с текущей особенностью
    path_new_data - the path to save angiographic studies
    sheet_name - sheet name in .xlsx file
    column_type - tag e.g. "occlusion", "atifacts"...
    """

    if not os.path.exists(os.path.join(path_new_data, cur_dirname)):
        os.mkdir(os.path.join(path_new_data, cur_dirname))

    #select tags from .xlsx
    if column_type == "True/nan":
        df_for_new_dir = df[df[column_name] == 1].copy()
    elif column_type == "strings":
        df_for_new_dir = df[df[column_name].notna()].copy()
    else:
        df_for_new_dir = df.copy()

    #values to save in .npz
    if not df_for_new_dir.empty:
        for idx in df_for_new_dir.index:

            cur_filename = df_for_new_dir.loc[idx, "путь"]
            cur_filename = cur_filename[cur_filename.find("Файл") :]
            cur_path = df_for_new_dir.loc[idx, "путь"].replace("home", "/home")
            new_filename = df_for_new_dir.loc[idx, "StudyInstanceUID"]

            if sheet_name == "ЛЕВЫЕ ТИПЫ 112":
                cur_path = cur_path.replace("Файл", "Файл")

            # проверка корректности таблицы
            assert (
                "левый" in df_for_new_dir.loc[idx, "Истинный тип"]
                or "правый" in df_for_new_dir.loc[idx, "Истинный тип"]
            ), f"Некорректные поле Истинный лист {df_for_new_dir.loc[idx, 'Истинный тип']}"

            dominance_dir = (
                "Left_Dominance"
                if "левый" in df_for_new_dir.loc[idx, "Истинный тип"]
                else "Right_Dominance"
            )

            path_to_files = os.path.join(path_new_data, cur_dirname, dominance_dir)
            #cur_path = cur_path.replace("", "")

            if not os.path.exists(path_to_files):
                os.mkdir(path_to_files)
            shutil.copytree(cur_path, os.path.join(path_to_files, cur_filename))

            # new patient folder name
            new_filename = cur_filename + "_" + new_filename

            os.rename(
                os.path.join(path_to_files, cur_filename),
                os.path.join(path_to_files, new_filename),
            )

            # new name for an artery folder
            artery_type = rename_to_lca_rca(path_to_files, new_filename, cur_path)

            # save .npz file
            path_to_studyid = os.path.join(path_to_files, new_filename)
            save_npz_from_dcm(path_to_studyid, cur_path, df_for_new_dir, artery_type)
            check_correctness_lca_rca(path_to_files, new_filename, cur_path)

        print(f"{column_name}: {len(df_for_new_dir)}")

        df = df.drop(index=df_for_new_dir.index)
    return df


def rename_to_lca_rca(path_to_files, new_filename, cur_path):
    """
    to replace folders rca, lca to RCA and LCA
    path_to_files - path to studies
    new_filename - new file name
    """

    if len(os.listdir(os.path.join(path_to_files, new_filename))) == 1:
        raise ValueError(
            f"Нет RCA или LCA в папке {os.path.join(path_to_files, new_filename)}"
        )

    for left_right_dirname in os.listdir(os.path.join(path_to_files, new_filename)):

        # Унифицирование названий папок для артериий и удаление лишних файлов
        if "лев" in left_right_dirname or left_right_dirname.lower() == "лка":
            os.rename(
                os.path.join(path_to_files, new_filename, left_right_dirname),
                os.path.join(path_to_files, new_filename, "LCA"),
            )
            artery_type = "LCA"

        elif "прав" in left_right_dirname or left_right_dirname.lower() == "пка":
            os.rename(
                os.path.join(path_to_files, new_filename, left_right_dirname),
                os.path.join(path_to_files, new_filename, "RCA"),
            )
            artery_type = "RCA"

        elif left_right_dirname[0] == ".":
            os.remove(os.path.join(path_to_files, new_filename, left_right_dirname))

        elif left_right_dirname.startswith("Left_Coro"):
            raise ValueError(
                f"Left coro встетилось в папке {os.path.join(path_to_files, new_filename)}"
            )

        else:
            print(cur_path)
            print(
                os.path.join(path_new_data, cur_dirname, new_filename),
                left_right_dirname,
            )
            raise ValueError

    return artery_type


def save_npz_from_dcm(path_to_studyid, cur_path, df_for_new_dir, artery_type):
    """
    From .dcm to .npz
    path_to_studyid - path to the study in a new dataset
    cur_path - current path
    df_for_new_dir - table with additional tags
    """

    for lca_rca_dirname in os.listdir(path_to_studyid):
        # all angiographic views
        for left_coro_dirname in os.listdir(
            os.path.join(path_to_studyid, lca_rca_dirname)
        ):
            projection_path = os.path.join(
                path_to_studyid, lca_rca_dirname, left_coro_dirname
            )
            file_flag = False
            not_one_file_flag = False

            if not os.path.isdir(projection_path):
                file_flag = True

            elif len(os.listdir(projection_path)) > 1:

                lengths_files = []

                for file in os.listdir(projection_path):
                    assert file.endswith(".dcm"), f"Неправильный файл в папке, {file}"
                    dcm_file = os.path.join(projection_path, file)
                    ds = dcm.dcmread(dcm_file)
                    pixel_array = ds.pixel_array
                    if len(pixel_array.shape) == 2:
                        len_pixel_arr = 1
                    else:
                        len_pixel_arr = len(pixel_array)

                lengths_files.append([len_pixel_arr, dcm_file])
                lengths_files.sort(key=lambda x: x[0])

                target_dcm_file = lengths_files[-1][1]
                not_one_file_flag = True

            # path to dcm
            if file_flag:
                dcm_file = projection_path

            elif not_one_file_flag:
                dcm_file = target_dcm_file

            else:
                dcm_file = os.listdir(projection_path)[0]
                dcm_file = os.path.join(projection_path, dcm_file)

            # extract matadata from .dcm
            ds = dcm.dcmread(dcm_file)
            SeriesInstanceUID = ds["SeriesInstanceUID"].value
            StudyInstanceUID = ds["StudyInstanceUID"].value

            try:
                PositionerPrimaryAngle = ds["PositionerPrimaryAngle"].value
            except:  # noqa: E722
                PositionerPrimaryAngle = None
                print("PositionerPrimaryAngle is None")

            try:
                PositionerSecondaryAngle = ds["PositionerSecondaryAngle"].value
            except:  # noqa: E722
                PositionerSecondaryAngle = None
                print("PositionerSecondaryAngle is None")

            # to avoid duplicates
            if SeriesInstanceUID in SeriesInstanceUID_arr:
                print("Встретилась та же серия")
                return

            SeriesInstanceUID_arr.add(SeriesInstanceUID)

            # video extcration
            pixel_for_npz = ds.pixel_array
            if sheet_name == "ЛЕВЫЕ ТИПЫ 112":
                filename_in_table = (
                    cur_path[cur_path.find("Файл") :]
                    .replace(" ", "")
                    .replace("Файл", "Файл")
                )
            else:
                filename_in_table = cur_path[cur_path.find("Файл") :].replace(" ", "")
            dcm_file_oldfilename = [os.path.basename(dcm_file)]

            for path_list in os.walk(cur_path):
                if dcm_file_oldfilename in path_list:
                    break

            real_path_to_dcm = []
            for part_of_path_list in path_list:
                if part_of_path_list == []:
                    continue
                if isinstance(part_of_path_list, list):
                    real_path_to_dcm.append(part_of_path_list[0])
                else:
                    real_path_to_dcm.append(part_of_path_list)

            if real_path_to_dcm is []:
                raise AssertionError("Проблемы с поиском dcm файла в исходном каталоге")
            real_path_to_dcm = os.path.join(*real_path_to_dcm)

            path_for_npz = real_path_to_dcm.replace("..", "home")
            is_collaterals_for_npz = df_for_new_dir[
                df_for_new_dir["имя папки"] == filename_in_table
            ]["коллатерали"]
            is_collaterals_for_npz = (
                True if is_collaterals_for_npz.tolist()[0] == 1 else False
            )

            if len(pixel_for_npz.shape) == 2:
                delete_left_coro_data(file_flag, projection_path)
                continue

            save_npz_path = os.path.join(
                path_to_studyid, lca_rca_dirname, SeriesInstanceUID
            )

            assert artery_type in ["LCA", "RCA"]

            if file_flag:
                series_number = "1"
            else:
                series_number = dcm_file.split("/")[-2].split("_")[-1]
            np.savez(
                save_npz_path,
                pixel_array=pixel_for_npz,
                pixel_path=path_for_npz,
                seriesid=SeriesInstanceUID,
                studyid=StudyInstanceUID,
                is_collaterals=is_collaterals_for_npz,
                series_number=int(series_number),
                primary_angle=PositionerPrimaryAngle,
                secondary_angle=PositionerSecondaryAngle,
                artery_type=artery_type,
            )

            delete_left_coro_data(file_flag, projection_path)

        assert len(os.listdir(path_to_studyid)) > 0, f"Пустая папка, {path_to_studyid}"


def delete_left_coro_data(file_flag, projection_path):
    """
    Delete a file or a folder
    file_flag - is file
    projection_path - path to
    """

    if file_flag:
        os.remove(projection_path)
    else:
        shutil.rmtree(projection_path)


def check_correctness_lca_rca(path_to_files, new_filename, cur_path):
    """
    Check correctness of LCA and RCA
    file_flag - if is a file
    projection_path - angiograpic view path
    path_to_files - path to studies
    new_filename - folder name
    """

    if sorted(os.listdir(os.path.join(path_to_files, new_filename))) != sorted(
        ["RCA", "LCA"]
    ):
        print(os.path.join(path_to_files, new_filename))
        raise AssertionError("Файл не соответствует структуре [RCA, LCA]")
    for left_right_dirname in os.listdir(os.path.join(path_to_files, new_filename)):
        if (
            len(
                os.listdir(
                    os.path.join(path_to_files, new_filename, left_right_dirname)
                )
            )
            < 1
        ):
            raise AssertionError(
                f"В папке присутствует пустая папка {left_right_dirname}"
            )


def save_npz_from_dcm_new_left_dominance(
    projection_path, filename, artery_type, new_filename_path
):
    """
    To add additional Left dominant studies
    path_to_studyid - path to study_id
    cur_path - path in the old dataset
    df_for_new_dir - data frame with additional tags
    """

    not_one_file_flag = False
    file_flag = False

    # to find out dublicates
    global study_set
    global filename_set

    if not os.path.isdir(projection_path):
        file_flag = True

    #if there are more than one .dcm file
    elif len(os.listdir(projection_path)) > 1:

        lengths_files = []

        for file in os.listdir(projection_path):
            assert file.endswith(".dcm"), f"Неправильный файл в папке, {file}"
            dcm_file = os.path.join(projection_path, file)
            ds = dcm.dcmread(dcm_file)
            if "Image Storage" not in ds.SOPClassUID.name:
                print("ERROR in dcm storage")
                continue
            pixel_array = ds.pixel_array
            if len(pixel_array.shape) == 2:
                len_pixel_arr = 1
            else:
                len_pixel_arr = len(pixel_array)

        #if more then one file, choose with the largest number of frames
        lengths_files.append([len_pixel_arr, dcm_file])
        lengths_files.sort(key=lambda x: x[0])

        target_dcm_file = lengths_files[-1][1]
        not_one_file_flag = True

    # path to dcm
    if file_flag:
        dcm_file = projection_path

    elif not_one_file_flag:
        dcm_file = target_dcm_file

    else:
        dcm_file = os.listdir(projection_path)[0]
        dcm_file = os.path.join(projection_path, dcm_file)

    ds = dcm.dcmread(dcm_file)
    if "Image Storage" not in ds.SOPClassUID.name:
        print("2 ERROR in dcm storage")
        return

    SeriesInstanceUID = ds["SeriesInstanceUID"].value
    StudyInstanceUID = ds["StudyInstanceUID"].value
    pixel_for_npz_arr = ds.pixel_array

    try:
        PositionerPrimaryAngle = ds["PositionerPrimaryAngle"].value
    except:  # noqa: E722
        PositionerPrimaryAngle = None
        print("PositionerPrimaryAngle is None")

    try:
        PositionerSecondaryAngle = ds["PositionerSecondaryAngle"].value
    except:  # noqa: E722
        PositionerSecondaryAngle = None
        print("PositionerSecondaryAngle is None")

    # Проверка, что данная серия не встречалась ранее
    if SeriesInstanceUID in SeriesInstanceUID_arr:
        print("Встретилась та же серия")
        return

    SeriesInstanceUID_arr.add(SeriesInstanceUID)

    if StudyInstanceUID in study_set and filename not in filename_set:
        print(StudyInstanceUID)
    else:
        study_set.add(StudyInstanceUID)
        filename_set.add(filename)

    if len(pixel_for_npz_arr.shape) == 2:
        return

    # save npz in a new folder
    new_filename_path = new_filename_path + "_" + StudyInstanceUID
    new_artery_path = os.path.join(new_filename_path, artery_type)
    assert os.path.basename(new_artery_path) in [
        "LCA",
        "RCA",
    ], f"{os.path.basename(new_artery_path)}"
    if file_flag:
        left_coro_idx = "1"
    else:
        left_coro_idx = os.path.basename(projection_path).split("_")[-1]
    npz_path = os.path.join(new_artery_path, left_coro_idx + "_" + SeriesInstanceUID)

    for dir_to_create in [new_filename_path, new_artery_path]:
        if not os.path.exists(dir_to_create):
            os.mkdir(dir_to_create)

    if "occlusion" in npz_path:
        is_collaterals = True
    else:
        is_collaterals = False

    assert artery_type in ["LCA", "RCA"]

    np.savez(
        npz_path,
        pixel_array=pixel_for_npz_arr,
        pixel_path=dcm_file,
        series_number=int(left_coro_idx),
        seriesid=SeriesInstanceUID,
        studyid=StudyInstanceUID,
        is_collaterals=is_collaterals,
        primary_angle=PositionerPrimaryAngle,
        secondary_angle=PositionerSecondaryAngle,
        artery_type=artery_type,
    )


def save_npz_from_dcm_new_occlusion(
    projection_path, dominant_type, filename, artery_type, new_dataset_path
):
    """
    Occlusion dataset
    path_to_studyid - path to study_id
    cur_path - path in the old dataset
    df_for_new_dir - data frame with additional tags    """

    not_one_file_flag = False

    # to avoid dublicates
    global study_set
    global filename_set

    # if a file or dir
    if not os.path.isdir(projection_path):
        print("Dir is a file")
        return

    if dominant_type == "left":
        new_dominance_path = os.path.join(new_dataset_path, "Left_Dominance")
    elif dominant_type == "right":
        new_dominance_path = os.path.join(new_dataset_path, "Right_Dominance")

    #if more then one .dcm file
    if len(os.listdir(projection_path)) > 1:

        lengths_files = []

        for file in os.listdir(projection_path):
            assert file.endswith(".dcm"), f"Неправильный файл в папке, {file}"
            dcm_file = os.path.join(projection_path, file)
            ds = dcm.dcmread(dcm_file)
            pixel_array = ds.pixel_array
            if len(pixel_array.shape) == 2:
                len_pixel_arr = 1
            else:
                len_pixel_arr = len(pixel_array)

        #choose a serie with the largest number of frames
        lengths_files.append([len_pixel_arr, dcm_file])
        lengths_files.sort(key=lambda x: x[0])

        target_dcm_file = lengths_files[-1][1]
        not_one_file_flag = True

    #dcm path
    if not_one_file_flag:
        dcm_file = target_dcm_file

    else:
        dcm_file = os.listdir(projection_path)[0]
        dcm_file = os.path.join(projection_path, dcm_file)

    # metadata from .dcm
    ds = dcm.dcmread(dcm_file)
    SeriesInstanceUID = ds["SeriesInstanceUID"].value
    StudyInstanceUID = ds["StudyInstanceUID"].value

    try:
        PositionerPrimaryAngle = ds["PositionerPrimaryAngle"].value
    except:  # noqa: E722
        PositionerPrimaryAngle = None
        print("PositionerPrimaryAngle is None")

    try:
        PositionerSecondaryAngle = ds["PositionerSecondaryAngle"].value
    except:  # noqa: E722
        PositionerSecondaryAngle = None
        print("PositionerSecondaryAngle is None")

    if SeriesInstanceUID in SeriesInstanceUID_arr:
        print("Встретилась та же серия")
        return

    #add new serie ID to the list
    SeriesInstanceUID_arr.add(SeriesInstanceUID)

    if StudyInstanceUID in study_set and filename not in filename_set:
        print(StudyInstanceUID)
    else:
        study_set.add(StudyInstanceUID)
        filename_set.add(filename)

    pixel_for_npz_arr = ds.pixel_array

    if len(pixel_for_npz_arr.shape) == 2:
        return

    # save npz
    new_filename_path = os.path.join(
        new_dominance_path, filename + "_" + StudyInstanceUID
    )
    new_artery_path = os.path.join(new_filename_path, artery_type)
    npz_path = os.path.join(new_artery_path, SeriesInstanceUID)

    left_coro_idx = os.path.basename(projection_path).split("_")[-1]

    for dir_to_create in [new_dominance_path, new_filename_path, new_artery_path]:
        if not os.path.exists(dir_to_create):
            os.mkdir(dir_to_create)

    assert artery_type in ["LCA", "RCA"]

    np.savez(
        npz_path,
        pixel_array=pixel_for_npz_arr,
        pixel_path=dcm_file,
        seriesid=SeriesInstanceUID,
        studyid=StudyInstanceUID,
        series_number=int(left_coro_idx),
        is_collaterals=True,
        primary_angle=PositionerPrimaryAngle,
        secondary_angle=PositionerSecondaryAngle,
        artery_type=artery_type,
    )


def create_occlusion_dataset(dominant_type, new_dataset_path):
    """
    to add new portion of occlusions.
    dominant_type - left / right dominance
    new_dataset_path - new data set path
    """

    print(dominant_type)

    if dominant_type == "left":
        dataset_path = "/home/mazanov/data/datasets/occlusions/Окклюзии (5 отправка)/Окклюзии левый тип_аноним"
    elif dominant_type == "right":
        dataset_path = "/home/mazanov/data/datasets/occlusions/Окклюзии (5 отправка)/Окклюзии правый тип_ аноним"
    assert os.path.exists(dataset_path), dataset_path

    for filename in tqdm(os.listdir(dataset_path)):
        filename_path = os.path.join(dataset_path, filename)
        for lca_rca_dir in os.listdir(filename_path):
            lca_rca_dir_path = os.path.join(filename_path, lca_rca_dir)
            if "ЛКА" == lca_rca_dir:
                artery_type = "LCA"
            elif "ПКА" == lca_rca_dir:
                artery_type = "RCA"
            elif lca_rca_dir.startswith("."):
                continue
            else:
                raise AssertionError(f"Wrong artery type: {lca_rca_dir}")

            for left_coro_dir in os.listdir(lca_rca_dir_path):
                left_coro_dir_path = os.path.join(lca_rca_dir_path, left_coro_dir)

                # создание и сохранение npz в новую папку
                save_npz_from_dcm_new_occlusion(
                    left_coro_dir_path,
                    dominant_type,
                    filename,
                    artery_type,
                    new_dataset_path,
                )


#to avoid dublicates
SeriesInstanceUID_arr = set()
StudyInstanceUID_arr = set()

#pathes to the data
path_table = "*****.xlsx"
path_data = "*****"
path_new_data = "*****"
data_dirs = [
    "Разметка КАГ по типам/левый тип (50)",
    "Разметка КАГ по типам/правый тип (230)",
    "Разметка КАГ по типам 2/Левый тип",
    "Разметка КАГ по типам 2/Правый тип",
    "ЛЕВЫЕ ТИПЫ 112/No_Name/Левый тип",
]

#new folder creation
if not os.path.exists(path_new_data):
    os.mkdir(path_new_data)
else:
    shutil.rmtree(path_new_data)
    os.mkdir(path_new_data)

#additional tags are in .xlsx file
xls = pd.ExcelFile(path_table)

# iterate through all sheets
for i, (sheet_name, dirname) in enumerate(zip(xls.sheet_names[:-2], data_dirs)):
    if True:

        print("\n", sheet_name)
        df = pd.read_excel(xls, sheet_name)

        #poor quality
        column_name = "плохое качество"
        cur_dirname = "bad_quality"
        df = create_dir_for_column(
            column_name,
            df,
            cur_dirname,
            path_new_data,
            sheet_name=sheet_name,
            column_type="True/nan",
        )

        #challenging case
        column_name = "сложный случай"
        cur_dirname = "challenging_case"
        df = create_dir_for_column(
            column_name,
            df,
            cur_dirname,
            path_new_data,
            sheet_name=sheet_name,
            column_type="strings",
        )

        # aftifacts
        column_name = "артефакт"
        cur_dirname = "artifact"
        df = create_dir_for_column(
            column_name,
            df,
            cur_dirname,
            path_new_data,
            sheet_name=sheet_name,
            column_type="strings",
        )

        # uncertainty
        column_name = "неопределенный тип"
        cur_dirname = "indefinite_type"
        df = create_dir_for_column(
            column_name,
            df,
            cur_dirname,
            path_new_data,
            sheet_name=sheet_name,
            column_type="True/nan",
        )

        # occlusions
        column_name = "окклюзия"
        cur_dirname = "occlusion"
        df = create_dir_for_column(
            column_name,
            df,
            cur_dirname,
            path_new_data,
            sheet_name=sheet_name,
            column_type="True/nan",
        )

        # normals
        column_name = "остаток"
        cur_dirname = "normal"
        df_for_new_dir = df.copy()
        df = create_dir_for_column(
            column_name, df, cur_dirname, path_new_data, sheet_name=sheet_name
        )

#to avoid dublicates
study_set = set()
filename_set = set()

new_dataset_path = "/home/mazanov/data/datasets/main_dataset/occlusion"

#new occlusions
create_occlusion_dataset("left", new_dataset_path)
create_occlusion_dataset("right", new_dataset_path)

#new left dominant studies
path_new_left_dominance = (
    "***"
)

new_left_dominance_occlusion = [
    "Cardiac - R201901161051224",
    "Cardiac - R201908290957213",
]

new_left_dominance_artifact = [
    "Cardiac - R201901251109171",
    "Cardiac - R201904041012016",
    "Cardiac - R201904111134257",
    "Cardiac - R201904121226381",
    "Cardiac - R201905141301110",
    "Cardiac - R201905141310219",
]

new_left_dominance_bad_quality = [
    "Cardiac - R201902131020100",
    "Cardiac - R201904031413456",
    "Cardiac - R201902071135146",
    "Cardiac - R201902251148585",
    "Cardiac - R201906211220249",
    "Cardiac - R201907091202272",
    "Cardiac - R201908231700096",
]

new_left_dominance_indefinite_type = ["Cardiac - R201907261156449"]

#seek an appropriate folder for .
for data_filename in tqdm(os.listdir(path_new_left_dominance)):

    filename_path = os.path.join(path_new_left_dominance, data_filename)
    if data_filename in new_left_dominance_occlusion:
        new_dataset_path = os.path.join(path_new_data, "occlusion", "Left_Dominance")
    elif data_filename in new_left_dominance_artifact:
        new_dataset_path = os.path.join(path_new_data, "artifact", "Left_Dominance")
    elif data_filename in new_left_dominance_bad_quality:
        new_dataset_path = os.path.join(path_new_data, "bad_quality", "Left_Dominance")
    elif data_filename in new_left_dominance_indefinite_type:
        new_dataset_path = os.path.join(
            path_new_data, "indefinite_type", "Left_Dominance"
        )
    else:
        new_dataset_path = os.path.join(path_new_data, "normal", "Left_Dominance")

    new_filename_path = os.path.join(
        new_dataset_path, "".join(data_filename.split(" - "))
    )

    for artery_type in os.listdir(filename_path):
        if artery_type.upper() not in ["ЛКА", "ПКА"]:
            print(f"Wrong dirname: {artery_type}, file {data_filename}")
            break

        artery_path = os.path.join(filename_path, artery_type)
        new_artery_type = "LCA" if artery_type.upper() == "ЛКА" else "RCA"
        new_artery_path = os.path.join(new_filename_path, new_artery_type)

        for left_coro_dirname in os.listdir(artery_path):
            assert left_coro_dirname.startswith(
                "Left_Coro_"
            ) or left_coro_dirname.endswith(".dcm"), left_coro_dirname
            left_coro_path = os.path.join(artery_path, left_coro_dirname)

            # создание и сохранение npz
            save_npz_from_dcm_new_left_dominance(
                left_coro_path, filename_path, new_artery_type, new_filename_path
            )

# The total number of studies
print("\n")
count_all_files = 0
for dirname in os.listdir(path_new_data):
    count_files_in_dir = 0
    if os.path.exists(os.path.join(path_new_data, dirname, "Left_Dominance")):
        count_files_in_dir += len(
            os.listdir(os.path.join(path_new_data, dirname, "Left_Dominance"))
        )
    if os.path.exists(os.path.join(path_new_data, dirname, "Right_Dominance")):
        count_files_in_dir += len(
            os.listdir(os.path.join(path_new_data, dirname, "Right_Dominance"))
        )
    count_all_files += count_files_in_dir
    print(f"Всего {dirname}: {count_files_in_dir}")

# Подсчет пациентов с левыми и правыми доминантностями в сете
left_count = 0
right_count = 0
for feature in os.listdir(path_new_data):
    feature_path = os.path.join(path_new_data, feature)
    for dominance in os.listdir(feature_path):
        dominance_path = os.path.join(feature_path, dominance)
        if "Left" in dominance:
            left_count += len(os.listdir(dominance_path))
        elif "Right" in dominance:
            right_count += len(os.listdir(dominance_path))
        else:
            raise AssertionError(dominance)
print("Всего Left dominance:", left_count)
print("Всего Right dominance:", right_count)

print("Всего файлов:", count_all_files)
print("Всего файлов:", count_all_files)
