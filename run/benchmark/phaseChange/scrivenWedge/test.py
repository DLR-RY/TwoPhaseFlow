from typing import Dict, List, Union
from pathlib import Path
import os

# excludeFolders = ['postProcessing',"Cases"]


def is_float(f):
    try:
        float(f)
    except:
        return False
    return True


def exclude_folder(dir, exclude_folders):
    if dir in exclude_folders:
        return True
    elif "processor" in dir:
        return True
    elif is_float(dir):
        return True
    return False


def exclude_file(f, exclude_files):
    if "log" in f:
        return True
    if any(f in exFile for exFile in exclude_files):
        return True
    return False


def _gen_file_filter(files, cat):
    if type(files) is str:
        file_gen = (
            p.relative_to(files)
            for p in Path(files).rglob("*")
            if str(p).endswith(tuple(cat))
        )
    elif type(files) is list:
        file_gen = (p for p in files if p.endswith(tuple(cat)))
    else:
        ValueError("only string and list of string are supported")
    return file_gen


def _categorize_by_dict(classify_dict: Dict[str, List[str]], files):
    """categorizes file based on the dictionary

    classify_dict = {
                    "data": [".csv"],
                    "media": [".png",".mp4"],
                    "doc": [".pdf"]
                }
    Args:
        d ([type]): [description]u
        files ([type]): [description]

    Returns:
        [type]: [description]
    """

    cat_files = []
    for key in classify_dict:
        # cat_files[key] = []
        file_gen = _gen_file_filter(files, classify_dict[key])
        for f in file_gen:
            cat_files.append({"category": key, "path": f})
            # cat_files[key].append(str(f))
    return cat_files


classify_dict = {"data": [".csv"], "media": [".png", ".mp4"], "doc": [".pdf"]}


def categorize_files(
    case_dir: Union[str, Path],
    exclude_folder,
    exclude_file,
    exclude_folders=["postProcessing", "Cases"],
    exclude_files=["log"],
):
    files_with_category = []
    for root, dirs, files in os.walk(case_dir):
        dirs[:] = [d for d in dirs if not exclude_folder(d, exclude_folders)]
        files[:] = [f for f in files if not exclude_file(f, exclude_files)]
        filter_files = [str(Path(root, f)) for f in files]
        cat_files = _categorize_by_dict(classify_dict, filter_files)
        if cat_files:
            files_with_category.append(cat_files)

    return files_with_category


path = "."



print(categorize_files(path, exclude_folder, exclude_file))

