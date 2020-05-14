# 2020/5/4
import os
import json
import tarfile
import zipfile
import shutil
import requests
from bs4 import BeautifulSoup
import re
from pymatgen import Element
from pathlib import Path
import subprocess


def predownload():
    # SSSP precision wget
    sssp_cutoff_table_url = "https://www.materialscloud.org/discover/data/discover/sssp/downloads/sssp_precision.json"
    with open("sssp_precision.json", "wb") as f:
        shutil.copyfileobj(requests.get(
            sssp_cutoff_table_url, stream=True).raw, f)
    sssp_pbe_url = "https://www.materialscloud.org/discover/data/discover/sssp/downloads/SSSP_precision_pseudos.tar.gz"
    with open("SSSP_precision_pseudos.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(sssp_pbe_url, stream=True).raw, f)
    with tarfile.open('SSSP_precision_pseudos.tar.gz', 'r') as tf:
        tf.extractall(path='SSSP_PBE_precision_pseudos')
    sssp_pbesol_url = "https://www.materialscloud.org/discover/data/discover/sssp/downloads/tarball_PBEsol_1.1/SSSP_PBEsol_pseudos.tar.gz"
    with open("SSSP_PBEsol_pseudos.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(sssp_pbesol_url, stream=True).raw, f)
    with tarfile.open('SSSP_PBEsol_pseudos.tar.gz', 'r') as tf:
        tf.extractall(path='SSSP_PBEsol_pseudos')
    shutil.move(
        "SSSP_PBEsol_pseudos/SSSP_PBEsol_pseudos/SSSP_PBEsol_precision_pseudos", "./SSSP_PBEsol_pseudos")
    shutil.rmtree("SSSP_PBEsol_pseudos")
    # Dojo wget
    dojo_cutoff_table_url = "https://raw.githubusercontent.com/abinit/pseudo_dojo/master/pseudo_dojo/pseudos/ONCVPSP-PBE-PDv0.4/stringent.djson"
    print(requests.get(dojo_cutoff_table_url, stream=True).text,
          file=open("stringent.djson", "w"))
    dojo_pbe_sr_url = "http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_stringent_upf.tgz"
    with open("nc-sr-04_pbe_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(dojo_pbe_sr_url, stream=True).raw, f)
    with tarfile.open('nc-sr-04_pbe_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-sr-04_pbe_stringent_upf')
    dojo_pbe_fr_url = "http://www.pseudo-dojo.org/pseudos/nc-fr-04_pbe_stringent_upf.tgz"
    with open("nc-fr-04_pbe_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(dojo_pbe_fr_url, stream=True).raw, f)
    with tarfile.open('nc-fr-04_pbe_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-fr-04_pbe_stringent_upf')
    dojo_pbesol_sr_url = "http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbesol_stringent_upf.tgz"
    with open("nc-sr-04_pbesol_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(
            dojo_pbesol_sr_url, stream=True).raw, f)
    with tarfile.open('nc-sr-04_pbesol_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-sr-04_pbesol_stringent_upf')
    dojo_pbesol_fr_url = "http://www.pseudo-dojo.org/pseudos/nc-fr-04_pbesol_stringent_upf.tgz"
    with open("nc-fr-04_pbesol_stringent_upf.tgz", "wb") as f:
        shutil.copyfileobj(requests.get(
            dojo_pbesol_fr_url, stream=True).raw, f)
    with tarfile.open('nc-fr-04_pbesol_stringent_upf.tgz', 'r') as tf:
        tf.extractall(path='nc-fr-04_pbesol_stringent_upf')
    # SG15 wget
    sg15_url = "http://www.quantum-simulation.org/potentials/sg15_oncv/sg15_oncv_upf_2020-02-06.tar.gz"
    with open("sg15_oncv_upf_2020-02-06.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(sg15_url, stream=True).raw, f)
    with tarfile.open('sg15_oncv_upf_2020-02-06.tar.gz', 'r') as tf:
        tf.extractall(path='sg15_oncv_upf_2020-02-06')
    sg15_url = "https://github.com/pipidog/ONCVPSP/archive/master.zip"
    with open("master.zip", "wb") as f:
        shutil.copyfileobj(requests.get(sg15_url, stream=True).raw, f)
    with zipfile.ZipFile('master.zip', 'r') as z:
        z.extractall('master')
    shutil.move("master/ONCVPSP-master/sg15", ".")
    shutil.rmtree("master")
    # GBRV wget
    gbrv_pbe_url = "https://www.physics.rutgers.edu/gbrv/all_pbe_UPF_v1.5.tar.gz"
    with open("all_pbe_UPF_v1.5.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(gbrv_pbe_url, stream=True).raw, f)
    with tarfile.open('all_pbe_UPF_v1.5.tar.gz', 'r') as tf:
        tf.extractall(path='all_pbe_UPF_v1.5')
    gbrv_pbesol_url = "https://www.physics.rutgers.edu/gbrv/all_pbesol_UPF_v1.5.tar.gz"
    with open("all_pbesol_UPF_v1.5.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(gbrv_pbesol_url, stream=True).raw, f)
    with tarfile.open('all_pbesol_UPF_v1.5.tar.gz', 'r') as tf:
        tf.extractall(path='all_pbesol_UPF_v1.5')
    # pslib 100 031 compile
    ps100_url = "https://people.sissa.it/dalcorso/pslibrary/pslibrary.1.0.0.tar.gz"
    with open("pslibrary.1.0.0.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(ps100_url, stream=True).raw, f)
    with tarfile.open('pslibrary.1.0.0.tar.gz', 'r') as tf:
        tf.extractall(path='pslibrary.1.0.0')
    ps031_url = "https://people.sissa.it/dalcorso/pslibrary/pslibrary.0.3.1.tar.gz"
    with open("pslibrary.0.3.1.tar.gz", "wb") as f:
        shutil.copyfileobj(requests.get(ps031_url, stream=True).raw, f)
    with tarfile.open('pslibrary.0.3.1.tar.gz', 'r') as tf:
        tf.extractall(path='pslibrary.0.3.1')
    with open("./paths.json", "r") as f:
        QEroot = json.load(f)["QE"][:-4]
    QE_path = ["#!/bin/bash\n", QEroot]
    print(*QE_path, sep="\n", file=open("pslibrary.1.0.0/QE_path", "w"))
    print(*QE_path, sep="\n", file=open("pslibrary.0.3.1/QE_path", "w"))
    os.mkdir("./pslib100")
    os.mkdir("./pslib031")
    target100 = ["./pslibrary.1.0.0/pbe", "./pslibrary.1.0.0/pbesol",
                 "./pslibrary.1.0.0/rel-pbe", "./pslibrary.1.0.0/rel-pbesol"]
    target031 = ["./pslibrary.0.3.1/pbe", "./pslibrary.0.3.1/pbesol",
                 "./pslibrary.0.3.1/rel-pbe", "./pslibrary.0.3.1/rel-pbesol"]
    for t100 in target100:
        subprocess.run("./make_ps", cwd=t100)
        shutil.move(f"{t100}/*", "./pslib100/")
    for t031 in target031:
        subprocess.run("./make_ps", cwd=t031)
        shutil.move(f"{t031}/*", "./pslib031/")
    # pslib dict
    pslib_list_html = requests.get(
        "https://dalcorso.github.io/pslibrary/PP_list.html")
    pslib_list_soup = BeautifulSoup(pslib_list_html.text, "html.parser")
    pslib_list_text = pslib_list_soup.find(class_="highlight").text
    pslib_list_tmp = [i.strip() for i in pslib_list_text.split("\n")]
    pslib_list = [re.sub(':.*_psl.', ':', i)
                  for i in pslib_list_tmp][:-1]  # the last is empty
    pslib_ver_dict = {}
    [pslib_ver_dict.update({i.split(":")[0]: "100" if (i.split(":")[1] in ["1.0.0", "1.0.1"]) else "031"})
     for i in pslib_list]
    pslib_ver_dict["W"] = "100"
    # https://www.materialscloud.org/discover/sssp/plot/precision/In
    pslib_ver_dict["In"] = "031"
    # download CoInMn from QE home page and swap 031
    os.mkdir("./pslib_CoInMn")
    pslib_url = "https://www.quantum-espresso.org/pseudopotentials/ps-library/"
    for special in ["co", "in", "mn"]:
        html = requests.get(pslib_url + special.lower())
        soup = BeautifulSoup(html.text, "html.parser")
        candidates = [i.getText()
                      for i in soup.select("#content-right td")]
        for func in ["PBE", "PBESOL"]:
            for pstype in ["PAW", "USPP"]:
                for soc in ["Scalar relativistic", "Full relativistic"]:
                    ps_file_indexs = [j for j, k in enumerate(
                        candidates) if f"Pseudopotential type: {pstype} \nFunctional type: {func}\nNon Linear Core Correction\n{soc}\n" in k]
                    if len(ps_file_indexs) < 1:
                        print(
                            f"I can't find a {func}-{soc}-{pstype} potential of {special} in PS Library")
                    else:
                        if len(ps_file_indexs) > 1:
                            # select older version
                            ps_file_index = sorted(
                                [(len(candidates[j]), j) for j in ps_file_indexs], reverse=True)[0][1]
                        else:
                            ps_file_index = ps_file_indexs[0]
                    href = "https://www.quantum-espresso.org" + \
                        soup.select(
                            "#content-right td")[ps_file_index].find("a")["href"]
                    ps_text = requests.get(href).text
                    ps_file = soup.select(
                        "#content-right td")[ps_file_index].find("a").getText().strip()
                    print(ps_text, sep='\n', file=open(
                        f"./pslib_CoInMn/{ps_file}", 'w'))
    shutil.copy("./pslib_CoInMn/*", "./pslib031/")


def pslib_read(psfile):
    # UPF ver2
    flag = True
    start = False
    wfc = []
    ecutwfc = None
    with open(psfile, "r") as f:
        while flag:
            line = f.readline()
            if ("cutoff" in line) and ("wavefunctions" in line):
                tmp = (re.search(r"\d+", line)).group()
                if tmp:
                    ecutwfc = int(tmp) if int(tmp) > 0 else None
            if "Valence configuration" in line:
                start = True
            if start:
                if line.split()[0][0] in ["1", "2", "3", "4", "5", "6", "7"]:
                    wfc.append([int(line.split()[0][0]), int(line.split()[2])])
                elif "Generation configuration" in line:
                    flag = False
    return wfc, ecutwfc


def gbrv_read(psfile):
    # UPF ver1
    flag = True
    start = False
    end = False
    wfc = []
    with open(psfile, "r") as f:
        while flag:
            line = f.readline()
            if "<PP_HEADER>" in line:
                start = True
            if start:
                if "Wavefunctions" in line.split()[0]:
                    end = True
            if end:
                if line.split()[0][0] in ["1", "2", "3", "4", "5", "6", "7"]:
                    wfc.append([int(line.split()[0][0]), int(line.split()[1])])
                elif "</PP_HEADER>" in line:
                    flag = False
    return wfc


def oncv_read(psfile, element):
    flag = True
    start = False
    end = False
    tmp = []
    with open(psfile, "r") as f:
        while flag:
            line = f.readline()
            if "PP_INPUTFILE" in line:
                start = True
            if start:
                if element in line.split():
                    nv = int(line.split()[3])
                elif "energy" in line.split():
                    start = False
                    end = True
            if end:
                if line.split()[0] in ["1", "2", "3", "4", "5", "6", "7"]:
                    tmp.append([int(line.split()[0]), int(line.split()[1])])
                elif "PSEUDOPOTENTIAL AND OPTIMIZATION" in line:
                    flag = False
    wfc = tmp[-1 * nv:]
    return wfc


if __name__ == "__main__":
    # pre definition
    predownload()  # download pseudo potentials
    os.makedir("pseudos")
    pstypes = {'100PAW': "PAW", 'GBRV-1.5': "US", 'SG15': "ONCV", 'GBRV-1.4': "US", '031US': "US",
               '100US': "US", '031PAW': "PAW", 'SG15-1.1': "ONCV", 'GBRV-1.2': "US", 'Dojo': "ONCV", 'Wentzcovitch': "PAW"}
    Hubbards = {'Sc': {'U': 2.9}, 'Ti': {'U': 4.4}, 'V': {'U': 2.7}, 'Cr': {'U': 3.5}, 'Mn': {'U': 4.0}, 'Fe': {'U': 4.6}, 'Co': {'U': 5.0}, 'Ni': {'U': 5.1}, 'Cu': {'U': 4.0}, 'Zn': {'U': 7.5}, 'Ga': {'U': 3.9}, 'Sn': {'U': 3.5}, 'Nb': {'U': 2.1}, 'Mo': {'U': 2.4}, 'Ta': {'U': 2.0}, 'W': {'U': 2.2}, 'Tc': {'U': 2.7}, 'Ru': {'U': 3.0}, 'Rh': {'U': 3.3}, 'Pd': {'U': 3.6}, 'Ag': {'U': 5.8}, 'Cd': {'U': 2.1}, 'In': {
        'U': 1.9}, 'Re': {'U': 2.4}, 'Os': {'U': 2.6}, 'Ir': {'U': 2.8}, 'Pt': {'U': 3.0}, 'Au': {'U': 4.0}, 'La': {'U': 8.1, 'J': 0.6}, 'Ce': {'U': 7.0, 'J': 0.7}, 'Pr': {'U': 6.5, 'J': 1.0}, 'Nd': {'U': 7.2, 'J': 1.0}, 'Sm': {'U': 7.4, 'J': 1.0}, 'Eu': {'U': 6.4, 'J': 1.0}, 'Gd': {'U': 6.7, 'J': 0.1}, 'Dy': {'U': 5.6, 'J': 0.0}, 'Tm': {'U': 7.0, 'J': 1.0}, 'Yb': {'U': 7.0, 'J': 0.67}, 'Lu': {'U': 4.8, 'J': 0.95}}
    with open("./sssp_precision.json", "r") as f:
        ssspp = json.load(f)
    with open("./stringent.djson", "r") as f:
        dojo = json.load(f)
    with open("./pslib_ver.json", "r") as f:
        pslib_ver_dict = json.load(f)
    # choose pseudo potentials
    elements = {}
    # for element in ["Ag"]:
    for element in ssspp:
        SOC = "sr"
        if (31 <= Element(element).Z <= 34) | (39 <= Element(element).Z <= 52) | (57 <= Element(element).Z <= 84):
            SOC = "fr"
        default = {"pstype": pstypes[ssspp[element]["pseudopotential"]],
                   "SOC": SOC, "Hubbard": Hubbards.get(element)}
        pseudopotential = {"PBE": {"sr": {}, "fr": {}},
                           "PBEsol": {"sr": {}, "fr": {}}}
        # Lanthanoid
        print(element)
        if 57 <= Element(element).Z <= 71:
            for func in ["PBE", "PBEsol"]:
                for pawus in [("PAW", "kjpaw"), ("US", "rrkjus")]:  # no ONCV
                    for rel in [("sr", ""), ("fr", "rel-")]:
                        candidates = [p.name for p in Path(
                            "./pslib" + pslib_ver_dict[element]).glob(f"{element}.{rel[1]}{func.lower()}-*{pawus[1]}_*")]
                        # take longer name which consider more semi-core states
                        psfile = max(candidates, key=len)
                        shutil.copy(
                            f"./pslib100/{psfile}", "./pseudos/")
                        wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
                        pseudopotential[func][rel[0]][pawus[0]] = {
                            "filename": psfile, "cutoff": ecutwfc, "dual": 8.0, "wfc": wfc, "cite": f"100{pawus[0]}", "resource": "PSLibrary100"}
        else:
            # PAW
            for func in ["PBE", "PBEsol"]:
                # sr
                if default["pstype"] == "PAW":
                    psfile = [s for s in os.listdir(
                        f"./SSSP_{func}_pseudos") if re.split("[._-]", s)[0].lower() == f"{element}".lower()][0]
                    shutil.copy(
                        f"./SSSP_{func}_pseudos/{psfile}", "./pseudos/")
                    wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
                    ecutwfc = ssspp[element]["cutoff"]
                    cite = ssspp[element]["pseudopotential"]
                    resource = "SSSP"
                else:
                    candidates = [p.name for p in Path(
                        "./pslib" + pslib_ver_dict[element]).glob(f"{element}.{func.lower()}-*kjpaw_*")]
                    # take longer name which consider more semi-core states
                    psfile = max(candidates, key=len)
                    shutil.copy(
                        f"./pslib{pslib_ver_dict[element]}/{psfile}", "./pseudos/")
                    wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
                    cite = pslib_ver_dict[element] + "PAW"
                    resource = "PSLibrary" + \
                        pslib_ver_dict[element] if (
                            element not in ["Co", "In", "Mn"]) else "Readey-to-use"
                dual = 8.0 if not (element in ["Mn", "Fe", "Co"]) else 12.0
                pseudopotential[func]["sr"]["PAW"] = {
                    "filename": psfile, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
                # fr
                psfile_fr = psfile.replace(
                    f"{element}.pbe", f"{element}.rel-pbe")
                shutil.copy(
                    list(Path("./").glob(f"pslib[0-9][0-9][0-9]/{psfile_fr[:-4]}*"))[0], "./pseudos/")
                pseudopotential[func]["fr"]["PAW"] = {
                    "filename": psfile_fr, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
            # US
            for func in ["PBE", "PBEsol"]:
                # sr
                if default["pstype"] == "US":
                    psfile = [s for s in os.listdir(
                        f"./SSSP_{func}_pseudos") if re.split("[._-]", s)[0].lower() == f"{element}".lower()][0]
                    shutil.copy(
                        f"./SSSP_{func}_pseudos/{psfile}", "./pseudos/")
                    ecutwfc = ssspp[element]["cutoff"]
                    cite = ssspp[element]["pseudopotential"]
                    resource = "SSSP"
                    if "GBRV" in cite:
                        wfc = gbrv_read(f"./pseudos/{psfile}")
                    else:
                        wfc = pslib_read(f"./pseudos/{psfile}")[0]
                # Po, At, Noble gas
                elif Element(element).Z in [2, 10, 18, 36, 54, 84, 85, 86]:
                    candidates = [p.name for p in Path(
                        "./pslib" + pslib_ver_dict[element]).glob(f"{element}.{func.lower()}-*rrkjus_*")]
                    # take longer name which consider more semi-core states
                    psfile = max(candidates, key=len)
                    shutil.copy(
                        f"./pslib{pslib_ver_dict[element]}/{psfile}", "./pseudos/")
                    wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
                    cite = pslib_ver_dict[element] + "US"
                    resource = "PSLibrary" + pslib_ver_dict[element]
                else:
                    candidates = [p.name for p in Path(
                        f"all_{func.lower()}_UPF_v1.5").glob(f"{element.lower()}_*")]
                    psfile = candidates[0]  # only one
                    shutil.copy(
                        f"all_{func.lower()}_UPF_v1.5/{psfile}", "./pseudos/")
                    ecutwfc = ssspp[element]["cutoff"] if default["pstype"] == "PAW" else max(
                        ssspp[element]["cutoff"] / 2, 40)
                    cite = "GBRV-1.5"
                    resource = "GBRV1.5"
                    wfc = gbrv_read(f"./pseudos/{psfile}")
                dual = 8.0 if not (element in ["Mn", "Fe", "Co"]) else 12.0
                pseudopotential[func]["sr"]["US"] = {
                    "filename": psfile, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
                # fr
                candidates = [p.name for p in Path(
                    "./pslib" + pslib_ver_dict[element]).glob(f"{element}.{func.lower()}-*rrkjus_*")]
                # take longer name which consider more semi-core states
                psfile = max(candidates, key=len)
                shutil.copy(
                    f"./pslib{pslib_ver_dict[element]}/{psfile}", "./pseudos/")
                wfc, ecutwfc = pslib_read(f"./pseudos/{psfile}")
                cite = pslib_ver_dict[element] + "US"
                resource = "PSLibrary" + \
                    pslib_ver_dict[element] if (
                        element not in ["Co", "In", "Mn"]) else "Readey-to-use"
                pseudopotential[func]["fr"]["US"] = {
                    "filename": psfile, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
            # ONCV
            # sr PBE
            dual = 4.0
            if default["pstype"] == "ONCV":
                cite = ssspp[element]["pseudopotential"]
                if cite == "Dojo":
                    psfile = f"{element}_ONCV_PBE-0.4.dojo.upf"  # rename
                    shutil.copy(
                        f"./SSSP_PBE_pseudos/{ssspp[element]['filename']}", f"./pseudos/{psfile}")
                else:
                    psfile = ssspp[element]["filename"]
                    shutil.copy(f"./SSSP_PBE_pseudos/{psfile}", "./pseudos/")
                wfc = oncv_read(f"./pseudos/{psfile}", element)
                ecutwfc = ssspp[element]["cutoff"]
                resource = "SSSP"
            elif element == "Rn":  # Dojo
                psfile = "Rn_ONCV_PBE_fr-0.4.dojo.upf"
                shutil.copy("./nc-fr-04_pbe_stringent_upf/Rn.upf",
                            f"./pseudos/{psfile}")
                wfc = oncv_read(f"./pseudos/{psfile}", "Rn")
                ecutwfc = dojo["pseudos_metadata"]["Rn"]["hints"]["high"]["ecut"]
                cite = "Dojo"
                resource = "PseudoDojo0.4"
            else:
                candidates = [p for p in Path("./").glob(f"./sg15/{element}_ONCV_PBE_sr.upf")] + [
                    p for p in Path("./").glob(f"./sg15_oncv_upf_2020-02-06/{element}_ONCV_PBE-*")]
                psfile = f"{element}_ONCV_PBE-1.0.oncvpsp.upf"  # rename
                # take newer
                shutil.copy(f"{max(candidates)}",
                            f"./pseudos/{psfile}")
                wfc = oncv_read(f"./pseudos/{psfile}", element)
                ecutwfc = dojo["pseudos_metadata"][element]["hints"]["high"]["ecut"]
                cite = "SG15"
                resource = "ONCVPSP"
            pseudopotential["PBE"]["sr"]["ONCV"] = {
                "filename": psfile, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
            # fr PBE
            if cite == "Dojo":
                psfile_fr = f"{element}_ONCV_PBE_fr-0.4.dojo.upf"
                shutil.copy(f"./nc-fr-04_pbe_stringent_upf/{element}.upf",
                            f"./pseudos/{psfile_fr}")
            else:
                candidates = [p for p in Path("./").glob(f"./sg15/{element}_ONCV_PBE_fr.upf")] + [
                    p for p in Path("./").glob(f"./sg15_oncv_upf_2020-02-06/{element}_ONCV_PBE_FR-*")]
                psfile_fr = f"{element}_ONCV_PBE_fr-1.0.oncvpsp.upf"
                shutil.copy(f"{max(candidates)}",
                            f"./pseudos/{psfile_fr}")
            pseudopotential["PBE"]["fr"]["ONCV"] = {
                "filename": psfile_fr, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
            # sr PBEsol
            if default["pstype"] == "ONCV":
                psfile = [s for s in os.listdir(
                    f"./SSSP_PBEsol_pseudos") if re.split("[._-]", s)[0].lower() == f"{element}".lower()][0]
                shutil.copy(f"./SSSP_PBEsol_pseudos/{psfile}", "./pseudos/")
                wfc = oncv_read(f"./pseudos/{psfile}", element)
                ecutwfc = ssspp[element]["cutoff"]
                cite = ssspp[element]["pseudopotential"]
                resource = "SSSP"
            else:
                psfile = f"{element}_ONCV_PBEsol-0.4.dojo.upf"  # rename
                shutil.copy(
                    f"./nc-sr-04_pbesol_stringent_upf/{element}.upf", f"./pseudos/{psfile}")
                wfc = oncv_read(f"./pseudos/{psfile}", element)
                ecutwfc = dojo["pseudos_metadata"][element]["hints"]["high"]["ecut"]
                cite = "Dojo"
                resource = "PseudoDojo0.4"
            pseudopotential["PBEsol"]["sr"]["ONCV"] = {
                "filename": psfile, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
            # fr PBEsol
            psfile = f"{element}_ONCV_PBEsol_fr-0.4.dojo.upf"
            shutil.copy(
                f"./nc-fr-04_pbesol_stringent_upf/{element}.upf", f"./pseudos/{psfile}")
            wfc = oncv_read(f"./pseudos/{psfile}", element)
            ecutwfc = dojo["pseudos_metadata"][element]["hints"]["high"]["ecut"]
            cite = "Dojo"
            resource = "PseudoDojo0.4"
            pseudopotential["PBEsol"]["fr"]["ONCV"] = {
                "filename": psfile, "cutoff": ecutwfc, "dual": dual, "wfc": wfc, "cite": cite, "resource": resource}
        # save to elements
        elements.update(
            {element: {"pseudopotential": pseudopotential, "default": default}})
    with open("./elements.json", "w") as f:
        json.dump(elements, f)
    # remove folder and move downloads to misc
    shutil.rmtree("SSSP_PBE_pseudos")
    shutil.rmtree("SSSP_PBEsol_pseudos")
    shutil.rmtree("all_pbe_UPF_v1.5")
    shutil.rmtree("all_pbesol_UPF_v1.5")
    shutil.rmtree("nc-fr-04_pbe_stringent_upf")
    shutil.rmtree("nc-fr-04_pbesol_stringent_upf")
    shutil.rmtree("nc-sr-04_pbe_stringent_upf")
    shutil.rmtree("nc-sr-04_pbesol_stringent_upf")
    shutil.rmtree("pslib031")
    shutil.rmtree("pslib100")
    shutil.rmtree("pslib_CoInMn")
    shutil.rmtree("sg15")
    shutil.rmtree("sg15_oncv_upf_2020-02-06")
    shutil.move("*.tar.gz", "misc/")
    shutil.move("*.tgz", "misc/")
    shutil.move("*.zip", "misc/")
    shutil("pslib_ver.json", "misc/")
    shutil("sssp_precision.json", "misc/")
    shutil("stringent.djson", "misc/")
